// Included inside the fit implementation's anonymous namespace; Windows process facilities only.
#undef near // Windows' legacy qualifier must not erase shared C++ identifiers.
struct Channel { HANDLE pipe; HANDLE process = nullptr; };
struct Handle final {
    HANDLE handle = nullptr;
    Handle() = default;
    explicit Handle(HANDLE h) : handle(h == INVALID_HANDLE_VALUE ? nullptr : h) {}
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
    Handle(Handle&& other) noexcept : handle(std::exchange(other.handle, nullptr)) {}
    Handle& operator=(Handle&& other) noexcept {
        if (handle) CloseHandle(handle);
        handle = std::exchange(other.handle, nullptr);
        return *this;
    }
    ~Handle() { if (handle) CloseHandle(handle); }
};
void windows_require(bool ok) {
    if (!ok) throw std::runtime_error("fit worker Windows error " + std::to_string(GetLastError()));
}
std::wstring wide(const std::string& value) {
    const int n = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, value.data(), static_cast<int>(value.size()), nullptr, 0);
    windows_require(n > 0);
    std::wstring out(n, L'\0');
    windows_require(MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, value.data(), static_cast<int>(value.size()), out.data(), n) == n);
    return out;
}
struct Image final {
    Handle file;
    std::wstring path;
    explicit operator bool() const { return file.handle != nullptr; }
};
Image open_image(const void* address) {
    HMODULE module = nullptr;
    windows_require(GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
        reinterpret_cast<LPCWSTR>(address), &module));
    std::wstring path(32768, L'\0');
    const DWORD n = GetModuleFileNameW(module, path.data(), static_cast<DWORD>(path.size()));
    windows_require(n > 0 && n < path.size());
    path.resize(n);
    // Deny write/delete sharing for as long as a worker can load the image.
    Handle file(CreateFileW(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr));
    windows_require(file.handle != nullptr);
    return {std::move(file), std::move(path)};
}
Image extension_image;
const Image native_image = open_image(reinterpret_cast<const void*>(&interrupted));

// One pending operation at a time; drain cancellation before its stack storage disappears.
void transfer(Channel fd, char* at, std::size_t left, bool out) {
    Handle event(CreateEventW(nullptr, TRUE, FALSE, nullptr));
    windows_require(event.handle != nullptr);
    while (left > 0) {
        OVERLAPPED operation{};
        operation.hEvent = event.handle;
        windows_require(ResetEvent(event.handle));
        DWORD moved = 0;
        const DWORD bytes = static_cast<DWORD>(std::min<std::size_t>(left, MAXDWORD));
        const BOOL ok = out ? WriteFile(fd.pipe, at, bytes, &moved, &operation) : ReadFile(fd.pipe, at, bytes, &moved, &operation);
        if (!ok) {
            windows_require(GetLastError() == ERROR_IO_PENDING);
            try {
                for (;;) {
                    const HANDLE waiting[] = {event.handle, fd.process};
                    const DWORD wait = WaitForMultipleObjects(fd.process ? 2 : 1, waiting, FALSE, 100);
                    if (wait == WAIT_OBJECT_0) break;
                    if (wait == WAIT_OBJECT_0 + 1) throw std::runtime_error("fit worker exited during pipe transfer");
                    windows_require(wait != WAIT_FAILED);
                    if (wait != WAIT_TIMEOUT) throw std::runtime_error("fit worker unexpected pipe wait result");
                    if (interrupted()) throw FitInterrupted{};
                }
                windows_require(GetOverlappedResult(fd.pipe, &operation, &moved, FALSE));
            } catch (...) {
                CancelIoEx(fd.pipe, &operation);
                GetOverlappedResult(fd.pipe, &operation, &moved, TRUE);
                throw;
            }
        }
        if (moved == 0) throw std::runtime_error("fit worker pipe closed");
        at += moved;
        left -= moved;
    }
}
// CommandLineToArgvW/CRT quoting, including backslashes immediately before a closing quote.
std::wstring quote_argument(const std::wstring& word) {
    std::wstring out = L"\"";
    std::size_t slashes = 0;
    for (const auto c : word) {
        if (c == L'\\') { ++slashes; continue; }
        out.append(slashes * (c == L'"' ? 2 : 1), L'\\');
        slashes = 0;
        if (c == L'"') out += L'\\';
        out += c;
    }
    out.append(2 * slashes, L'\\');
    return out + L'"';
}
struct Workers final {
    Handle job;
    std::vector<Channel> sockets;
    std::vector<Handle> processes;
    ~Workers() {
        if (job.handle) TerminateJobObject(job.handle, 1);
        for (auto fd : sockets) CloseHandle(fd.pipe);
        for (auto& process : processes) WaitForSingleObject(process.handle, INFINITE);
    }
    void start() {
        if (!job.handle) {
            job = Handle(CreateJobObjectW(nullptr, nullptr));
            windows_require(job.handle != nullptr);
            JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
            limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
            windows_require(SetInformationJobObject(job.handle, JobObjectExtendedLimitInformation, &limits, sizeof limits));
        }
        // A current-user-only, local pipe. Inherit the connected endpoint; never a discoverable channel name.
        PSECURITY_DESCRIPTOR descriptor = nullptr;
        windows_require(ConvertStringSecurityDescriptorToSecurityDescriptorW(L"D:P(A;;GA;;;OW)", SDDL_REVISION_1, &descriptor, nullptr));
        SECURITY_ATTRIBUTES security{sizeof(SECURITY_ATTRIBUTES), descriptor, FALSE};
        static std::atomic<unsigned long> sequence{0};
        const auto name = L"\\\\.\\pipe\\epcsaft-fit-" + std::to_wstring(GetCurrentProcessId()) + L"-" + std::to_wstring(sequence++);
        Handle server(CreateNamedPipeW(name.c_str(), PIPE_ACCESS_DUPLEX | FILE_FLAG_OVERLAPPED | FILE_FLAG_FIRST_PIPE_INSTANCE,
            PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT | PIPE_REJECT_REMOTE_CLIENTS, 1, 65536, 65536, 0, &security));
        LocalFree(descriptor);
        windows_require(server.handle != nullptr);
        Handle client(CreateFileW(name.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr, OPEN_EXISTING, FILE_FLAG_OVERLAPPED, nullptr));
        windows_require(client.handle != nullptr);
        // CreateFile connected the local client before this call.
        Handle connected(CreateEventW(nullptr, TRUE, FALSE, nullptr));
        windows_require(connected.handle != nullptr);
        OVERLAPPED connection{};
        connection.hEvent = connected.handle;
        windows_require(ConnectNamedPipe(server.handle, &connection) || GetLastError() == ERROR_PIPE_CONNECTED);
        Handle null(CreateFileW(L"NUL", GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
            nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr));
        windows_require(null.handle != nullptr);
        Handle inherited[6];
        HANDLE list[6];
        const HANDLE sources[] = {client.handle, extension_image.file.handle, native_image.file.handle,
            null.handle, GetStdHandle(STD_OUTPUT_HANDLE), GetStdHandle(STD_ERROR_HANDLE)};
        for (int i = 0; i < 6; ++i) {
            const HANDLE source = sources[i] && sources[i] != INVALID_HANDLE_VALUE ? sources[i] : null.handle;
            windows_require(DuplicateHandle(GetCurrentProcess(), source, GetCurrentProcess(), &list[i], 0, TRUE, DUPLICATE_SAME_ACCESS));
            inherited[i] = Handle(list[i]);
        }
        std::wstring command;
        for (const auto& word : fit_host.worker) command += quote_argument(wide(word)) + L' ';
        for (int i = 0; i < 3; ++i) command += quote_argument(std::to_wstring(reinterpret_cast<std::uintptr_t>(list[i]))) + L' ';
        command += quote_argument(extension_image.path) + L' ' + quote_argument(native_image.path);
        SIZE_T bytes = 0;
        InitializeProcThreadAttributeList(nullptr, 2, 0, &bytes);
        std::vector<unsigned char> storage(bytes);
        auto* attributes = reinterpret_cast<LPPROC_THREAD_ATTRIBUTE_LIST>(storage.data());
        windows_require(InitializeProcThreadAttributeList(attributes, 2, 0, &bytes));
        struct Attributes { LPPROC_THREAD_ATTRIBUTE_LIST p; ~Attributes() { DeleteProcThreadAttributeList(p); } } cleanup{attributes};
        windows_require(UpdateProcThreadAttribute(attributes, 0, PROC_THREAD_ATTRIBUTE_HANDLE_LIST, list, sizeof list, nullptr, nullptr));
        windows_require(UpdateProcThreadAttribute(attributes, 0, PROC_THREAD_ATTRIBUTE_JOB_LIST, &job.handle, sizeof(HANDLE), nullptr, nullptr));
        STARTUPINFOEXW startup{};
        startup.StartupInfo.cb = sizeof startup;
        startup.StartupInfo.dwFlags = STARTF_USESTDHANDLES;
        startup.StartupInfo.hStdInput = list[3];
        startup.StartupInfo.hStdOutput = list[4];
        startup.StartupInfo.hStdError = list[5];
        startup.lpAttributeList = attributes;
        PROCESS_INFORMATION process{};
        // Allocate tracking before creation: an allocation failure must not leave an untracked process.
        processes.reserve(processes.size() + 1);
        sockets.reserve(sockets.size() + 1);
        windows_require(CreateProcessW(wide(fit_host.worker.front()).c_str(), command.data(), nullptr, nullptr, TRUE,
            EXTENDED_STARTUPINFO_PRESENT | CREATE_NO_WINDOW, nullptr, nullptr, &startup.StartupInfo, &process));
        CloseHandle(process.hThread);
        processes.emplace_back(process.hProcess);
        sockets.push_back({std::exchange(server.handle, nullptr), process.hProcess});
    }
};
