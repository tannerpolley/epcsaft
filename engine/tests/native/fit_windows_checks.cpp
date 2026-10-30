// Exercise the production facilities without adding a numerical evaluator or runtime test switches.
#define NOMINMAX
#define _WIN32_WINNT 0x0A00
#include <windows.h>
#include <sddl.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
struct FitInterrupted {};
struct { std::vector<std::string> worker; } fit_host;
bool cancel = false;
bool interrupted() { return cancel; }
#include "regression/fit_windows.hpp"

void check(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
void image_locks() {
    wchar_t directory[MAX_PATH];
    const DWORD n = GetTempPathW(MAX_PATH, directory);
    windows_require(n > 0 && n < MAX_PATH);
    for (int i = 0; i < 2; ++i) {
        struct TemporaryImage {
            wchar_t path[MAX_PATH]{};
            std::wstring renamed;
            HMODULE module = nullptr;
            ~TemporaryImage() {
                if (module) FreeLibrary(module);
                DeleteFileW(path);
                DeleteFileW(renamed.c_str());
            }
        } temporary;
        windows_require(GetTempFileNameW(directory, L"fit", 0, temporary.path) != 0);
        temporary.renamed = std::wstring(temporary.path) + L".renamed";
        windows_require(CopyFileW(native_image.path.c_str(), temporary.path, FALSE));
        temporary.module = LoadLibraryExW(temporary.path, nullptr, 0);
        windows_require(temporary.module != nullptr);
        auto image = open_image(temporary.module);
        // Remove the loader's mapping: only the production file owner may deny these operations.
        windows_require(FreeLibrary(temporary.module));
        temporary.module = nullptr;
        check(!DeleteFileW(temporary.path), "production image lock allowed delete");
        check(GetLastError() == ERROR_SHARING_VIOLATION, "delete failed for another reason");
        check(!MoveFileW(temporary.path, temporary.renamed.c_str()), "production image lock allowed rename");
        check(GetLastError() == ERROR_SHARING_VIOLATION, "rename failed for another reason");
        image.file = Handle{};
        windows_require(MoveFileW(temporary.path, temporary.renamed.c_str()));
        windows_require(DeleteFileW(temporary.renamed.c_str()));
        std::cout << "temporary image " << i << ": delete/rename refused; both succeed after owner release\n";
    }
}
DWORD handles() {
    DWORD count = 0;
    windows_require(GetProcessHandleCount(GetCurrentProcess(), &count));
    return count;
}
void marker(Workers& workers) {
    char byte = 0;
    transfer(workers.sockets.front(), &byte, 1, false);
    check(byte == 'R', "child did not reach inheritance check");
}
}

int main(int argc, char** argv) {
    try {
        if (argc > 1) {
            // The numeric value is intentionally absent from the explicit inheritance list.
            SetEvent(reinterpret_cast<HANDLE>(std::stoull(argv[2])));
            char byte = 'R';
            transfer({reinterpret_cast<HANDLE>(std::stoull(argv[3]))}, &byte, 1, true);
            Sleep(INFINITE);
            return 0;
        }
        image_locks();
        char path[32768];
        windows_require(GetModuleFileNameA(nullptr, path, sizeof path) > 0);
        extension_image = open_image(reinterpret_cast<const void*>(&main));
        SECURITY_ATTRIBUTES security{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
        Handle unrelated(CreateEventW(&security, TRUE, FALSE, nullptr));
        windows_require(unrelated.handle != nullptr);
        fit_host.worker = {path, "hold", std::to_string(reinterpret_cast<std::uintptr_t>(unrelated.handle))};
        // Initialize Windows process-launch facilities before counting per-fit resources.
        { Workers warm; warm.start(); marker(warm); }
        const auto before = handles();
        for (bool out : {false, true}) {
            {
                Workers workers;
                workers.start();
                marker(workers);
                check(WaitForSingleObject(unrelated.handle, 0) == WAIT_TIMEOUT, "unrelated event inherited");
                std::vector<char> bytes(out ? 1024 * 1024 : 1);
                cancel = true;
                bool caught = false;
                try { transfer(workers.sockets.front(), bytes.data(), bytes.size(), out); }
                catch (const FitInterrupted&) { caught = true; }
                cancel = false;
                check(caught, "transfer did not cancel during pending I/O");
                // Callback is reached only after a 100 ms pending-operation timeout.
                std::cout << (out ? "write" : "read") << ": pending cancellation drained\n";
            }
            check(handles() == before, "pending cancellation leaked owned handles");
        }
        {
            Workers surviving;
            surviving.start();
            marker(surviving);
            {
                Workers cancelled;
                cancelled.start();
                marker(cancelled);
            }
            check(WaitForSingleObject(surviving.processes.front().handle, 0) == WAIT_TIMEOUT,
                  "one fit cleanup terminated the concurrent fit");
            windows_require(TerminateProcess(surviving.processes.front().handle, 17));
            windows_require(WaitForSingleObject(surviving.processes.front().handle, 30000) == WAIT_OBJECT_0);
            char byte = 0;
            bool refused = false;
            try {
                transfer(surviving.sockets.front(), &byte, 1, false);
            } catch (const std::runtime_error& error) {
                refused = true;
                check(std::string(error.what()).find("997") == std::string::npos, "stale pending error");
                std::cout << "dead peer: " << error.what() << '\n';
            }
            check(refused, "dead peer accepted read");
        }
        check(handles() == before, "worker-set cleanup leaked handles");
        check(WaitForSingleObject(unrelated.handle, 0) == WAIT_TIMEOUT, "unrelated handle changed");
        std::cout << "explicit inheritance and concurrent worker-set isolation passed; owned handle delta 0\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
