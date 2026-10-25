"""Guarded native fit worker: python -I fit_launch.py <image.json> <log.jsonl> -I -c <unchanged native command>.

The native fit host passes its extension and native images as descriptors 3/4. This adds the
campaign audit guard and descriptor hash checks, records the worker, then runs the native
command unchanged (it loads the images directly and never imports the editable package).
"""
import json
import os
import sys

sys.addaudithook(lambda event, args: (_ for _ in ()).throw(RuntimeError(f'campaign guard: {event} refused'))
                 if event in ('subprocess.Popen', 'os.system', 'os.posix_spawn', 'os.spawn', 'os.exec') else None)
image_path, log_path, flag, c, command = sys.argv[1:6]
assert (flag, c) == ('-I', '-c') and sys.flags.isolated, sys.argv
import hashlib
image = json.loads(open(image_path).read())
record = dict(pid=os.getpid(), ppid=os.getppid(), executable=sys.executable)
for fd, name in ((3, image['extension']), (4, 'libepcsaft_native.so')):
    with open(f'/proc/self/fd/{fd}', 'rb') as stream:
        found = hashlib.file_digest(stream, 'sha256').hexdigest()
    if found != image['libraries_sha256'][name]:
        os._exit(3)
    record[f'fd{fd}'] = dict(path=os.readlink(f'/proc/self/fd/{fd}'), sha256=found)
with open(log_path, 'a') as log:
    log.write(json.dumps(record) + '\n')
exec(compile(command, '<native fit worker>', 'exec'), {'__name__': '__main__'})
