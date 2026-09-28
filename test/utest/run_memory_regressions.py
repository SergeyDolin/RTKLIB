#!/usr/bin/env python3
"""Run production-code memory regressions out of tree, optionally with sanitizers."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ('rtkcmn rinex rtkpos postpos solution lambda geoid sbas preceph pntpos '
           'ephemeris options ppp ppp_ar rtcm rtcm2 rtcm3 rtcm3e ionex tides ambupd gpt3').split()
TESTS = ['t_excsats']

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sanitize', action='store_true')
    parser.add_argument('--build-dir', type=Path)
    args = parser.parse_args()
    out = (args.build_dir or Path(tempfile.mkdtemp(prefix='rtk-memory-'))).resolve()
    out.mkdir(parents=True, exist_ok=True)
    flags = ['-O1', '-g', '-ansi', '-ffunction-sections', '-fdata-sections',
             '-fno-omit-frame-pointer', '-I'+str(ROOT/'src'), '-DTRACE',
             '-DENAGLO', '-DENAQZS', '-DENAGAL', '-DENACMP', '-DENAIRN',
             '-DNFREQ=5', '-DNEXOBS=2']
    if args.sanitize:
        flags += ['-fsanitize=address,undefined']
    for test in TESTS:
        command = [os.environ.get('CC', 'gcc')]+flags
        command += [str(ROOT/'test/utest'/(test+'.c'))]
        command += [str(ROOT/'src'/(name+'.c')) for name in SOURCES]
        command += ['-Wl,--gc-sections', '-lm', '-o', str(out/test)]
        with (out/(test+'.build.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        env = os.environ.copy()
        if args.sanitize:
            env['ASAN_OPTIONS'] = 'detect_leaks=1:halt_on_error=1'
            env['UBSAN_OPTIONS'] = 'print_stacktrace=1:halt_on_error=1'
        with (out/(test+'.stdout')).open('w') as stdout, (out/(test+'.stderr')).open('w') as stderr:
            result = subprocess.run([str(out/test)], cwd=out, env=env,
                                    stdout=stdout, stderr=stderr)
        print((out/(test+'.stdout')).read_text(), end='', flush=True)
        if result.returncode:
            print((out/(test+'.stderr')).read_text(), flush=True)
            raise SystemExit(result.returncode)
        print(test, 'PASS', flush=True)
    print('Logs:', out)

if __name__ == '__main__':
    main()
