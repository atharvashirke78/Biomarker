"""Every script in nextflow/bin with a shebang must be executable: processes call them by name, and a non-executable one fails with
'bad interpreter: Permission denied' (internal_validation.R did, silently, because its shell fallback hides the error; the stub run never runs the script)."""
import os
import stat
import unittest

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'nextflow', 'bin')


class T(unittest.TestCase):
    def test_scripts_with_shebang_are_executable(self):
        bad = []
        for n in sorted(os.listdir(BIN)):
            p = os.path.join(BIN, n)
            if os.path.isfile(p) and not n.endswith(('.bak', '.md')) and '.bak' not in n:
                with open(p, 'rb') as f:
                    if f.read(2) == b'#!' and not os.stat(p).st_mode & stat.S_IXUSR:
                        bad.append(n)
        self.assertEqual(bad, [], 'not executable: %s' % bad)


if __name__ == '__main__':
    unittest.main()
