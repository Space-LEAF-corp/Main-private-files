import os
import tempfile
import unittest

from jarvondis.jarvondis import Jarvondis, Personality


class TestJarvondis(unittest.TestCase):
    def test_learn_and_persistence_csv(self):
        fd, path = tempfile.mkstemp(suffix=".csv")
        os.close(fd)
        try:
            j = Jarvondis(personality=Personality(tone="witty"), memory_file=path, memory_format="csv")
            out = j.respond("hello world")
            self.assertIn("Echoing back", out)
            j.save_memory()

            # load into a fresh instance
            j2 = Jarvondis(personality=Personality(tone="witty"), memory_file=path, memory_format="csv")
            mem = j2.memory
            self.assertTrue(len(mem) >= 1)
            last = mem[-1]
            self.assertEqual(last.get("input"), "hello world")
        finally:
            try:
                os.remove(path)
            except Exception:
                pass


if __name__ == "__main__":
    unittest.main()
