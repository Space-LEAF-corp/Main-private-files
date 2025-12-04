"""Tests for CodeValidator and FibonacciResequencer."""
import unittest

from fibonacci_resequencer import CodeValidator, FibonacciResequencer, fib_reverse_sequence


class TestFibonacciSequence(unittest.TestCase):
    def test_fib_reverse_sequence(self):
        seq = fib_reverse_sequence()
        self.assertEqual(seq, [21, 13, 8, 5, 3, 2, 1])


class TestCodeValidator(unittest.TestCase):
    def test_valid_python_forward(self):
        code = "x = 1\ny = 2\nprint(x + y)"
        validator = CodeValidator(code)
        ok, msg = validator.validate_forward()
        self.assertTrue(ok)

    def test_invalid_python_forward(self):
        code = "x = 1\nif x:\n  (missing close)"
        validator = CodeValidator(code)
        ok, msg = validator.validate_forward()
        self.assertFalse(ok)

    def test_validate_both(self):
        code = "def foo():\n    return 42"
        validator = CodeValidator(code)
        ok, msg = validator.validate_both()
        self.assertTrue(ok)


class TestFibonacciResequencer(unittest.TestCase):
    def setUp(self):
        # Create a code sample with enough lines to demonstrate chunking
        lines = [f"# Line {i}" for i in range(50)]
        self.code = '\n'.join(lines)
        self.resequencer = FibonacciResequencer(self.code)

    def test_split_into_chunks(self):
        chunks = self.resequencer.split_into_chunks()
        chunk_sizes = [len(c) for c in chunks]
        # Should have chunks of sizes ~21, 13, 8, 5, 3, 2, 1 + remainder
        self.assertGreater(len(chunks), 0)
        self.assertEqual(chunk_sizes[0], 21)  # first chunk should be 21
        if len(chunks) > 1:
            self.assertEqual(chunk_sizes[1], 13)  # second should be 13

    def test_resequence_reverse(self):
        chunks = self.resequencer.split_into_chunks()
        reversed_code = self.resequencer.resequence_reverse(chunks)
        reversed_lines = reversed_code.split('\n')
        # Reversed should start with the last chunk's first line
        self.assertEqual(len(reversed_lines), len(self.code.split('\n')))

    def test_split_and_resequence(self):
        resequenced, chunk_info = self.resequencer.split_and_resequence(reverse=True)
        self.assertIsNotNone(resequenced)
        self.assertGreater(len(chunk_info), 0)

    def test_analyze_chunks(self):
        analysis = self.resequencer.analyze_chunks()
        self.assertIn("num_chunks", analysis)
        self.assertIn("fib_sizes", analysis)
        self.assertIn("chunk_sizes", analysis)
        self.assertEqual(analysis["fib_sizes"], [21, 13, 8, 5, 3, 2, 1])


if __name__ == "__main__":
    unittest.main()
