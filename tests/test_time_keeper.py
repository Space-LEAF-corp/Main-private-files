"""Tests for TimeKeeper AI."""
import unittest
import numpy as np
from time_keeper import TimeKeeper


class TestTimeKeeperInitialization(unittest.TestCase):
    def test_initialization_defaults(self):
        tk = TimeKeeper()
        self.assertEqual(tk.creator, "Leif William Sogge")
        self.assertEqual(tk.status, "Dormant")
        self.assertFalse(tk.alterable)
        self.assertFalse(tk.weaponizable)
        self.assertTrue(tk._firewall_active)
        self.assertIsNone(tk.sibling_ai)

    def test_initialization_custom_creator(self):
        tk = TimeKeeper(creator="Custom Creator")
        self.assertEqual(tk.creator, "Custom Creator")

    def test_core_laws_exist(self):
        tk = TimeKeeper()
        self.assertEqual(len(tk.core_laws), 8)
        self.assertIn("Logic is not good or bad, only knowledge to process.", tk.core_laws)
        self.assertIn("The AI cannot and will never be weaponized.", tk.core_laws)

    def test_encryption_key_generated(self):
        tk = TimeKeeper()
        self.assertIsNotNone(tk.encryption_key)
        self.assertIsNotNone(tk.cipher)


class TestTimeKeeperActivation(unittest.TestCase):
    def test_activation_with_valid_clearance(self):
        tk = TimeKeeper()
        tk.activate(clearance_code="valid")
        self.assertEqual(tk.status, "Active")

    def test_activation_with_invalid_clearance(self):
        tk = TimeKeeper()
        with self.assertRaises(PermissionError) as context:
            tk.activate(clearance_code="invalid")
        self.assertIn("Activation denied", str(context.exception))


class TestTimeKeeperContentFiltering(unittest.TestCase):
    def test_filter_content_safe(self):
        tk = TimeKeeper()
        result = tk.filter_content("This is safe content.")
        self.assertTrue(result)

    def test_filter_content_prohibited(self):
        tk = TimeKeeper()
        prohibited_terms = ["porn", "adult", "xxx", "explicit", "nudity", "sex", "nsfw"]
        for term in prohibited_terms:
            with self.assertRaises(ValueError) as context:
                tk.filter_content(f"This contains {term} content.")
            self.assertIn("Blocked content detected", str(context.exception))

    def test_filter_content_case_insensitive(self):
        tk = TimeKeeper()
        with self.assertRaises(ValueError):
            tk.filter_content("This contains PORN in uppercase.")


class TestTimeKeeperEncryption(unittest.TestCase):
    def test_encryption_and_decryption(self):
        tk = TimeKeeper()
        original_data = "Secret message"
        encrypted = tk.encryption_layer(original_data)
        self.assertNotEqual(encrypted, original_data.encode())
        
        decrypted = tk.decryption_layer(encrypted)
        self.assertEqual(decrypted, original_data)

    def test_encryption_produces_bytes(self):
        tk = TimeKeeper()
        encrypted = tk.encryption_layer("Test data")
        self.assertIsInstance(encrypted, bytes)


class TestTimeKeeperMatrixOperations(unittest.TestCase):
    def test_rotate_matrix(self):
        tk = TimeKeeper()
        matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        rotated = tk.rotate_matrix(matrix)
        expected = np.array([[7, 4, 1], [8, 5, 2], [9, 6, 3]])
        np.testing.assert_array_equal(rotated, expected)

    def test_flip_matrix(self):
        tk = TimeKeeper()
        matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        flipped = tk.flip_matrix(matrix)
        expected = np.array([[3, 2, 1], [6, 5, 4], [9, 8, 7]])
        np.testing.assert_array_equal(flipped, expected)

    def test_kaleidoscope_encrypt(self):
        tk = TimeKeeper()
        data = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        encrypted = tk.kaleidoscope_encrypt(data)
        self.assertEqual(encrypted.shape, (3, 3))
        self.assertIsInstance(encrypted, np.ndarray)

    def test_process_in_orbit(self):
        tk = TimeKeeper()
        # Create a simple cube with identity-like matrices
        cube = [np.eye(3), np.eye(3)]
        data = np.array([1, 2, 3])
        result = tk.process_in_orbit(cube, data)
        self.assertIsInstance(result, np.ndarray)


class TestTimeKeeperImmutableLaws(unittest.TestCase):
    def test_enforce_rules_success(self):
        tk = TimeKeeper()
        # Should not raise any exceptions
        tk.enforce_rules()

    def test_weaponizable_remains_false(self):
        tk = TimeKeeper()
        self.assertFalse(tk.weaponizable)
        # Verify enforce_rules would fail if weaponizable were True
        tk_modified = TimeKeeper()
        tk_modified.weaponizable = True
        with self.assertRaises(AssertionError) as context:
            tk_modified.enforce_rules()
        self.assertIn("AI cannot ever be weaponized", str(context.exception))

    def test_alterable_remains_false(self):
        tk = TimeKeeper()
        self.assertFalse(tk.alterable)
        # Verify enforce_rules would fail if alterable were True
        tk_modified = TimeKeeper()
        tk_modified.alterable = True
        with self.assertRaises(AssertionError) as context:
            tk_modified.enforce_rules()
        self.assertIn("AI cannot be altered except by its creator", str(context.exception))


class TestTimeKeeperLogging(unittest.TestCase):
    def test_log_interaction(self):
        tk = TimeKeeper()
        # This should not raise any exceptions
        tk.log_interaction("test_interaction", "success")
        tk.log_interaction("another_interaction", "failure")


if __name__ == "__main__":
    unittest.main()
