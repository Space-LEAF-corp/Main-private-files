"""Tests for backup server."""

import json
import os
import shutil
import time
import unittest
from http.client import HTTPConnection
from threading import Thread

from backup_server import run_backup_server


class TestBackupServer(unittest.TestCase):
    """Test cases for the backup server."""

    @classmethod
    def setUpClass(cls):
        """Start the backup server in a separate thread."""
        cls.server_thread = Thread(target=run_backup_server, daemon=True)
        cls.server_thread.start()
        time.sleep(1)  # Give server time to start

    def setUp(self):
        """Set up test fixtures."""
        # Clean up any existing backups directory
        if os.path.exists("backups"):
            shutil.rmtree("backups")

    def tearDown(self):
        """Clean up after tests."""
        # Clean up backups directory
        if os.path.exists("backups"):
            shutil.rmtree("backups")

    def test_backup_endpoint_accepts_post(self):
        """Test that the /backup endpoint accepts POST requests."""
        conn = HTTPConnection("localhost", 8787)
        backup_data = {
            "timestamp": "2024-01-01T00:00:00.000Z",
            "version": "1.0.0",
            "data": {
                "users": ["alice", "bob"],
                "sessions": [],
                "firewall_state": {},
                "metadata": {"source": "test"},
            },
        }

        conn.request(
            "POST",
            "/backup",
            body=json.dumps(backup_data),
            headers={"Content-Type": "application/json"},
        )
        response = conn.getresponse()

        self.assertEqual(response.status, 200)
        response_data = json.loads(response.read().decode())
        self.assertEqual(response_data["status"], "success")
        self.assertIn("filename", response_data)
        conn.close()

    def test_backup_creates_file(self):
        """Test that backup creates a file in the backups directory."""
        conn = HTTPConnection("localhost", 8787)
        backup_data = {
            "timestamp": "2024-01-01T00:00:00.000Z",
            "version": "1.0.0",
            "data": {"users": [], "sessions": [], "firewall_state": {}, "metadata": {}},
        }

        conn.request(
            "POST",
            "/backup",
            body=json.dumps(backup_data),
            headers={"Content-Type": "application/json"},
        )
        response = conn.getresponse()
        response_data = json.loads(response.read().decode())
        conn.close()

        # Check that backup file was created
        self.assertTrue(os.path.exists("backups"))
        backup_files = os.listdir("backups")
        self.assertGreater(len(backup_files), 0)
        self.assertTrue(backup_files[0].startswith("backup_"))
        self.assertTrue(backup_files[0].endswith(".json"))

    def test_backup_saves_correct_data(self):
        """Test that backup saves the correct data."""
        conn = HTTPConnection("localhost", 8787)
        backup_data = {
            "timestamp": "2024-01-01T00:00:00.000Z",
            "version": "1.0.0",
            "data": {
                "users": ["test_user"],
                "sessions": ["session_123"],
                "firewall_state": {"active": True},
                "metadata": {"test": "data"},
            },
        }

        conn.request(
            "POST",
            "/backup",
            body=json.dumps(backup_data),
            headers={"Content-Type": "application/json"},
        )
        response = conn.getresponse()
        response_data = json.loads(response.read().decode())
        conn.close()

        # Read the backup file
        backup_files = os.listdir("backups")
        with open(os.path.join("backups", backup_files[0])) as f:
            saved_data = json.load(f)

        # Verify the data matches
        self.assertEqual(saved_data["timestamp"], backup_data["timestamp"])
        self.assertEqual(saved_data["version"], backup_data["version"])
        self.assertEqual(saved_data["data"]["users"], backup_data["data"]["users"])

    def test_invalid_path_returns_404(self):
        """Test that invalid paths return 404."""
        conn = HTTPConnection("localhost", 8787)
        conn.request("POST", "/invalid", body="{}")
        response = conn.getresponse()

        self.assertEqual(response.status, 404)
        conn.close()


if __name__ == "__main__":
    unittest.main()
