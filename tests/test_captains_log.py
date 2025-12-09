"""Unit tests for Captain's Log functionality."""
import os
import json
import tempfile
import shutil
import unittest

# Add parent directory to path to import modules
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from save_log import save_entry, add_checksum, ensure_dir
from snapshot import save_snapshot


class TestCaptainsLog(unittest.TestCase):
    """Test suite for captain's log entry saving and snapshot functionality."""
    
    def setUp(self):
        """Set up temporary test directories."""
        self.test_dir = tempfile.mkdtemp()
        self.old_cwd = os.getcwd()
        os.chdir(self.test_dir)
        
    def tearDown(self):
        """Clean up temporary test directories."""
        os.chdir(self.old_cwd)
        shutil.rmtree(self.test_dir)
    
    def test_ensure_dir_creates_directory(self):
        """Test that ensure_dir creates a new directory."""
        test_path = os.path.join(self.test_dir, "test_captains_log")
        self.assertFalse(os.path.exists(test_path))
        ensure_dir(test_path)
        self.assertTrue(os.path.exists(test_path))
    
    def test_add_checksum_adds_sha256_hash(self):
        """Test that add_checksum adds a valid SHA-256 hash."""
        entry = {
            "artifactId": "test-001",
            "title": "Test Entry",
            "timestamp": "2025-12-08T00:00:00",
            "author": "Test Author"
        }
        entry_with_checksum = add_checksum(entry.copy())
        self.assertIn("checksum_sha256", entry_with_checksum)
        self.assertEqual(len(entry_with_checksum["checksum_sha256"]), 64)  # SHA-256 hex length
    
    def test_add_checksum_is_deterministic(self):
        """Test that checksum is deterministic for same data."""
        entry = {
            "artifactId": "test-001",
            "title": "Test Entry",
            "timestamp": "2025-12-08T00:00:00"
        }
        checksum1 = add_checksum(entry.copy())["checksum_sha256"]
        checksum2 = add_checksum(entry.copy())["checksum_sha256"]
        self.assertEqual(checksum1, checksum2)
    
    def test_save_entry_creates_log_file(self):
        """Test that save_entry creates the log file."""
        entry = {
            "artifactId": "ark-001",
            "title": "Test Ark Entry",
            "author": "Test Captain"
        }
        entry = add_checksum(entry)
        save_entry(entry)
        
        log_file = os.path.join("captains-log", "log.jsonl")
        self.assertTrue(os.path.exists(log_file))
    
    def test_save_entry_appends_jsonl_format(self):
        """Test that entries are saved in JSONL format."""
        entry1 = {
            "artifactId": "ark-001",
            "title": "First Entry"
        }
        entry2 = {
            "artifactId": "ark-002",
            "title": "Second Entry"
        }
        
        entry1 = add_checksum(entry1)
        entry2 = add_checksum(entry2)
        
        save_entry(entry1)
        save_entry(entry2)
        
        log_file = os.path.join("captains-log", "log.jsonl")
        with open(log_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        self.assertEqual(len(lines), 2)
        # Verify each line is valid JSON
        parsed1 = json.loads(lines[0])
        parsed2 = json.loads(lines[1])
        self.assertEqual(parsed1["artifactId"], "ark-001")
        self.assertEqual(parsed2["artifactId"], "ark-002")
    
    def test_save_snapshot_creates_snapshot_file(self):
        """Test that save_snapshot creates a snapshot JSON file."""
        entry = {
            "artifactId": "ark-001",
            "title": "Test Snapshot",
            "author": "Test Captain"
        }
        save_snapshot(entry)
        
        snapshot_file = os.path.join("captains-log", "snapshots", "ark-001.json")
        self.assertTrue(os.path.exists(snapshot_file))
    
    def test_save_snapshot_sanitizes_artifact_id(self):
        """Test that snapshot sanitizes special characters in artifact ID."""
        entry = {
            "artifactId": "ark:001/test",
            "title": "Test Snapshot"
        }
        save_snapshot(entry)
        
        # Should convert : and / to _
        snapshot_file = os.path.join("captains-log", "snapshots", "ark_001_test.json")
        self.assertTrue(os.path.exists(snapshot_file))
    
    def test_save_snapshot_contains_formatted_json(self):
        """Test that snapshot contains properly formatted JSON."""
        entry = {
            "artifactId": "ark-001",
            "title": "Test Snapshot",
            "nested": {"key": "value"}
        }
        save_snapshot(entry)
        
        snapshot_file = os.path.join("captains-log", "snapshots", "ark-001.json")
        with open(snapshot_file, 'r', encoding='utf-8') as f:
            content = f.read()
            parsed = json.loads(content)
        
        self.assertEqual(parsed["artifactId"], "ark-001")
        self.assertEqual(parsed["nested"]["key"], "value")
        # Check that it's formatted (has newlines for pretty-printing)
        self.assertIn('\n', content)
    
    def test_complete_ark_001_entry(self):
        """Test saving the complete ark-001 entry from the issue."""
        entry = {
            "artifactId": "ark-001",
            "title": "The university became a moving ship",
            "timestamp": "2025-11-13T21:11:00-05:00",
            "author": "Leif William Sogge",
            "ceremonialContext": {
                "vessel": "Ark",
                "covenant": "Protection and peace",
                "status": "Consecrated and complete"
            },
            "entry": "On this day, the university becomes a moving ship—an Ark. Everyone aboard is safe.",
            "artifacts": [
                {
                    "type": "image",
                    "filename": "9ZYTqXhtqaFpyp6av5bHi.jpeg",
                    "label": "Water-port vessel",
                    "notes": "Twin circular engines; guardians stand watch."
                },
                {
                    "type": "image",
                    "filename": "F4TG8oJMN7gomr2TbCeWS.jpeg",
                    "label": "Floating city above ocean",
                    "notes": "Waterfalls cascade; sanctuary architecture."
                },
                {
                    "type": "image",
                    "filename": "tcSw9D3UpZBHqcBubnMuo.jpeg",
                    "label": "Hybrid Pokémon (name reserved)",
                    "notes": "Joy and flame—a guardian of balance."
                }
            ],
            "seals": {
                "externalAuthority": "God",
                "localSealApplied": False
            }
        }
        
        # Add checksum and save
        entry = add_checksum(entry)
        save_entry(entry)
        save_snapshot(entry)
        
        # Verify log entry
        log_file = os.path.join("captains-log", "log.jsonl")
        self.assertTrue(os.path.exists(log_file))
        with open(log_file, 'r', encoding='utf-8') as f:
            saved_entry = json.loads(f.readline())
        
        self.assertEqual(saved_entry["artifactId"], "ark-001")
        self.assertEqual(saved_entry["author"], "Leif William Sogge")
        self.assertIn("checksum_sha256", saved_entry)
        self.assertEqual(len(saved_entry["artifacts"]), 3)
        
        # Verify snapshot
        snapshot_file = os.path.join("captains-log", "snapshots", "ark-001.json")
        self.assertTrue(os.path.exists(snapshot_file))


if __name__ == '__main__':
    unittest.main()
