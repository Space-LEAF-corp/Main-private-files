"""
Captain's Log Snapshot Creator - Python Version
Saves individual entry snapshots as formatted JSON files in captains-log/snapshots/
"""
import os
import json
import re


SNAP_DIR = "captains-log/snapshots"


def save_snapshot(entry):
    """
    Save a snapshot of an entry to a separate JSON file.
    
    Args:
        entry (dict): Entry object to snapshot
    """
    if not os.path.exists(SNAP_DIR):
        os.makedirs(SNAP_DIR, exist_ok=True)
    
    # Sanitize artifactId for use as filename
    safe_id = re.sub(r'[^a-zA-Z0-9_-]', '_', entry['artifactId'])
    file_path = os.path.join(SNAP_DIR, f"{safe_id}.json")
    
    # Write formatted JSON (pretty-printed with 2-space indentation)
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(entry, f, indent=2, ensure_ascii=False)
    
    print(f"Snapshot saved: {file_path}")
