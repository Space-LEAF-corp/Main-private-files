"""
Retrieve Logs API
-----------------
Provides access to stored GQI&EA logs.
"""

import os

LOG_DIR = "data/logs/"

def retrieve_logs() -> list:
    """
    Return a list of all log entries.
    Placeholder implementation.
    """
    logs = []
    for filename in os.listdir(LOG_DIR):
        path = os.path.join(LOG_DIR, filename)
        with open(path, "r") as f:
            logs.append(f.read())
    return logs
