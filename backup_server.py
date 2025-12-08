"""Backup server for Main Private Files system.

This server runs on port 8787 and provides a backup endpoint
for the JavaScript/TypeScript backup client.
"""
from __future__ import annotations

import json
import os
from datetime import datetime
from http.server import BaseHTTPRequestHandler, HTTPServer
from typing import Dict, Any


class BackupHandler(BaseHTTPRequestHandler):
    """HTTP handler for backup requests."""

    def do_POST(self) -> None:
        """Handle POST requests to /backup endpoint."""
        if self.path == "/backup":
            try:
                # Read the request body
                content_length = int(self.headers.get("Content-Length", 0))
                body = self.rfile.read(content_length)
                backup_data = json.loads(body.decode("utf-8"))

                # Save backup to file
                result = self._save_backup(backup_data)

                # Send success response
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps(result).encode("utf-8"))

            except Exception as e:
                # Log the actual error for debugging
                print(f"[Backup Server Error] {str(e)}")
                # Send generic error response to avoid leaking sensitive information
                self.send_response(500)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                error_response = {"error": "Internal server error", "status": "failed"}
                self.wfile.write(json.dumps(error_response).encode("utf-8"))
        else:
            # Unknown path
            self.send_response(404)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(json.dumps({"error": "Not found"}).encode("utf-8"))

    def _save_backup(self, backup_data: Dict[str, Any]) -> Dict[str, Any]:
        """Save backup data to file.
        
        Args:
            backup_data: The backup bundle to save
            
        Returns:
            Result dictionary with status and backup info
        """
        # Create backups directory if it doesn't exist
        backup_dir = "backups"
        os.makedirs(backup_dir, exist_ok=True)

        # Generate filename with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"backup_{timestamp}.json"
        filepath = os.path.join(backup_dir, filename)

        # Write backup to file
        with open(filepath, "w") as f:
            json.dump(backup_data, f, indent=2)

        return {
            "status": "success",
            "message": "Backup saved successfully",
            "filename": filename,
            "timestamp": timestamp,
            "size": os.path.getsize(filepath),
        }

    def log_message(self, format: str, *args: Any) -> None:
        """Custom log message handler."""
        print(f"[Backup Server] {format % args}")


def run_backup_server(host: str = "localhost", port: int = 8787) -> None:
    """Run the backup server.
    
    Args:
        host: Host to bind to (default: localhost)
        port: Port to listen on (default: 8787)
    """
    server = HTTPServer((host, port), BackupHandler)
    print(f"🟢 Backup server running on http://{host}:{port}")
    print(f"   Backup endpoint: POST http://{host}:{port}/backup")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n🔴 Backup server stopped")
        server.shutdown()


if __name__ == "__main__":
    run_backup_server()
