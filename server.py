"""Local dual-layer server: Socket + HTTP with authentication.

This server provides two connection methods:
1. Socket (TCP) layer — low-level, binary protocol
2. HTTP layer — REST-like, JSON responses

Both require authentication via the AuthManager.
"""
from __future__ import annotations

import json
import socket
import threading
import time
from http.server import BaseHTTPRequestHandler, HTTPServer
from typing import Dict, Optional
from urllib.parse import parse_qs, urlparse

from auth import AuthManager
from secured_firewall import SecuredFirewall


class DualLayerServer:
    def __init__(self, host: str = "localhost", socket_port: int = 9000, http_port: int = 9001):
        self.host = host
        self.socket_port = socket_port
        self.http_port = http_port
        self.firewall = SecuredFirewall()
        self.running = False

    def start(self) -> None:
        self.running = True
        # Start socket server in a thread
        socket_thread = threading.Thread(target=self._socket_server, daemon=True)
        socket_thread.start()
        # Start HTTP server in a thread
        http_thread = threading.Thread(target=self._http_server, daemon=True)
        http_thread.start()
        print(f"🟢 Dual-layer server online")
        print(f"   Socket layer: {self.host}:{self.socket_port}")
        print(f"   HTTP layer: http://{self.host}:{self.http_port}")

    def _socket_server(self) -> None:
        """Simple socket server (TCP) for direct binary communication."""
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        sock.bind((self.host, self.socket_port))
        sock.listen(5)
        print(f"Socket server listening on {self.host}:{self.socket_port}")

        while self.running:
            try:
                sock.settimeout(1)
                conn, addr = sock.accept()
                print(f"Socket connection from {addr}")
                self._handle_socket_client(conn)
            except socket.timeout:
                continue
            except Exception as e:
                print(f"Socket error: {e}")

    def _handle_socket_client(self, conn: socket.socket) -> None:
        try:
            # Simple protocol: send command as JSON, receive response as JSON
            data = conn.recv(1024).decode("utf-8")
            if not data:
                conn.close()
                return

            req = json.loads(data)
            cmd = req.get("cmd")
            res = self._process_command(cmd, req)
            conn.sendall(json.dumps(res).encode("utf-8"))
            conn.close()
        except Exception as e:
            conn.sendall(json.dumps({"error": str(e)}).encode("utf-8"))
            conn.close()

    def _http_server(self) -> None:
        """HTTP server for REST-like API."""
        handler = self._make_http_handler()
        server = HTTPServer((self.host, self.http_port), handler)
        print(f"HTTP server listening on http://{self.host}:{self.http_port}")
        while self.running:
            server.handle_request()

    def _make_http_handler(self):
        server = self

        class DualLayerHTTPHandler(BaseHTTPRequestHandler):
            def do_GET(self):
                parsed = urlparse(self.path)
                path = parsed.path
                query = parse_qs(parsed.query)

                if path == "/register":
                    user_id = query.get("user_id", [""])[0]
                    password = query.get("password", [""])[0]
                    dna = query.get("dna", [""])[0]
                    res = server.firewall.register_user(user_id, password, dna)
                elif path == "/login":
                    user_id = query.get("user_id", [""])[0]
                    password = query.get("password", [""])[0]
                    dna = query.get("dna", [""])[0]
                    res = server.firewall.login_step1(user_id, password, dna)
                elif path == "/login_otp":
                    user_id = query.get("user_id", [""])[0]
                    otp = query.get("otp", [""])[0]
                    res = server.firewall.login_step2(user_id, otp)
                elif path == "/access":
                    session = query.get("session_token", [""])[0]
                    res = server.firewall.access_firewall(session)
                else:
                    res = {"error": "unknown path"}

                self.send_response(200)
                self.send_header("Content-type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps(res).encode("utf-8"))

            def log_message(self, format, *args):
                pass  # suppress default logging

        return DualLayerHTTPHandler

    def _process_command(self, cmd: str, req: Dict) -> Dict:
        """Process socket commands."""
        if cmd == "register":
            return self.firewall.register_user(
                req.get("user_id", ""),
                req.get("password", ""),
                req.get("dna", ""),
            )
        elif cmd == "login":
            return self.firewall.login_step1(
                req.get("user_id", ""),
                req.get("password", ""),
                req.get("dna", ""),
            )
        elif cmd == "login_otp":
            return self.firewall.login_step2(req.get("user_id", ""), req.get("otp", ""))
        elif cmd == "access":
            return self.firewall.access_firewall(req.get("session_token", ""))
        else:
            return {"error": f"unknown command: {cmd}"}

    def stop(self) -> None:
        self.running = False
