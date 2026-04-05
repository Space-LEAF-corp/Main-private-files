import ui
import requests
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

# -------------------------
# Mock local server (runs inside Pythonista)
# -------------------------

class MockHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/scan":
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b"OK")
        else:
            self.send_response(404)
            self.end_headers()

    def do_POST(self):
        if self.path == "/alert":
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b"ALERT SENT")
        else:
            self.send_response(404)
            self.end_headers()

def start_server():
    server = HTTPServer(("localhost", 8080), MockHandler)
    server.serve_forever()

# Start server in background thread
threading.Thread(target=start_server, daemon=True).start()

# -------------------------
# Pythonista UI
# -------------------------

class ReversePulseScanShield(ui.View):
    def __init__(self):
        self.name = "Reverse Pulse Scan Shield"
        self.background_color = "black"

        self.scan_button = ui.Button(
            title="Scan",
            bg_color="#4CAF50",
            tint_color="white",
            action=self.scan_action
        )
        self.scan_button.frame = (20, 50, 200, 50)
        self.add_subview(self.scan_button)

        self.alert_button = ui.Button(
            title="Send Alert",
            bg_color="#F44336",
            tint_color="white",
            action=self.alert_action
        )
        self.alert_button.frame = (20, 120, 200, 50)
        self.add_subview(self.alert_button)

    def scan_action(self, sender):
        try:
            r = requests.get("http://localhost:8080/scan", timeout=2)
            if r.status_code == 200:
                ui.alert("Scan Results", "All clear!", "OK")
            else:
                ui.alert("Scan Results", "Threat detected!", "OK")
        except Exception as e:
            ui.alert("Error", f"Could not reach server:\n{e}", "OK")

    def alert_action(self, sender):
        try:
            r = requests.post("http://localhost:8080/alert", timeout=2)
            if r.status_code == 200:
                ui.alert("Alert Sent", "Alert sent successfully!", "OK")
            else:
                ui.alert("Alert Sent", "Failed to send alert!", "OK")
        except Exception as e:
            ui.alert("Error", f"Could not reach server:\n{e}", "OK")


view = ReversePulseScanShield()
view.present("sheet")