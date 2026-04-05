import ui
import requests
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
import time

# -------------------------
# Mock local server
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

threading.Thread(target=start_server, daemon=True).start()

# -------------------------
# Cyber aesthetic theme system
# -------------------------

THEMES = {
    "Nebula Blue": {"glow": "#00C8FF", "bg": "#000000"},
    "Solar Amber": {"glow": "#FFB300", "bg": "#000000"},
    "Panther Violet": {"glow": "#C060FF", "bg": "#000000"},
    "Ghostline White": {"glow": "#E0E0E0", "bg": "#000000"},
}

current_theme = "Panther Violet"

# -------------------------
# Cyber Terminal UI
# -------------------------

class CyberLauncher(ui.View):
    def __init__(self):
        self.name = "Reverse Pulse Scan Shield"
        self.flex = "WH"
        self.background_color = THEMES[current_theme]["bg"]

        self.header = ui.Label(
            text="REVERSE PULSE SCAN SHIELD — G‐NET ACCESS",
            text_color=THEMES[current_theme]["glow"],
            alignment=1,
            font=("Menlo-Bold", 14)
        )
        self.header.frame = (0, 20, self.width, 30)
        self.header.flex = "W"
        self.add_subview(self.header)

        self.line = ui.View(
            background_color=THEMES[current_theme]["glow"]
        )
        self.line.frame = (20, 55, self.width - 40, 1)
        self.line.flex = "W"
        self.add_subview(self.line)

        self.buttons = []
        labels = ["SCAN", "SEND ALERT", "SYSTEM STATUS", "SETTINGS"]
        for i, label in enumerate(labels):
            btn = ui.Button(
                title=label,
                bg_color="#111111",
                tint_color=THEMES[current_theme]["glow"],
                font=("Menlo-Bold", 16),
                action=self.button_action
            )
            btn.frame = (40, 120 + i * 70, self.width - 80, 50)
            btn.corner_radius = 8
            btn.flex = "W"
            self.add_subview(btn)
            self.buttons.append(btn)

        self.footer = ui.Label(
            text=f"Theme: {current_theme}   |   v1.0",
            text_color=THEMES[current_theme]["glow"],
            alignment=1,
            font=("Menlo", 10)
        )
        self.footer.frame = (0, self.height - 40, self.width, 20)
        self.footer.flex = "WT"
        self.add_subview(self.footer)

        ui.animate(self.pulse_header, duration=0.8)

    def pulse_header(self):
        self.header.alpha = 0.4
        ui.animate(lambda: setattr(self.header, "alpha", 1.0), duration=0.8)

    def button_action(self, sender):
        if sender.title == "SCAN":
            self.run_scan()
        elif sender.title == "SEND ALERT":
            self.run_alert()
        elif sender.title == "SYSTEM STATUS":
            ui.alert("System Status", "All systems nominal.", "OK")
        elif sender.title == "SETTINGS":
            self.open_settings()

    def run_scan(self):
        try:
            r = requests.get("http://localhost:8080/scan", timeout=2)
            if r.status_code == 200:
                ui.alert("Scan", "All clear.", "OK")
            else:
                ui.alert("Scan", "Threat detected.", "OK")
        except Exception as e:
            ui.alert("Error", str(e), "OK")

    def run_alert(self):
        try:
            r = requests.post("http://localhost:8080/alert", timeout=2)
            if r.status_code == 200:
                ui.alert("Alert", "Alert sent.", "OK")
            else:
                ui.alert("Alert", "Failed.", "OK")
        except Exception as e:
            ui.alert("Error", str(e), "OK")

    def open_settings(self):
        theme_names = list(THEMES.keys())
        def theme_chosen(i):
            global current_theme
            if i >= 0:
                current_theme = theme_names[i]
                self.apply_theme()

        ui.ActionSheet(
            title="Select Theme",
            actions=theme_names,
            callback=theme_chosen
        ).present()

    def apply_theme(self):
        glow = THEMES[current_theme]["glow"]
        self.background_color = THEMES[current_theme]["bg"]
        self.header.text_color = glow
        self.line.background_color = glow
        self.footer.text_color = glow
        for b in self.buttons:
            b.tint_color = glow
        self.footer.text = f"Theme: {current_theme}   |   v1.0"

view = CyberLauncher()
view.present("fullscreen")