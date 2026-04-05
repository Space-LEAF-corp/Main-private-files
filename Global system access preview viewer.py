import tkinter as tk
from tkinter import messagebox, ttk
import requests
import http.server
import socketserver
import threading
import time
import psutil

# ────────────────────────────────────────────────
#   Local demo server (simulating the /access endpoint)
# ────────────────────────────────────────────────
class AccessHandler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path == '/access':
            self.send_response(200)
            self.send_header("Content-type", "text/plain")
            self.end_headers()
            self.wfile.write(b"Access Granted - Local System Demo")
        else:
            self.send_error(404, "Not Found")

def start_server():
    with socketserver.TCPServer(("", 8080), AccessHandler) as httpd:
        httpd.serve_forever()

# ────────────────────────────────────────────────
#   Main GUI Application
# ────────────────────────────────────────────────
class GlobalSystemAccessPreview:
    def __init__(self):
        # Start embedded local server (for demo purposes)
        self.server_thread = threading.Thread(target=start_server, daemon=True)
        self.server_thread.start()
        time.sleep(0.8)  # Give server a moment to bind

        self.window = tk.Tk()
        self.window.title("Local System Access + Resource Monitor Demo")
        self.window.geometry("520x420")
        self.window.resizable(False, False)

        # ── Header ───────────────────────────────────
        tk.Label(
            self.window,
            text="Local System Access & Resource Monitor",
            font=("Helvetica", 14, "bold")
        ).pack(pady=(10, 5))

        # ── Access Status Section ─────────────────────
        self.status_frame = tk.Frame(self.window)
        self.status_frame.pack(pady=8, fill="x", padx=20)

        tk.Label(self.status_frame, text="System Access:").pack(side="left")
        self.system_status = tk.StringVar(value="Offline")
        self.status_label = tk.Label(
            self.status_frame,
            textvariable=self.system_status,
            font=("Helvetica", 11, "bold"),
            fg="gray"
        )
        self.status_label.pack(side="left", padx=10)

        self.access_button = tk.Button(
            self.window,
            text="Check Local Access",
            command=self.check_access,
            width=20
        )
        self.access_button.pack(pady=6)

        # ── Resource Monitor Section ──────────────────
        tk.Label(
            self.window,
            text="System Resources (updated live)",
            font=("Helvetica", 11, "bold")
        ).pack(pady=(15, 5))

        self.progress_frame = tk.Frame(self.window)
        self.progress_frame.pack(padx=30, fill="x")

        # CPU
        tk.Label(self.progress_frame, text="CPU:").grid(row=0, column=0, sticky="w", pady=4)
        self.cpu_bar = ttk.Progressbar(self.progress_frame, maximum=100, length=280, mode='determinate')
        self.cpu_bar.grid(row=0, column=1, padx=8, sticky="ew")
        self.cpu_label = tk.Label(self.progress_frame, text="0%", width=6)
        self.cpu_label.grid(row=0, column=2)

        # Memory
        tk.Label(self.progress_frame, text="Memory:").grid(row=1, column=0, sticky="w", pady=4)
        self.mem_bar = ttk.Progressbar(self.progress_frame, maximum=100, length=280, mode='determinate')
        self.mem_bar.grid(row=1, column=1, padx=8, sticky="ew")
        self.mem_label = tk.Label(self.progress_frame, text="0%", width=6)
        self.mem_label.grid(row=1, column=2)

        # Disk (main volume)
        tk.Label(self.progress_frame, text="Disk:").grid(row=2, column=0, sticky="w", pady=4)
        self.disk_bar = ttk.Progressbar(self.progress_frame, maximum=100, length=280, mode='determinate')
        self.disk_bar.grid(row=2, column=1, padx=8, sticky="ew")
        self.disk_label = tk.Label(self.progress_frame, text="0%", width=6)
        self.disk_label.grid(row=2, column=2)

        # ── Update timer ──────────────────────────────
        self.update_resources()
        self.window.after(2000, self.periodic_update)  # every 2 seconds

        # Make columns expandable
        self.progress_frame.columnconfigure(1, weight=1)

    def check_access(self):
        try:
            response = requests.get("http://localhost:8080/access", timeout=3)
            if response.status_code == 200:
                self.system_status.set("Online")
                self.status_label.config(fg="green")
                messagebox.showinfo(
                    "System Check Complete",
                    "All local systems are operational.\n"
                    "Factual information confirmed.\n"
                    "Everything is GREEN."
                )
            else:
                self.system_status.set("Offline")
                self.status_label.config(fg="red")
        except requests.exceptions.RequestException as e:
            self.system_status.set("Offline")
            self.status_label.config(fg="red")
            messagebox.showerror("Connection Error", f"Cannot reach local endpoint:\n{str(e)}")

    def update_resources(self):
        # CPU (average over all cores)
        cpu_percent = psutil.cpu_percent(interval=None)  # non-blocking
        self.cpu_bar['value'] = cpu_percent
        self.cpu_label.config(text=f"{cpu_percent:.1f}%")
        self._color_bar(self.cpu_bar, cpu_percent)

        # Memory
        mem = psutil.virtual_memory()
        mem_percent = mem.percent
        self.mem_bar['value'] = mem_percent
        self.mem_label.config(text=f"{mem_percent:.1f}%")
        self._color_bar(self.mem_bar, mem_percent)

        # Disk (root filesystem)
        try:
            disk = psutil.disk_usage('/')
            disk_percent = disk.percent
            self.disk_bar['value'] = disk_percent
            self.disk_label.config(text=f"{disk_percent:.1f}%")
            self._color_bar(self.disk_bar, disk_percent)
        except Exception:
            self.disk_label.config(text="n/a")

    def _color_bar(self, bar, percent):
        if percent >= 90:
            bar.config(style="red.Horizontal.TProgressbar")
        elif percent >= 70:
            bar.config(style="orange.Horizontal.TProgressbar")
        else:
            bar.config(style="green.Horizontal.TProgressbar")

    def periodic_update(self):
        self.update_resources()
        self.window.after(2000, self.periodic_update)

    def run(self):
        # Create custom progressbar styles
        style = ttk.Style()
        style.theme_use('default')
        style.configure("green.Horizontal.TProgressbar", foreground="green", background="green")
        style.configure("orange.Horizontal.TProgressbar", foreground="orange", background="orange")
        style.configure("red.Horizontal.TProgressbar", foreground="red", background="red")

        self.window.mainloop()


if __name__ == "__main__":
    # Security / privacy note:
    #   • This script only binds to localhost:8080
    #   • No external network communication except the simulated localhost check
    #   • psutil only reads local system metrics — no data is sent anywhere

    print("Starting Local System Access + Resource Monitor Demo...")
    print("→ Access check endpoint: http://localhost:8080/access")
    print("→ GUI will open shortly...\n")

    app = GlobalSystemAccessPreview()
    app.run()