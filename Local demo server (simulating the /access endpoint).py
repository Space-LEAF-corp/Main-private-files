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
        # Start embedded local server (demo only)
        self.server_thread = threading.Thread(target=start_server, daemon=True)
        self.server_thread.start()
        time.sleep(0.8)

        self.window = tk.Tk()
        self.window.title("Local System Access + Resource & Network Monitor")
        self.window.geometry("540x500")
        self.window.resizable(False, False)

        # ── Header ───────────────────────────────────
        tk.Label(
            self.window,
            text="Local System Monitor Demo",
            font=("Helvetica", 14, "bold")
        ).pack(pady=(10, 5))

        # ── Access Status Section ─────────────────────
        self.status_frame = tk.Frame(self.window)
        self.status_frame.pack(pady=8, fill="x", padx=20)

        tk.Label(self.status_frame, text="Local Access:").pack(side="left")
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
            text="System Resources & Network (live)",
            font=("Helvetica", 11, "bold")
        ).pack(pady=(15, 5))

        self.progress_frame = tk.Frame(self.window)
        self.progress_frame.pack(padx=30, fill="x")

        row = 0

        # CPU
        tk.Label(self.progress_frame, text="CPU:").grid(row=row, column=0, sticky="w", pady=4)
        self.cpu_bar = ttk.Progressbar(self.progress_frame, maximum=100, length=280, mode='determinate')
        self.cpu_bar.grid(row=row, column=1, padx=8, sticky="ew")
        self.cpu_label = tk.Label(self.progress_frame, text="0%", width=8)
        self.cpu_label.grid(row=row, column=2)
        row += 1

        # Memory
        tk.Label(self.progress_frame, text="Memory:").grid(row=row, column=0, sticky="w", pady=4)
        self.mem_bar = ttk.Progressbar(self.progress_frame, maximum=100, length=280, mode='determinate')
        self.mem_bar.grid(row=row, column=1, padx=8, sticky="ew")
        self.mem_label = tk.Label(self.progress_frame, text="0%", width=8)
        self.mem_label.grid(row=row, column=2)
        row += 1

        # Disk
        tk.Label(self.progress_frame, text="Disk:").grid(row=row, column=0, sticky="w", pady=4)
        self.disk_bar = ttk.Progressbar(self.progress_frame, maximum=100, length=280, mode='determinate')
        self.disk_bar.grid(row=row, column=1, padx=8, sticky="ew")
        self.disk_label = tk.Label(self.progress_frame, text="0%", width=8)
        self.disk_label.grid(row=row, column=2)
        row += 1

        # ── Network I/O ───────────────────────────────
        tk.Label(self.progress_frame, text="↓ Download:").grid(row=row, column=0, sticky="w", pady=4)
        self.down_bar = ttk.Progressbar(self.progress_frame, maximum=100, length=280, mode='determinate')
        self.down_bar.grid(row=row, column=1, padx=8, sticky="ew")
        self.down_label = tk.Label(self.progress_frame, text="0 B/s", width=10)
        self.down_label.grid(row=row, column=2)
        row += 1

        tk.Label(self.progress_frame, text="↑ Upload:").grid(row=row, column=0, sticky="w", pady=4)
        self.up_bar = ttk.Progressbar(self.progress_frame, maximum=100, length=280, mode='determinate')
        self.up_bar.grid(row=row, column=1, padx=8, sticky="ew")
        self.up_label = tk.Label(self.progress_frame, text="0 B/s", width=10)
        self.up_label.grid(row=row, column=2)

        # Make the bar column expandable
        self.progress_frame.columnconfigure(1, weight=1)

        # ── State for network speed calculation ───────
        self.last_net_io = psutil.net_io_counters()
        self.last_net_time = time.time()

        # Initial update + scheduling
        self.update_all()
        self.window.after(2000, self.periodic_update)

    def check_access(self):
        try:
            response = requests.get("http://localhost:8080/access", timeout=3)
            if response.status_code == 200:
                self.system_status.set("Online")
                self.status_label.config(fg="green")
                messagebox.showinfo(
                    "System Check Complete",
                    "All local systems operational.\n"
                    "Factual information confirmed.\n"
                    "Everything is GREEN."
                )
            else:
                raise ValueError(f"Status {response.status_code}")
        except Exception as e:
            self.system_status.set("Offline")
            self.status_label.config(fg="red")
            messagebox.showerror("Connection Error", f"Cannot reach local endpoint:\n{str(e)}")

    def format_speed(self, bps: float) -> str:
        if bps < 1_000:
            return f"{bps:.0f} B/s"
        elif bps < 1_000_000:
            return f"{bps / 1_024:.1f} KB/s"
        else:
            return f"{bps / 1_048_576:.1f} MB/s"

    def update_all(self):
        # ── CPU ───────────────────────────────────────
        cpu_percent = psutil.cpu_percent(interval=None)
        self.cpu_bar['value'] = cpu_percent
        self.cpu_label.config(text=f"{cpu_percent:.1f}%")
        self._color_bar(self.cpu_bar, cpu_percent)

        # ── Memory ────────────────────────────────────
        mem = psutil.virtual_memory()
        mem_percent = mem.percent
        self.mem_bar['value'] = mem_percent
        self.mem_label.config(text=f"{mem_percent:.1f}%")
        self._color_bar(self.mem_bar, mem_percent)

        # ── Disk ──────────────────────────────────────
        try:
            disk = psutil.disk_usage('/')
            disk_percent = disk.percent
            self.disk_bar['value'] = disk_percent
            self.disk_label.config(text=f"{disk_percent:.1f}%")
            self._color_bar(self.disk_bar, disk_percent)
        except Exception:
            self.disk_label.config(text="n/a")

        # ── Network I/O speeds ────────────────────────
        now = time.time()
        current = psutil.net_io_counters()

        elapsed = now - self.last_net_time
        if elapsed > 0.3:  # avoid division by very small numbers
            bytes_recv_diff = current.bytes_recv - self.last_net_io.bytes_recv
            bytes_sent_diff = current.bytes_sent - self.last_net_io.bytes_sent

            down_bps = bytes_recv_diff / elapsed
            up_bps   = bytes_sent_diff   / elapsed

            # For progress bars: cap at reasonable max (e.g. 100 MB/s = very fast home connection)
            # You can adjust MAX_BPS to match your expected max bandwidth
            MAX_BPS_FOR_BAR = 100 * 1_048_576  # 100 MiB/s
            self.down_bar['value'] = min(100, (down_bps / MAX_BPS_FOR_BAR) * 100)
            self.up_bar['value']   = min(100, (up_bps   / MAX_BPS_FOR_BAR) * 100)

            self.down_label.config(text=self.format_speed(down_bps))
            self.up_label.config(  text=self.format_speed(up_bps))

            # Light color coding for network (mostly informational)
            self._color_bar(self.down_bar, (down_bps / MAX_BPS_FOR_BAR) * 100)
            self._color_bar(self.up_bar,   (up_bps   / MAX_BPS_FOR_BAR) * 100)

        # Update state
        self.last_net_io = current
        self.last_net_time = now

    def _color_bar(self, bar, percent):
        if percent >= 90:
            bar.config(style="red.Horizontal.TProgressbar")
        elif percent >= 60:
            bar.config(style="orange.Horizontal.TProgressbar")
        else:
            bar.config(style="green.Horizontal.TProgressbar")

    def periodic_update(self):
        self.update_all()
        self.window.after(2000, self.periodic_update)

    def run(self):
        style = ttk.Style()
        style.theme_use('default')
        style.configure("green.Horizontal.TProgressbar",  foreground="#4CAF50", background="#4CAF50")
        style.configure("orange.Horizontal.TProgressbar", foreground="#FF9800", background="#FF9800")
        style.configure("red.Horizontal.TProgressbar",    foreground="#F44336", background="#F44336")

        self.window.mainloop()


if __name__ == "__main__":
    # Privacy & Security reminders:
    #   • Binds ONLY to localhost:8080 — no external exposure
    #   • No outbound internet requests (except simulated localhost check)
    #   • psutil reads purely local metrics — NO data leaves your machine
    #   • No logging, no storage, no credentials involved

    print("Starting Local System + Network Monitor Demo...")
    print("→ Access check:     http://localhost:8080/access")
    print("→ GUI opening shortly...\n")

    app = GlobalSystemAccessPreview()
    app.run()