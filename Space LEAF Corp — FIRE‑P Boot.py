import http.server
import socketserver
import webbrowser
import threading

PORT = 8080

HTML = r"""
<!DOCTYPE html>
<html>
<head>
<title>Space LEAF Corp — FIRE‑P Boot</title>
<style>
body {
    font-family: Arial, sans-serif;
    background: #000;
    color: #eee;
    padding: 40px;
}
h1 { color: #ff6600; }
h2 { color: #66ccff; }
input, select {
    padding: 10px;
    font-size: 18px;
    width: 320px;
    margin-top: 10px;
}
button {
    padding: 10px 20px;
    font-size: 18px;
    margin-top: 15px;
    cursor: pointer;
}
.screen {
    display: none;
}
#screen-boot {
    display: block;
}
.panel {
    background: #111;
    padding: 20px;
    border-radius: 8px;
    margin-top: 20px;
}
ul {
    margin-top: 10px;
}
</style>
<script>
let selectedMode = null;

function speak(text) {
    if (!("speechSynthesis" in window)) return;
    const utter = new SpeechSynthesisUtterance(text);
    speechSynthesis.speak(utter);
}

function showScreen(id) {
    document.getElementById("screen-boot").style.display = "none";
    document.getElementById("screen-mode").style.display = "none";
    document.getElementById("screen-shell").style.display = "none";
    document.getElementById(id).style.display = "block";
}

function startFireP() {
    const val = document.getElementById("code").value.trim().toUpperCase();
    if (val === "FIRE-P" || val === "FIREP") {
        showScreen("screen-mode");
        speak(
            "Welcome to Space LEAF Corp. You have activated FIRE P, the Fragile Infrastructure and Residence Ecosystem Pulse. " +
            "Accessibility startup mode is now enabled. Select your operating mode to continue."
        );
    } else {
        alert("Invalid startup code. Type FIRE-P to begin.");
    }
}

function onModeChange() {
    const sel = document.getElementById("mode-select");
    selectedMode = sel.value || null;
    document.getElementById("enter-btn").disabled = !selectedMode;
}

function enterShell() {
    if (!selectedMode) {
        alert("Please select an operating mode first.");
        return;
    }
    const title = document.getElementById("shell-title");
    const body = document.getElementById("shell-body");

    if (selectedMode === "firep") {
        title.textContent = "🔥 FIRE‑P Environmental Pulse Dashboard";
        body.innerHTML = `
            <p>Monitoring fragile infrastructure and residence ecosystem pulse.</p>
            <ul>
                <li>Air quality: Stable (GREEN)</li>
                <li>Power grid: Normal load</li>
                <li>Water systems: No alerts</li>
                <li>Wildlife / ecosystem: No disruptions detected</li>
            </ul>
        `;
        speak("FIRE P mode loaded. Environmental pulse dashboard active.");
    } else if (selectedMode === "spaceleaf") {
        title.textContent = "🌿 Space LEAF Corp OS Shell";
        body.innerHTML = `
            <p>Welcome to the ceremonial, privacy-first Space LEAF Corp environment.</p>
            <ul>
                <li>Captain's Log</li>
                <li>Ceremonies & Oaths</li>
                <li>Whisker's Guidance System</li>
                <li>Council Charter & Flags</li>
            </ul>
        `;
        speak("Space LEAF Corp operating shell loaded.");
    } else if (selectedMode === "apple") {
        title.textContent = "🍎 Apple Companion Mode (Mock)";
        body.innerHTML = `
            <p>This mock mode would open Apple account or settings on a real device.</p>
            <ul>
                <li>Open iCloud portal</li>
                <li>Accessibility shortcuts</li>
                <li>Device settings</li>
            </ul>
        `;
        speak("Apple companion mode mock loaded.");
    } else if (selectedMode === "microsoft") {
        title.textContent = "🪟 Microsoft Companion Mode (Mock)";
        body.innerHTML = `
            <p>This mock mode would open Microsoft account or accessibility center.</p>
            <ul>
                <li>Microsoft account portal</li>
                <li>Xbox / PC accessibility</li>
                <li>Cloud services</li>
            </ul>
        `;
        speak("Microsoft companion mode mock loaded.");
    } else if (selectedMode === "tesla") {
        title.textContent = "🚗 Tesla Companion Mode (Mock)";
        body.innerHTML = `
            <p>This mock mode would open the Tesla app or web portal.</p>
            <ul>
                <li>Vehicle status</li>
                <li>Charging overview</li>
                <li>Location & climate controls</li>
            </ul>
        `;
        speak("Tesla companion mode mock loaded.");
    }

    showScreen("screen-shell");
}

window.addEventListener("DOMContentLoaded", () => {
    document.getElementById("enter-btn").disabled = true;
});
</script>
</head>
<body>

<!-- Screen 1: FIRE-P Boot -->
<div id="screen-boot" class="screen">
    <h1>🔥 FIRE‑P Startup</h1>
    <p>Type <strong>FIRE‑P</strong> to begin system initialization.</p>
    <input id="code" placeholder="Enter FIRE-P" />
    <br>
    <button onclick="startFireP()">Start</button>

    <div class="panel">
        <h2>Accessibility</h2>
        <p>On startup, FIRE‑P enables:</p>
        <ul>
            <li>Audio narration (where supported)</li>
            <li>High‑contrast visual mode</li>
            <li>Keyboard and touch input</li>
        </ul>
    </div>
</div>

<!-- Screen 2: Mode Selection -->
<div id="screen-mode" class="screen">
    <h1>Mode Selection — Space LEAF Corp</h1>
    <p>Select your operating mode, then press ENTER.</p>

    <select id="mode-select" onchange="onModeChange()">
        <option value="">-- Choose a mode --</option>
        <option value="firep">🔥 FIRE‑P Only (Environmental Pulse)</option>
        <option value="spaceleaf">🌿 Space LEAF Corp OS Shell</option>
        <option value="apple">🍎 Apple Companion Mode (Mock)</option>
        <option value="microsoft">🪟 Microsoft Companion Mode (Mock)</option>
        <option value="tesla">🚗 Tesla Companion Mode (Mock)</option>
    </select>
    <br>
    <button id="enter-btn" onclick="enterShell()">ENTER</button>

    <div class="panel">
        <h2>Accessibility Enabled</h2>
        <ul>
            <li>Audio narration ON (if supported)</li>
            <li>Text instructions ON</li>
            <li>High‑contrast mode ON</li>
        </ul>
    </div>
</div>

<!-- Screen 3: OS Shell Mock -->
<div id="screen-shell" class="screen">
    <h1 id="shell-title">Space LEAF Corp Shell</h1>
    <div id="shell-body" class="panel">
        <p>Shell content will appear here based on selected mode.</p>
    </div>
    <button onclick="showScreen('screen-mode')">Back to Mode Selection</button>
</div>

</body>
</html>
"""

class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        self.wfile.write(HTML.encode("utf-8"))

def start_server():
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        httpd.serve_forever()

if __name__ == "__main__":
    threading.Thread(target=start_server, daemon=True).start()
    webbrowser.open(f"http://localhost:{PORT}")
    print("FIRE‑P / Space LEAF boot server running on port", PORT)
