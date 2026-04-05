import ui
import http.server
import socketserver
import webbrowser
import threading

PORT = 8080

HTML = """
<!DOCTYPE html>
<html>
<head>
<title>FIRE‐P Boot</title>
<style>
body {
    font-family: Arial, sans-serif;
    background: #000;
    color: #eee;
    padding: 40px;
}
h1 { color: #ff6600; }
input {
    padding: 10px;
    font-size: 18px;
    width: 300px;
}
button {
    padding: 10px 20px;
    font-size: 18px;
    margin-top: 10px;
}
#instructions {
    margin-top: 30px;
    padding: 20px;
    background: #111;
    border-radius: 8px;
    display: none;
}
</style>
<script>
function startSystem() {
    const val = document.getElementById("code").value.trim().toUpperCase();
    if (val === "FIRE-P" || val === "FIREP") {
        document.getElementById("instructions").style.display = "block";
        const audio = new SpeechSynthesisUtterance(
            "Welcome to Space LEAF Corp. You have activated FIRE P. " +
            "This startup mode ensures clarity, safety, and accessibility for all users. " +
            "Follow the on screen instructions to continue."
        );
        speechSynthesis.speak(audio);
    } else {
        alert("Invalid startup code. Type FIRE-P to begin.");
    }
}
</script>
</head>
<body>
<h1>🔥 FIRE‐P Startup</h1>
<p>Type <strong>FIRE‐P</strong> to begin system initialization.</p>

<input id="code" placeholder="Enter FIRE-P" />
<button onclick="startSystem()">Start</button>

<div id="instructions">
<h2>FIRE‐P Startup Mode Activated</h2>
<p>Accessibility Enabled:</p>
<ul>
    <li>Audio narration ON</li>
    <li>Text instructions ON</li>
    <li>High‐contrast mode ON</li>
    <li>Gesture + keyboard input accepted</li>
</ul>
<p>Press ENTER to proceed to login setup.</p>
</div>

</body>
</html>
"""

class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        self.wfile.write(HTML.encode())

def start_server():
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        httpd.serve_forever()

threading.Thread(target=start_server, daemon=True).start()

webbrowser.open(f"http://localhost:{PORT}")
print("FIRE‐P boot server running on port", PORT)
v = ui.load_view()
v.present('sheet')
