import qrcode
import hashlib
import uuid
from flask import Flask, request, jsonify

app = Flask(__name__)

# In-memory store (replace with a database in production)
valid_tokens = {}

def generate_qr(user_id: str):
    # Create unique token
    raw = f"{user_id}-{uuid.uuid4()}"
    token = hashlib.sha256(raw.encode()).hexdigest()
    valid_tokens[token] = True
    
    # Generate QR
    qr = qrcode.make(token)
    qr.save(f"{user_id}_qr.png")
    return token

@app.route("/verify", methods=["POST"])
def verify_qr():
    token = request.json.get("token")
    if token in valid_tokens:
        return jsonify({"status": "valid", "access": "granted"})
    return jsonify({"status": "invalid", "access": "denied"}), 403

if __name__ == "__main__":
    # Example: generate QR for user123
    print("Generated token:", generate_qr("user123"))
    app.run(port=5000)
