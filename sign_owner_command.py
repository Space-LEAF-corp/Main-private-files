# sign_owner_command.py
import time, hashlib, hmac

def hmac_sha256(secret: str, message: str) -> str:
    return hmac.new(secret.encode("utf-8"), message.encode("utf-8"), hashlib.sha256).hexdigest()

def make_signature(owner_id: str, secret: str, payload: str) -> tuple:
    ts = int(time.time())
    msg = f"{owner_id}:{ts}:{payload}"
    return hmac_sha256(secret, msg), ts

if __name__ == "__main__":
    owner_id = "leif.w.sogge"
    secret = "CHANGE_ME_TO_A_STRONG_SECRET"
    payload = "set_lockdown:False"
    sig, ts = make_signature(owner_id, secret, payload)
    print({"signature": sig, "timestamp": ts})
