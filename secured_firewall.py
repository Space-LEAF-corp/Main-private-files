"""Example: diamond_firewall integrated with multi-layer authentication.

This module shows how to use AuthManager with DiamondFirewall to:
1. Register a user with QR-DNA + password + user_id
2. Login with multi-layer auth (password + dna + otp)
3. Grant access to the firewall only after successful auth
"""
from auth import AuthManager
from diamond_firewall import DiamondFirewall


class SecuredFirewall:
    """Wraps DiamondFirewall with authentication."""

    def __init__(self, users_file: str = "auth_users.json"):
        self.auth_manager = AuthManager(users_file)
        # Initialize firewall with demo guardians and captains
        self.firewall = DiamondFirewall(
            guardians=["Miko", "JD", "Guardian_A", "Guardian_B"],
            captains=["Captain_1", "Captain_2", "Captain_3", "Captain_4"],
            un_consent=True,
        )
        self._authenticated_sessions = {}  # session_token -> user_id

    def register_user(self, user_id: str, password: str, dna_code: str) -> dict:
        """Register a new user for first sign-in."""
        try:
            self.auth_manager.register(user_id, password, dna_code)
            return {"status": "ok", "message": f"User {user_id} registered successfully"}
        except ValueError as e:
            return {"status": "error", "message": str(e)}

    def login_step1(self, user_id: str, password: str, dna_code: str) -> dict:
        """Initial login with password and dna_code. Returns OTP challenge."""
        res = self.auth_manager.start_login(user_id, password, dna_code)
        if res.get("status") == "challenge":
            return {"status": "challenge", "otp_token": res.get("otp")}
        return res

    def login_step2(self, user_id: str, otp: str) -> dict:
        """Complete login with OTP. Returns session token if successful."""
        res = self.auth_manager.complete_login(user_id, otp)
        if res.get("status") == "ok":
            session_token = res.get("session_token")
            self._authenticated_sessions[session_token] = user_id
            return {"status": "ok", "session_token": session_token}
        return res

    def access_firewall(self, session_token: str) -> dict:
        """Attempt to access the firewall. Requires valid session token."""
        user_id = self._authenticated_sessions.get(session_token)
        if not user_id:
            return {"status": "denied", "reason": "invalid or expired session"}

        # Try to access the firewall with a valid dna code
        try:
            result = self.firewall.access_request("LINEAGE_SAFE_QR123")
            return {"status": "ok", "message": result}
        except Exception as e:
            return {"status": "denied", "reason": str(e)}

    def intrude_attempt(self, attacker_id: str, payload: str) -> dict:
        """Simulate an intrusion attempt (firewall traps attacker)."""
        result = self.firewall.intrusion_attempt(attacker_id, payload)
        return {"status": "trapped", "message": result}
