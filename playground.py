"""Interactive Playground for Security Features.

This module provides a safe, interactive environment for experimenting with
the authentication, firewall, and MJ protocol features. Perfect for learning,
testing, and prototyping security flows.

Features:
- Interactive authentication demos
- Firewall testing sandbox
- MJ protocol experimentation
- Guided tutorials for new users
"""
from __future__ import annotations

import json
import os
import tempfile
from typing import Dict, List, Any

from auth import AuthManager
from diamond_firewall import DiamondFirewall
from secured_firewall import SecuredFirewall
from mj_protocol import MJLocalLayer, CeremonialsManager


class PlaygroundSession:
    """Represents a playground session with isolated state."""
    
    def __init__(self, session_id: str):
        self.session_id = session_id
        self.users = {}
        self.firewall_logs = []
        self.commands_history = []
        
    def log_command(self, command: str, result: Any):
        """Log a command and its result."""
        self.commands_history.append({
            "command": command,
            "result": result
        })
        
    def get_history(self) -> List[Dict]:
        """Get command history."""
        return self.commands_history


class Playground:
    """Interactive playground for security features experimentation.
    
    This class provides a sandboxed environment where users can:
    - Test authentication flows
    - Experiment with firewall rules
    - Learn about MJ protocol
    - Run security demonstrations
    """
    
    def __init__(self, isolated: bool = True):
        """Initialize playground.
        
        Args:
            isolated: If True, uses temporary files for complete isolation
        """
        self.isolated = isolated
        if isolated:
            fd, self.auth_file = tempfile.mkstemp(suffix=".json", prefix="playground_")
            os.close(fd)
            try:
                os.remove(self.auth_file)
            except Exception:
                pass
        else:
            self.auth_file = "playground_users.json"
            
        self.secured_firewall = SecuredFirewall(users_file=self.auth_file)
        self.sessions = {}
        self.tutorial_step = 0
        
    def cleanup(self):
        """Clean up temporary files."""
        if self.isolated:
            try:
                os.remove(self.auth_file)
            except Exception:
                pass
                
    def create_session(self, session_id: str) -> PlaygroundSession:
        """Create a new playground session."""
        session = PlaygroundSession(session_id)
        self.sessions[session_id] = session
        return session
        
    def get_session(self, session_id: str) -> PlaygroundSession | None:
        """Get existing session."""
        return self.sessions.get(session_id)
        
    # ==================== Authentication Demos ====================
    
    def demo_registration(self, user_id: str = "demo_user") -> Dict:
        """Demo user registration with guided explanation.
        
        Args:
            user_id: User identifier for the demo
            
        Returns:
            Dictionary with demo results and explanations
        """
        password = "SecureDemo123!"
        dna_code = f"LINEAGE_SAFE_{user_id.upper()}_001"
        
        result = self.secured_firewall.register_user(user_id, password, dna_code)
        
        return {
            "demo": "registration",
            "explanation": "User registration requires: user_id, strong password (8+ chars), and QR-DNA code",
            "user_id": user_id,
            "password": password,
            "dna_code": dna_code,
            "result": result,
            "next_step": "Try demo_login() to authenticate"
        }
        
    def demo_login(self, user_id: str = "demo_user") -> Dict:
        """Demo full login flow with OTP challenge.
        
        Args:
            user_id: User identifier (must be registered first)
            
        Returns:
            Dictionary with login flow results
        """
        password = "SecureDemo123!"
        dna_code = f"LINEAGE_SAFE_{user_id.upper()}_001"
        
        # Step 1: Initial login
        step1 = self.secured_firewall.login_step1(user_id, password, dna_code)
        
        if step1.get("status") != "challenge":
            return {
                "demo": "login",
                "error": "Login step 1 failed",
                "details": step1,
                "hint": "Did you run demo_registration() first?"
            }
            
        otp = step1.get("otp_token")
        
        # Step 2: Complete with OTP
        step2 = self.secured_firewall.login_step2(user_id, otp)
        
        return {
            "demo": "login",
            "explanation": "Login requires 2 steps: password+DNA verification, then OTP challenge",
            "step1_status": step1.get("status"),
            "otp_received": otp[:8] + "..." if otp else None,
            "step2_status": step2.get("status"),
            "session_token": step2.get("session_token", "N/A"),
            "next_step": "Try demo_firewall_access() to access the firewall"
        }
        
    def demo_firewall_access(self, user_id: str = "demo_user") -> Dict:
        """Demo firewall access after successful authentication.
        
        Args:
            user_id: User identifier (must be logged in)
            
        Returns:
            Dictionary with firewall access results
        """
        # Get session token by doing a fresh login
        password = "SecureDemo123!"
        dna_code = f"LINEAGE_SAFE_{user_id.upper()}_001"
        
        step1 = self.secured_firewall.login_step1(user_id, password, dna_code)
        if step1.get("status") != "challenge":
            return {"error": "Login failed", "details": step1}
            
        otp = step1.get("otp_token")
        step2 = self.secured_firewall.login_step2(user_id, otp)
        
        if step2.get("status") != "ok":
            return {"error": "OTP verification failed", "details": step2}
            
        session_token = step2.get("session_token")
        
        # Access firewall
        access = self.secured_firewall.access_firewall(session_token)
        
        return {
            "demo": "firewall_access",
            "explanation": "Firewall access requires valid session token from successful login",
            "session_token": session_token[:8] + "..." if session_token else None,
            "access_result": access,
            "firewall_status": access.get("status"),
            "guardians": access.get("guardians_notified", []),
            "next_step": "Try demo_intrusion_detection() to see how intrusions are handled"
        }
        
    def demo_intrusion_detection(self) -> Dict:
        """Demo intrusion detection and mirror layer trapping."""
        result = self.secured_firewall.intrude_attempt(
            "malicious_actor",
            "unauthorized_payload"
        )
        
        return {
            "demo": "intrusion_detection",
            "explanation": "Intrusions are automatically detected and trapped in mirror layer",
            "intruder": "malicious_actor",
            "result": result,
            "protection": "Diamond Firewall with Mirror Layer Protection"
        }
        
    # ==================== Guided Tutorial ====================
    
    def start_tutorial(self) -> Dict:
        """Start guided tutorial from step 1."""
        self.tutorial_step = 1
        return {
            "tutorial": "Welcome to Security Playground!",
            "step": 1,
            "instruction": "Let's start by registering a user. Run: tutorial_next()",
            "overview": [
                "Step 1: User Registration",
                "Step 2: Multi-layer Login",
                "Step 3: Firewall Access",
                "Step 4: Intrusion Detection"
            ]
        }
        
    def tutorial_next(self) -> Dict:
        """Execute next tutorial step."""
        if self.tutorial_step == 0:
            return {"error": "Tutorial not started. Run start_tutorial() first"}
            
        if self.tutorial_step == 1:
            result = self.demo_registration("tutorial_user")
            self.tutorial_step = 2
            result["tutorial_step"] = 2
            result["next_instruction"] = "Great! Now run tutorial_next() again for login"
            return result
            
        elif self.tutorial_step == 2:
            result = self.demo_login("tutorial_user")
            self.tutorial_step = 3
            result["tutorial_step"] = 3
            result["next_instruction"] = "Excellent! Run tutorial_next() for firewall access"
            return result
            
        elif self.tutorial_step == 3:
            result = self.demo_firewall_access("tutorial_user")
            self.tutorial_step = 4
            result["tutorial_step"] = 4
            result["next_instruction"] = "Perfect! Run tutorial_next() to see intrusion detection"
            return result
            
        elif self.tutorial_step == 4:
            result = self.demo_intrusion_detection()
            self.tutorial_step = 0
            result["tutorial_step"] = "Complete!"
            result["congratulations"] = "You've completed the security tutorial!"
            result["next_steps"] = [
                "Try custom_experiment() for your own tests",
                "Explore science_lab.py for deeper security analysis"
            ]
            return result
            
        return {"error": "Tutorial already complete"}
        
    # ==================== Custom Experiments ====================
    
    def custom_experiment(
        self,
        user_id: str,
        password: str,
        dna_code: str,
        test_intrusion: bool = False
    ) -> Dict:
        """Run custom security experiment.
        
        Args:
            user_id: Custom user identifier
            password: Custom password to test
            dna_code: Custom DNA code (must start with LINEAGE_SAFE_)
            test_intrusion: Whether to simulate intrusion attempt
            
        Returns:
            Experiment results
        """
        results = {"experiment": "custom", "steps": []}
        
        # Register
        reg = self.secured_firewall.register_user(user_id, password, dna_code)
        results["steps"].append({"action": "register", "result": reg})
        
        if reg.get("status") != "ok":
            return results
            
        # Login
        step1 = self.secured_firewall.login_step1(user_id, password, dna_code)
        results["steps"].append({"action": "login_step1", "result": step1})
        
        if step1.get("status") != "challenge":
            return results
            
        otp = step1.get("otp_token")
        step2 = self.secured_firewall.login_step2(user_id, otp)
        results["steps"].append({"action": "login_step2", "result": step2})
        
        if step2.get("status") != "ok":
            return results
            
        # Access
        session = step2.get("session_token")
        access = self.secured_firewall.access_firewall(session)
        results["steps"].append({"action": "firewall_access", "result": access})
        
        # Optional intrusion test
        if test_intrusion:
            intrusion = self.demo_intrusion_detection()
            results["steps"].append({"action": "intrusion_test", "result": intrusion})
            
        results["summary"] = "Experiment completed successfully"
        return results
        
    def get_help(self) -> Dict:
        """Get help information about playground features."""
        return {
            "playground_help": {
                "demos": {
                    "demo_registration()": "Test user registration",
                    "demo_login()": "Test multi-layer login",
                    "demo_firewall_access()": "Test firewall access",
                    "demo_intrusion_detection()": "Test intrusion handling"
                },
                "tutorial": {
                    "start_tutorial()": "Begin guided tutorial",
                    "tutorial_next()": "Advance to next tutorial step"
                },
                "custom": {
                    "custom_experiment()": "Run custom security experiments"
                },
                "utility": {
                    "cleanup()": "Clean up temporary files",
                    "get_help()": "Show this help message"
                }
            },
            "tip": "Start with start_tutorial() for a guided experience!"
        }


def quick_start():
    """Quick start demo - runs all basic demos."""
    print("🎮 Security Playground Quick Start\n")
    
    playground = Playground(isolated=True)
    
    try:
        print("=" * 60)
        print("DEMO 1: User Registration")
        print("=" * 60)
        result = playground.demo_registration("quickstart_user")
        print(json.dumps(result, indent=2))
        
        print("\n" + "=" * 60)
        print("DEMO 2: Multi-Layer Login")
        print("=" * 60)
        result = playground.demo_login("quickstart_user")
        print(json.dumps(result, indent=2))
        
        print("\n" + "=" * 60)
        print("DEMO 3: Firewall Access")
        print("=" * 60)
        result = playground.demo_firewall_access("quickstart_user")
        print(json.dumps(result, indent=2))
        
        print("\n" + "=" * 60)
        print("DEMO 4: Intrusion Detection")
        print("=" * 60)
        result = playground.demo_intrusion_detection()
        print(json.dumps(result, indent=2))
        
        print("\n" + "=" * 60)
        print("✅ Quick Start Complete!")
        print("=" * 60)
        print("\nNext steps:")
        print("- Import Playground and run start_tutorial()")
        print("- Try custom_experiment() with your own parameters")
        print("- Explore science_lab.py for security analysis")
        
    finally:
        playground.cleanup()


if __name__ == "__main__":
    quick_start()
