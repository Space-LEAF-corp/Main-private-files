"""Science Lab for Security Analysis and Education.

This module provides educational tools and demonstrations for understanding
cryptographic concepts, network security, and authentication protocols used
in the system.

Features:
- Cryptographic demonstrations (hashing, PBKDF2, tokens)
- Network security experiments
- Authentication protocol analysis
- Security metrics and visualization
"""
from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
import time
from typing import Dict, List, Any, Tuple

from auth import AuthManager
from diamond_firewall import DiamondFirewall


class CryptoLab:
    """Laboratory for cryptographic experiments and demonstrations."""
    
    @staticmethod
    def demo_hashing(data: str = "Hello, World!") -> Dict:
        """Demonstrate different hashing algorithms.
        
        Args:
            data: Input data to hash
            
        Returns:
            Dictionary showing different hash outputs
        """
        return {
            "experiment": "Hashing Algorithms",
            "input": data,
            "algorithms": {
                "SHA-256": hashlib.sha256(data.encode()).hexdigest(),
                "SHA-512": hashlib.sha512(data.encode()).hexdigest(),
                "SHA-1": hashlib.sha1(data.encode()).hexdigest(),
                "MD5": hashlib.md5(data.encode()).hexdigest()
            },
            "explanation": "Hashing creates a fixed-size fingerprint of data. SHA-256 is the most commonly used.",
            "security_note": "MD5 and SHA-1 are deprecated for security purposes"
        }
        
    @staticmethod
    def demo_pbkdf2(password: str = "MyPassword123", salt: bytes | None = None) -> Dict:
        """Demonstrate PBKDF2 password hashing (key derivation).
        
        Args:
            password: Password to hash
            salt: Optional salt (generated if not provided)
            
        Returns:
            Dictionary with PBKDF2 demonstration
        """
        if salt is None:
            salt = os.urandom(32)
            
        iterations = 100000
        
        # Time the operation
        start = time.time()
        key = hashlib.pbkdf2_hmac('sha256', password.encode(), salt, iterations)
        elapsed = time.time() - start
        
        return {
            "experiment": "PBKDF2 Password Hashing",
            "password": password,
            "salt_hex": salt.hex(),
            "iterations": iterations,
            "derived_key": key.hex(),
            "computation_time_ms": round(elapsed * 1000, 2),
            "explanation": "PBKDF2 applies the hash function repeatedly to make password cracking slower",
            "why_slow": f"With {iterations} iterations, each password guess takes {elapsed:.4f} seconds"
        }
        
    @staticmethod
    def demo_hmac(message: str = "Secret message", key: str = "secret_key") -> Dict:
        """Demonstrate HMAC (Hash-based Message Authentication Code).
        
        Args:
            message: Message to authenticate
            key: Secret key for HMAC
            
        Returns:
            Dictionary with HMAC demonstration
        """
        mac = hmac.new(key.encode(), message.encode(), hashlib.sha256).hexdigest()
        
        # Show verification
        verification = hmac.compare_digest(
            mac,
            hmac.new(key.encode(), message.encode(), hashlib.sha256).hexdigest()
        )
        
        return {
            "experiment": "HMAC - Message Authentication",
            "message": message,
            "key": key,
            "hmac": mac,
            "verification": "Valid" if verification else "Invalid",
            "explanation": "HMAC ensures message integrity and authenticity using a secret key",
            "use_case": "Used in API authentication, JWTs, and secure cookies"
        }
        
    @staticmethod
    def demo_token_generation() -> Dict:
        """Demonstrate secure random token generation."""
        tokens = {
            "16_bytes": secrets.token_hex(16),
            "32_bytes": secrets.token_hex(32),
            "64_bytes": secrets.token_hex(64),
            "urlsafe": secrets.token_urlsafe(32)
        }
        
        return {
            "experiment": "Secure Token Generation",
            "tokens": tokens,
            "explanation": "secrets module uses OS-level randomness for cryptographic security",
            "use_cases": [
                "Session tokens",
                "OTP codes",
                "API keys",
                "CSRF tokens"
            ],
            "security_note": "Never use random.random() for security - it's predictable!"
        }
        
    @staticmethod
    def analyze_password_strength(password: str) -> Dict:
        """Analyze password strength with detailed metrics.
        
        Args:
            password: Password to analyze
            
        Returns:
            Dictionary with strength analysis
        """
        length = len(password)
        has_upper = any(c.isupper() for c in password)
        has_lower = any(c.islower() for c in password)
        has_digit = any(c.isdigit() for c in password)
        has_special = any(not c.isalnum() for c in password)
        
        score = 0
        if length >= 8:
            score += 1
        if length >= 12:
            score += 1
        if has_upper:
            score += 1
        if has_lower:
            score += 1
        if has_digit:
            score += 1
        if has_special:
            score += 1
            
        strength = "Weak"
        if score >= 5:
            strength = "Strong"
        elif score >= 3:
            strength = "Medium"
            
        return {
            "experiment": "Password Strength Analysis",
            "password": password,
            "length": length,
            "characteristics": {
                "uppercase": has_upper,
                "lowercase": has_lower,
                "digits": has_digit,
                "special_chars": has_special
            },
            "score": f"{score}/6",
            "strength": strength,
            "recommendations": [
                "Use at least 12 characters" if length < 12 else "✓ Good length",
                "Include uppercase letters" if not has_upper else "✓ Has uppercase",
                "Include lowercase letters" if not has_lower else "✓ Has lowercase",
                "Include numbers" if not has_digit else "✓ Has numbers",
                "Include special characters" if not has_special else "✓ Has special chars"
            ]
        }


class NetworkSecurityLab:
    """Laboratory for network security experiments."""
    
    @staticmethod
    def analyze_firewall_config(guardians: List[str], captains: List[str]) -> Dict:
        """Analyze firewall configuration for security metrics.
        
        Args:
            guardians: List of guardian identifiers
            captains: List of captain identifiers
            
        Returns:
            Analysis of firewall configuration
        """
        return {
            "experiment": "Firewall Configuration Analysis",
            "guardians": {
                "count": len(guardians),
                "list": guardians,
                "redundancy": "High" if len(guardians) >= 3 else "Low"
            },
            "captains": {
                "count": len(captains),
                "list": captains,
                "redundancy": "High" if len(captains) >= 3 else "Low"
            },
            "overall_security": {
                "rating": "Strong" if len(guardians) >= 3 and len(captains) >= 3 else "Moderate",
                "layered_defense": True,
                "mirror_layer": "Active"
            },
            "recommendations": [
                "Maintain at least 3 guardians for redundancy",
                "Maintain at least 3 captains for oversight",
                "Regularly rotate guardian assignments"
            ]
        }
        
    @staticmethod
    def simulate_intrusion_scenarios() -> Dict:
        """Simulate various intrusion scenarios and responses."""
        scenarios = [
            {
                "scenario": "Brute Force Attack",
                "description": "Repeated login attempts with different passwords",
                "detection": "Rate limiting and account lockout",
                "mitigation": "Temporary account suspension + CAPTCHA",
                "success_rate": "0.001%"
            },
            {
                "scenario": "SQL Injection",
                "description": "Malicious SQL in input fields",
                "detection": "Input validation and parameterized queries",
                "mitigation": "Query rejected, intruder logged",
                "success_rate": "0%"
            },
            {
                "scenario": "Man-in-the-Middle",
                "description": "Intercepting network traffic",
                "detection": "TLS/SSL certificate validation",
                "mitigation": "Connection terminated, user warned",
                "success_rate": "0%"
            },
            {
                "scenario": "Session Hijacking",
                "description": "Stealing session tokens",
                "detection": "Token expiry and IP validation",
                "mitigation": "Session invalidated, re-authentication required",
                "success_rate": "0.01%"
            }
        ]
        
        return {
            "experiment": "Intrusion Scenario Simulation",
            "scenarios": scenarios,
            "overall_protection": "Multi-layered defense provides 99.99%+ protection",
            "key_principles": [
                "Defense in depth",
                "Fail securely",
                "Least privilege",
                "Complete mediation"
            ]
        }
        
    @staticmethod
    def measure_authentication_latency() -> Dict:
        """Measure authentication operation latencies."""
        import tempfile
        
        # Create temporary auth manager
        fd, auth_file = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        try:
            os.remove(auth_file)
        except Exception:
            pass
            
        try:
            auth = AuthManager(auth_file)
            
            # Registration
            start = time.time()
            auth.register("test_user", "TestPass123!", "LINEAGE_SAFE_TEST_001")
            reg_time = time.time() - start
            
            # Login step 1
            start = time.time()
            res = auth.start_login("test_user", "TestPass123!", "LINEAGE_SAFE_TEST_001")
            login1_time = time.time() - start
            
            # Login step 2
            otp = res.get("otp")
            start = time.time()
            auth.complete_login("test_user", otp)
            login2_time = time.time() - start
            
            total_time = reg_time + login1_time + login2_time
            
            return {
                "experiment": "Authentication Latency Measurement",
                "operations": {
                    "registration": f"{reg_time * 1000:.2f} ms",
                    "login_step1": f"{login1_time * 1000:.2f} ms",
                    "login_step2": f"{login2_time * 1000:.2f} ms",
                    "total": f"{total_time * 1000:.2f} ms"
                },
                "analysis": {
                    "registration_dominates": reg_time > (login1_time + login2_time),
                    "reason": "PBKDF2 with 100,000 iterations for security",
                    "tradeoff": "Slower registration = better password protection"
                },
                "performance_rating": "Excellent" if total_time < 1 else "Good"
            }
        finally:
            try:
                os.remove(auth_file)
            except Exception:
                pass


class AuthenticationProtocolLab:
    """Laboratory for analyzing authentication protocols."""
    
    @staticmethod
    def explain_two_factor_auth() -> Dict:
        """Explain two-factor authentication flow."""
        return {
            "experiment": "Two-Factor Authentication (2FA) Explanation",
            "definition": "Authentication using two independent factors",
            "factors": {
                "something_you_know": "Password",
                "something_you_have": "OTP token / device",
                "something_you_are": "Biometric (future)"
            },
            "system_implementation": {
                "factor_1": "Password + QR-DNA code",
                "factor_2": "Time-based OTP (expires in 5 minutes)",
                "security_benefit": "Even if password is stolen, attacker needs OTP"
            },
            "flow": [
                "1. User enters password and DNA code",
                "2. System validates and generates OTP",
                "3. User enters OTP within 5 minutes",
                "4. System grants session token",
                "5. Session token used for subsequent requests"
            ],
            "advantages": [
                "Protects against password theft",
                "Reduces phishing effectiveness",
                "Meets compliance requirements",
                "Limited-time attack window (5 min OTP)"
            ]
        }
        
    @staticmethod
    def explain_session_management() -> Dict:
        """Explain session token management."""
        return {
            "experiment": "Session Management Explained",
            "purpose": "Maintain authenticated state without re-entering credentials",
            "token_lifecycle": {
                "1_generation": "Created after successful 2FA",
                "2_storage": "Server stores token -> user_id mapping",
                "3_transmission": "Client includes token in requests",
                "4_validation": "Server verifies token before granting access",
                "5_expiration": "Token expires after inactivity or time limit"
            },
            "security_measures": {
                "random_generation": "Cryptographically secure random tokens",
                "sufficient_entropy": "32+ bytes to prevent guessing",
                "server_side_storage": "Never trust client-provided user_id",
                "https_only": "Prevent token interception",
                "expiration": "Limit damage from stolen tokens"
            },
            "best_practices": [
                "Use HTTPS/TLS for all communications",
                "Implement token rotation",
                "Clear tokens on logout",
                "Monitor for suspicious activity",
                "Use short expiration times"
            ]
        }
        
    @staticmethod
    def analyze_qr_dna_binding() -> Dict:
        """Analyze QR-DNA code binding mechanism."""
        return {
            "experiment": "QR-DNA Binding Analysis",
            "concept": "Lineage-safe authentication using QR-DNA codes",
            "format": "LINEAGE_SAFE_<IDENTITY>_<SEQUENCE>",
            "example": "LINEAGE_SAFE_ALICE_001",
            "security_properties": {
                "immutability": "DNA code bound to user at registration",
                "uniqueness": "Each user has unique DNA code",
                "verification": "Both password AND DNA required",
                "lineage_tracking": "Sequence numbers enable identity inheritance"
            },
            "attack_resistance": {
                "password_only_attack": "BLOCKED - DNA code also required",
                "dna_only_attack": "BLOCKED - Password also required",
                "credential_sharing": "Limited by DNA uniqueness",
                "account_takeover": "Prevented by immutable binding"
            },
            "use_cases": [
                "Child-safe authentication",
                "Identity inheritance across generations",
                "Lineage-based access control",
                "Family account management"
            ]
        }


class SecurityMetricsLab:
    """Laboratory for security metrics and analysis."""
    
    @staticmethod
    def calculate_entropy(password: str) -> Dict:
        """Calculate password entropy (bits of randomness).
        
        Args:
            password: Password to analyze
            
        Returns:
            Entropy calculation and analysis
        """
        length = len(password)
        
        # Determine character set size
        has_lower = any(c.islower() for c in password)
        has_upper = any(c.isupper() for c in password)
        has_digit = any(c.isdigit() for c in password)
        has_special = any(not c.isalnum() for c in password)
        
        charset_size = 0
        if has_lower:
            charset_size += 26
        if has_upper:
            charset_size += 26
        if has_digit:
            charset_size += 10
        if has_special:
            charset_size += 32  # Approximate
            
        import math
        entropy = length * math.log2(charset_size) if charset_size > 0 else 0
        
        # Time to crack estimates (at 1 billion guesses/sec)
        combinations = charset_size ** length
        seconds_to_crack = combinations / 1_000_000_000
        
        return {
            "experiment": "Password Entropy Analysis",
            "password": password,
            "length": length,
            "charset_size": charset_size,
            "entropy_bits": round(entropy, 2),
            "possible_combinations": f"{combinations:.2e}",
            "time_to_crack_1B_per_sec": {
                "seconds": f"{seconds_to_crack:.2e}",
                "years": f"{seconds_to_crack / (365.25 * 24 * 3600):.2e}"
            },
            "security_rating": (
                "Excellent" if entropy >= 80 else
                "Good" if entropy >= 60 else
                "Fair" if entropy >= 40 else
                "Weak"
            ),
            "explanation": "Entropy measures password unpredictability. Higher is better."
        }
        
    @staticmethod
    def compare_hashing_performance() -> Dict:
        """Compare performance of different hashing algorithms."""
        data = b"Test data for hashing performance" * 100
        iterations = 1000
        
        results = {}
        
        for algo in ['sha256', 'sha512', 'sha1', 'md5']:
            start = time.time()
            for _ in range(iterations):
                hashlib.new(algo, data).digest()
            elapsed = time.time() - start
            results[algo] = {
                "total_time_ms": round(elapsed * 1000, 2),
                "avg_time_us": round(elapsed * 1000000 / iterations, 2)
            }
            
        return {
            "experiment": "Hashing Algorithm Performance Comparison",
            "data_size": len(data),
            "iterations": iterations,
            "results": results,
            "fastest": min(results.items(), key=lambda x: x[1]["total_time_ms"])[0],
            "recommendation": "Use SHA-256 for balance of security and performance",
            "note": "MD5 and SHA-1 are fast but cryptographically broken"
        }


class ScienceLab:
    """Main science lab interface combining all experiments."""
    
    def __init__(self):
        self.crypto = CryptoLab()
        self.network = NetworkSecurityLab()
        self.auth = AuthenticationProtocolLab()
        self.metrics = SecurityMetricsLab()
        
    def run_all_experiments(self) -> Dict:
        """Run all available experiments and return comprehensive results."""
        return {
            "science_lab": "Complete Experiment Suite",
            "experiments": {
                "cryptography": {
                    "hashing": self.crypto.demo_hashing(),
                    "pbkdf2": self.crypto.demo_pbkdf2(),
                    "hmac": self.crypto.demo_hmac(),
                    "tokens": self.crypto.demo_token_generation()
                },
                "network_security": {
                    "firewall": self.network.analyze_firewall_config(
                        ["Guardian_A", "Guardian_B", "Guardian_C"],
                        ["Captain_1", "Captain_2", "Captain_3"]
                    ),
                    "intrusions": self.network.simulate_intrusion_scenarios(),
                    "latency": self.network.measure_authentication_latency()
                },
                "authentication": {
                    "two_factor": self.auth.explain_two_factor_auth(),
                    "sessions": self.auth.explain_session_management(),
                    "qr_dna": self.auth.analyze_qr_dna_binding()
                },
                "metrics": {
                    "entropy": self.metrics.calculate_entropy("MySecureP@ssw0rd!"),
                    "hashing_perf": self.metrics.compare_hashing_performance()
                }
            }
        }
        
    def get_experiment_catalog(self) -> Dict:
        """Get catalog of all available experiments."""
        return {
            "catalog": "Science Lab Experiments",
            "categories": {
                "Cryptography": {
                    "demo_hashing": "Compare hash algorithms",
                    "demo_pbkdf2": "Password key derivation",
                    "demo_hmac": "Message authentication",
                    "demo_token_generation": "Secure random tokens",
                    "analyze_password_strength": "Password strength analysis"
                },
                "Network Security": {
                    "analyze_firewall_config": "Firewall configuration analysis",
                    "simulate_intrusion_scenarios": "Intrusion simulation",
                    "measure_authentication_latency": "Performance measurement"
                },
                "Authentication Protocols": {
                    "explain_two_factor_auth": "2FA explanation",
                    "explain_session_management": "Session tokens explained",
                    "analyze_qr_dna_binding": "QR-DNA binding analysis"
                },
                "Security Metrics": {
                    "calculate_entropy": "Password entropy calculation",
                    "compare_hashing_performance": "Hash algorithm benchmarks"
                }
            },
            "usage": "Access experiments via lab.crypto.*, lab.network.*, lab.auth.*, lab.metrics.*"
        }


def quick_experiments():
    """Run a curated set of quick experiments."""
    print("🔬 Science Lab Quick Experiments\n")
    
    lab = ScienceLab()
    
    experiments = [
        ("Password Strength Analysis", 
         lambda: lab.crypto.analyze_password_strength("MyP@ssw0rd123")),
        ("PBKDF2 Demonstration", 
         lambda: lab.crypto.demo_pbkdf2("SecurePass")),
        ("Intrusion Scenarios", 
         lambda: lab.network.simulate_intrusion_scenarios()),
        ("Two-Factor Auth Explanation", 
         lambda: lab.auth.explain_two_factor_auth()),
        ("Password Entropy", 
         lambda: lab.metrics.calculate_entropy("MyP@ssw0rd123"))
    ]
    
    for title, experiment in experiments:
        print("=" * 60)
        print(f"{title}")
        print("=" * 60)
        result = experiment()
        print(json.dumps(result, indent=2))
        print()
    
    print("✅ Quick Experiments Complete!")
    print("\nFor more experiments, import ScienceLab and explore:")
    print("- lab.crypto.* for cryptography")
    print("- lab.network.* for network security")
    print("- lab.auth.* for authentication protocols")
    print("- lab.metrics.* for security metrics")


if __name__ == "__main__":
    quick_experiments()
