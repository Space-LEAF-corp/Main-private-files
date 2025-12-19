"""Diamond Firewall — cleaned and importable version.

This file is a cleaned copy of the original `# diamond_firewall.py` (renamed
to remove the leading #). It implements a minimal lineage-safe firewall
simulation with thread-safety and simple unit tests at the bottom.
"""
from __future__ import annotations
import hashlib
import logging
import threading
import time
import os
from typing import Dict, List, Optional, Any

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DiamondFirewall")


class AccessDeniedError(Exception):
    pass


class DiamondFirewall:
    """
    Minimal lineage-safe firewall simulation:
      - mirror_layer stores attacker traps (fake hashes)
      - diamond_layer holds renewal fingerprints for self-renewal cycles
      - dna_qr_check enforces collective consent
    Thread-safe for concurrent intrusion checks.
    """

    def __init__(
        self,
        guardians: List[str],
        captains: List[str],
        un_consent: bool,
        required_captains: int = 4,
        hash_algo: str = "sha256",
    ):
        self.guardians = list(guardians)
        self.captains = list(captains)
        self.un_consent = bool(un_consent)
        self.required_captains = int(required_captains)
        self.hash_algo = hash_algo
        self.mirror_layer: Dict[str, str] = {}
        self.diamond_layer: Dict[str, str] = {}
        self.active = True
        self._lock = threading.RLock()
        # Network policies (wifi, hifi, ethernet, cellular, satellite, internet)
        # Each policy is a dict with simple enforcement flags and thresholds.
        self.network_policies: Dict[str, Dict[str, Any]] = {
            "wifi": {"require_qr": True, "require_mfa": True, "trusted_ssids": []},
            "hifi": {"require_qr": False, "require_mfa": False, "trusted_devices": []},
            "ethernet": {"require_qr": False, "require_mfa": False},
            "cellular": {"require_qr": True, "uplink_strict": True, "allowed_carriers": []},
            "satellite": {"require_qr": True, "uplink_strict": True, "allowed_satellites": []},
            "internet": {"require_qr": True},
        }

        # Global uplink control: if True, firewall enforces strict checks for any uplinked networks
        self.global_uplink_protection: bool = True
        # Lockdown state and keys (primary: Space Leaf Corp main key, secondary: UN key)
        self.lockdown_active: bool = False
        self.lockdown_keys: Dict[str, str] = {"primary": "", "secondary_un": ""}

        # Intruder warning tracking: attacker_id -> count
        self.intruder_warnings: Dict[str, int] = {}
        # Optional notifier callback for external integration (callable(attacker_id, message, step))
        self.notifier_callback = None
        # Webhook notifier (URL + optional headers)
        self.notifier_webhook_url: Optional[str] = None
        self.notifier_webhook_headers: Dict[str, str] = {}
        # If UN activates lockdown, it will be limited flag
        self.lockdown_un_limited: bool = True
        # Recipients for official notifications (Space Leaf, UN, President of the United States)
        # Each recipient maps to a webhook URL (or None). Configure via `set_notification_recipients`.
        self.notification_recipients: Dict[str, Optional[str]] = {
            "space_leaf": None,
            "un": None,
            "president_of_the_united_states": None,
        }
        # Core admin controls: immutable core flag and core admins set
        # Leif William Sogge is a co-owner and core admin by default.
        self.immutable_core: bool = True
        self.core_admins: set[str] = {"Leif William Sogge"}
        # number of guardian approvals required for core changes (besides a core admin)
        self.required_guardian_approvals: int = max(1, min(len(self.guardians), 1))
        # Blocklist of attacker ids
        self.blocklist: set[str] = set()
        # Sponge service thread
        self._sponge_thread = None
        self._sponge_running = False

    def _make_hash(self, data: str) -> str:
        h = hashlib.new(self.hash_algo)
        h.update(data.encode("utf-8"))
        return h.hexdigest()

    def dna_qr_check(self, dna_code: str, required_captains: Optional[int] = None) -> bool:
        """
        Verify DNA QR code and collective consent.
        dna_code must begin with "LINEAGE_SAFE" to be considered valid.
        """
        required = self.required_captains if required_captains is None else required_captains
        with self._lock:
            ok = (
                dna_code.startswith("LINEAGE_SAFE")
                and len(self.captains) >= required
                and self.un_consent is True
                and len(self.guardians) > 0
            )
            logger.debug("dna_qr_check result=%s", ok)
            return ok

    def set_network_policy(self, network: str, policy: Dict[str, Any]) -> None:
        """Set or update policy for a named network type."""
        with self._lock:
            self.network_policies[network] = {**self.network_policies.get(network, {}), **policy}
            logger.info("Network policy updated for %s: %s", network, self.network_policies[network])

    def _check_network_policy(self, network: str, metadata: Optional[Dict[str, object]] = None) -> bool:
        """Evaluate network-specific protections. Returns True if allowed to proceed to dna check."""
        metadata = metadata or {}
        policy = self.network_policies.get(network, self.network_policies.get("internet", {}))

        # If policy requires QR and the caller didn't provide it, deny early (caller should perform qr check)
        if policy.get("require_qr", False) and not metadata.get("dna_code"):
            logger.warning("Access denied by network policy: dna_code required for %s", network)
            return False

        # Wi-Fi specific checks: trusted SSID list
        if network == "wifi" and policy.get("trusted_ssids"):
            ssid = metadata.get("ssid")
            trusted_ssids = policy.get("trusted_ssids", [])
            # Ensure trusted_ssids is a list, set, or tuple before using 'in'
            if ssid and isinstance(trusted_ssids, (list, set, tuple)) and ssid in trusted_ssids:
                logger.debug("Wi-Fi SSID %s is trusted", ssid)
                return True
            # otherwise require further checks (dna/mfa)

        # Cellular/satellite uplink checks
        if network in ("cellular", "satellite") and self.global_uplink_protection:
            # if strict uplink, require allowed carrier/satellite meta
            if policy.get("uplink_strict", False):
                allowed = policy.get("allowed_carriers") or policy.get("allowed_satellites")
                carrier = metadata.get("carrier") or metadata.get("satellite_id")
                # Ensure allowed is a list or set before using 'not in'
                if allowed and carrier:
                    if isinstance(allowed, (list, set, tuple)):
                        if carrier not in allowed:
                            logger.warning("Uplink denied for %s carrier %s", network, carrier)
                            return False
                    else:
                        logger.warning("Allowed carriers/satellites is not a list/set/tuple: %s", allowed)
                        return False

        return True

    def intrusion_attempt(self, attacker_id: str, payload: str) -> str:
        """
        Redirect attacker into mirror layer (phantom system).
        Returns a summary message. Trapped values are fake hashes.
        """
        # If firewall inactive, ignore
        if not self.active:
            logger.warning("Firewall inactive; intrusion attempt ignored.")
            return "Firewall inactive."

        # Issue progressive warnings first
        warning_msg = self.notify_intruder(attacker_id)
        if warning_msg:
            return warning_msg

        # After warnings exhausted, trap attacker
        salt = f"{time.time_ns()}:{os.getpid()}"
        fake_hash = self._make_hash(payload + salt)
        with self._lock:
            self.mirror_layer[attacker_id] = fake_hash
            logger.info("Attacker %s trapped in mirror layer.", attacker_id)
        return f"Attacker {attacker_id} trapped in mirror layer."

    def notify_intruder(self, attacker_id: str) -> Optional[str]:
        """Send progressive warnings to an intruder.

        Returns a warning message if a warning step is to be sent; returns None when
        warnings are exhausted and normal intrusion handling should proceed.
        """
        with self._lock:
            count = self.intruder_warnings.get(attacker_id, 0) + 1
            self.intruder_warnings[attacker_id] = count

        # Compose staged warnings
        if count == 1:
            msg = (
                f"WARNING 1 to {attacker_id}: Do you really want to do this?"
            )
        elif count == 2:
            msg = (
                f"WARNING 2 to {attacker_id}: This is illegal activity. Stop now or be recorded."
            )
        elif count == 3:
            # Final warning and escalate
            msg = (
                f"WARNING 3 to {attacker_id}: Authorities have been notified. Cease immediately."
            )
            # Escalate to authorities
            self._escalate_to_authorities(attacker_id)
        else:
            # count > 3, warnings exhausted, proceed to trap
            return None

        logger.warning(msg)
        # Optionally invoke external notifier (callable, webhook, syslog, etc.)
        try:
            if callable(self.notifier_callback):
                self.notifier_callback(attacker_id, msg, count)
        except Exception:
            logger.exception("Notifier callback failed for attacker=%s", attacker_id)

        # If a webhook is configured, send a non-blocking POST
        try:
            if self.notifier_webhook_url:
                self._send_webhook_notification(attacker_id, msg, count)
        except Exception:
            logger.exception("Webhook notifier failed for attacker=%s", attacker_id)

        # Return warning message to caller (so tests and integrations can surface it)
        return msg

    def _escalate_to_authorities(self, attacker_id: str) -> None:
        """Handle escalation after final warning.

        This function logs the escalation and marks the attacker in the mirror layer
        as 'escalated' so downstream systems can take action.
        """
        with self._lock:
            self.mirror_layer[attacker_id] = "ESCALATED_TO_AUTHORITIES"
        logger.error("Attacker %s escalated to authorities.", attacker_id)
        # Placeholder for further integration: syslog, webhook, SIEM, etc.
        # Also send a final webhook if configured
        try:
            if self.notifier_webhook_url:
                self._send_webhook_notification(attacker_id, "ESCALATED_TO_AUTHORITIES", 99)
        except Exception:
            logger.exception("Webhook final escalation failed for attacker=%s", attacker_id)

    def set_notification_recipients(self, space_leaf: Optional[str] = None, un: Optional[str] = None, president_of_the_united_states: Optional[str] = None) -> None:
        """Configure official notification webhooks for alerts.

        Provide webhook URLs (or None) for each recipient. These are used by the
        sponge service and escalation flow.
        """
        with self._lock:
            if space_leaf is not None:
                self.notification_recipients["space_leaf"] = space_leaf
            if un is not None:
                self.notification_recipients["un"] = un
            if president_of_the_united_states is not None:
                self.notification_recipients["president_of_the_united_states"] = president_of_the_united_states
        logger.info("Notification recipients updated: %s", self.notification_recipients)

    # Core admin methods
    def is_core_admin(self, user: str) -> bool:
        """Return True if `user` is a core admin (immutable core)."""
        return user in self.core_admins

    def update_core_config(self, user: str, key: str, value: object, consenting_guardians: Optional[List[str]] = None) -> bool:
        """Update a core configuration item.

        Rules:
         - If `user` is in `core_admins`, change requires explicit consenting_guardians
           that include at least one guardian (not just the core admin) unless
           `immutable_core` is False.
         - If `immutable_core` is False, core admins may change without extra consent.

        Returns True on success, False otherwise.
        """
        consenting_guardians = consenting_guardians or []

        # Only known keys allowed
        allowed_keys = {"lockdown_keys", "notification_recipients", "required_guardian_approvals"}
        if key not in allowed_keys:
            logger.warning("Attempt to update unknown core key=%s by user=%s", key, user)
            return False

        if not self.is_core_admin(user):
            logger.warning("User %s is not a core admin and cannot update core config", user)
            return False

        if not self.immutable_core:
            # apply change directly
            with self._lock:
                self._apply_core_change(key, value)
            logger.info("Core config key=%s updated by core admin=%s (immutable_core=False)", key, user)
            return True

        # immutable core: require consent from at least one guardian (excluding the admin if named)
        valid_consent = 0
        for g in set(consenting_guardians):
            if g in self.guardians and g != user:
                valid_consent += 1

        if valid_consent >= 1:
            with self._lock:
                self._apply_core_change(key, value)
            logger.info("Core config key=%s updated by core admin=%s with guardian consent=%s", key, user, consenting_guardians)
            return True

        logger.warning("Core config update denied for key=%s by user=%s: insufficient guardian consent", key, user)
        return False

    def _apply_core_change(self, key: str, value: object) -> None:
        """Internal helper to apply allowed core changes."""
        if key == "lockdown_keys":
            if isinstance(value, dict):
                from typing import cast
                value_dict: Dict[str, Any] = cast(Dict[str, Any], value)
                primary = value_dict.get("primary")
                secondary = value_dict.get("secondary_un")
                if primary is not None and secondary is not None:
                    self.set_lockdown_keys(str(primary), str(secondary))
            else:
                logger.warning("lockdown_keys value must be a dict")
        elif key == "notification_recipients":
            if isinstance(value, dict):
                value_dict: Dict[str, Any] = dict(value.items())
                def to_optional_str(val: Any) -> Optional[str]:
                    if val is None:
                        return None
                    return str(val)
                self.set_notification_recipients(
                    space_leaf=to_optional_str(value_dict.get("space_leaf")),
                    un=to_optional_str(value_dict.get("un")),
                    president_of_the_united_states=to_optional_str(value_dict.get("president_of_the_united_states")),
                )
            else:
                logger.warning("notification_recipients value must be a dict")
        elif key == "required_guardian_approvals":
            try:
                n = int(str(value))
                self.required_guardian_approvals = max(1, n)
            except Exception:
                logger.exception("Invalid value for required_guardian_approvals: %s", value)

    from typing import Dict, List, Any
    def run_sponge_cycle(self) -> Dict[str, Any]:
        """Run a single sponge/analysis cycle synchronously.

        This function safely inspects mirror_layer metadata (no payload exfiltration),
        applies simple heuristics (repeated intrusions → blocklist), updates firewall
        policies (e.g., increase strictness), and sends notifications via configured
        recipients and webhook/callbacks. Returns a summary dict.
        """
        summary: Dict[str, Any] = {"processed": 0, "escalated": [], "blocked": []}  # type: ignore
        with self._lock:
            attackers = list(self.mirror_layer.keys())

        for attacker in attackers:
            # Skip already escalated markers
            with self._lock:
                marker = self.mirror_layer.get(attacker)
            if marker == "ESCALATED_TO_AUTHORITIES":
                continue

            # Heuristic: if warnings >= 3 or repeated entries > 5 → escalate/block
            warnings = self.intruder_warnings.get(attacker, 0)
            # count repeated presence in mirror (naive check)
            repeated = 1

            # If severe, add to blocklist
            if warnings >= 3 or repeated > 5:
                with self._lock:
                    self.blocklist.add(attacker)
                    self.mirror_layer[attacker] = "ESCALATED_TO_AUTHORITIES"
                summary["blocked"].append(attacker)
                summary["escalated"].append(attacker)
                # notify recipients
                msg = f"Attacker {attacker} escalated and blocked"
                # send notifications via callback and webhook
                try:
                    if callable(self.notifier_callback):
                        self.notifier_callback(attacker, msg, 999)
                except Exception:
                    logger.exception("Notifier callback failed during sponge for %s", attacker)
                # send to configured official recipients
                for who, url in self.notification_recipients.items():
                    if url:
                        # POST a simple alert
                        try:
                            self.notifier_webhook_url, old = url, self.notifier_webhook_url
                            self.set_notifier_webhook(url, self.notifier_webhook_headers)
                            self._send_webhook_notification(attacker, msg, 999)
                            # restore previous notifier_webhook_url
                            self.set_notifier_webhook(old or "", self.notifier_webhook_headers)
                        except Exception:
                            logger.exception("Failed to notify %s about attacker %s", who, attacker)

            # Otherwise, perform lightweight signature analysis (hash checks) — simulated
            else:
                # Increase policy strictness as defensive measure
                with self._lock:
                    self.network_policies["internet"]["require_mfa"] = True
                # Notify about defensive update
                try:
                    if callable(self.notifier_callback):
                        self.notifier_callback(attacker, "Defensive policy hardened", 1)
                except Exception:
                    logger.exception("Notifier callback failed during soft-harden for %s", attacker)
            summary["processed"] += 1

        return summary

    def start_sponge_service(self, interval_seconds: float = 10.0) -> None:
        """Start background sponge thread that runs `run_sponge_cycle` every `interval_seconds`.

        This service is defensive only and inspects metadata from `mirror_layer`.
        """
        import threading

        def _loop():
            logger.info("Sponge service started (interval=%s)", interval_seconds)
            while self._sponge_running:
                try:
                    self.run_sponge_cycle()
                except Exception:
                    logger.exception("Sponge cycle failed")
                time.sleep(interval_seconds)

        with self._lock:
            if self._sponge_running:
                logger.info("Sponge already running")
                return
            self._sponge_running = True
            self._sponge_thread = threading.Thread(target=_loop, daemon=True)
            self._sponge_thread.start()

    def stop_sponge_service(self) -> None:
        """Stop the sponge background service."""
        with self._lock:
            self._sponge_running = False
        if self._sponge_thread:
            self._sponge_thread.join(timeout=2.0)
            self._sponge_thread = None
        logger.info("Sponge service stopped")

    def set_notifier_webhook(self, url: str, headers: Optional[Dict[str, str]] = None) -> None:
        """Configure a webhook URL to receive intrusion/warning notifications."""
        with self._lock:
            self.notifier_webhook_url = url
            self.notifier_webhook_headers = headers or {}
        logger.info("Notifier webhook configured: %s", url)

    def _send_webhook_notification(self, attacker_id: str, message: str, step: int) -> None:
        """Send a webhook notification asynchronously. Uses urllib to avoid external deps."""
        import json
        import threading
        import urllib.request

        # Ensure webhook URL is a non-empty string
        webhook_url = self.notifier_webhook_url if self.notifier_webhook_url else None
        if not webhook_url:
            logger.warning("Webhook URL is not set or not a string; skipping webhook notification.")
            return

        def _worker():
            try:
                payload = json.dumps({
                    "attacker_id": attacker_id,
                    "message": message,
                    "step": step,
                    "timestamp": time.time(),
                }).encode("utf-8")
                req = urllib.request.Request(webhook_url, data=payload, method="POST")
                for k, v in self.notifier_webhook_headers.items():
                    req.add_header(k, v)
                req.add_header("Content-Type", "application/json")
                # Fire and forget; read response to ensure delivery attempt
                with urllib.request.urlopen(req, timeout=5) as resp:
                    logger.debug("Webhook delivered status=%s", resp.status)
            except Exception:
                logger.exception("Failed to send webhook notification for attacker=%s", attacker_id)

        # Start background thread
        try:
            t = threading.Thread(target=_worker, daemon=True)
            t.start()
        except Exception:
            logger.exception("Failed to spawn webhook thread for attacker=%s", attacker_id)

    def collapse_cycle(self, purge_mirror: bool = True) -> str:
        """
        Firewall collapse + inversion cycle (self-renewal).
        Optionally purges mirror_layer (simulated cleanup).
        """
        if not self.active:
            logger.info("Collapse cycle requested while inactive — re-activating.")
            self.active = True

        timestamp = time.time()
        renewal_hash = self._make_hash(str(timestamp))
        with self._lock:
            self.diamond_layer["renewal"] = renewal_hash
            if purge_mirror:
                self.mirror_layer.clear()
                logger.info("Mirror layer purged during collapse cycle.")
            # Reset intruder warnings after a collapse cycle to avoid persistent escalation
            self.intruder_warnings.clear()
        logger.info("Diamond firewall collapsed inward — shield renewed.")
        return "Diamond firewall collapsed inward, anomalies purged, shield renewed."

    def access_request(
        self,
        dna_code: Optional[str] = None,
        network: str = "internet",
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        """
        Attempt to access true firewall. Supports network-aware checks.
        - dna_code: optional DNA QR string
        - network: one of 'wifi','hifi','ethernet','cellular','satellite','internet'
        - metadata: dict with network metadata (ssid, carrier, satellite_id, etc.)

        Raises AccessDeniedError on refusal. Returns success message on grant.
        """
        metadata = metadata or {}

        # If lockdown active, check admin keys first (they bypass network policy)
        if self.lockdown_active:
            admin_key_raw = metadata.get("admin_key") if metadata else None
            admin_key: str = str(admin_key_raw) if admin_key_raw is not None else ""
            if admin_key == self.lockdown_keys.get("primary"):
                logger.info("Lockdown bypass granted to primary key for network=%s", network)
                return f"Access granted by primary lockdown key for network={network}."
            if admin_key == self.lockdown_keys.get("secondary_un"):
                # UN gets limited access when lockdown is active
                if self.lockdown_un_limited:
                    logger.info("Lockdown limited access granted to UN secondary key for network=%s", network)
                    return f"Access granted by UN secondary key (limited access) for network={network}."
                else:
                    logger.info("Lockdown bypass granted to UN secondary key for network=%s", network)
                    return f"Access granted by UN secondary key for network={network}."
            logger.warning("Access denied due to active lockdown for network=%s", network)
            raise AccessDeniedError("Access denied due to active lockdown")

        # Evaluate network policy; if it fails, deny
        if not self._check_network_policy(network, {**metadata, **({"dna_code": dna_code} if dna_code else {})}):
            logger.warning("Access denied by network policy for network=%s", network)
            raise AccessDeniedError("Access denied by network policy")

        # If policy allows proceeding, perform dna check if provided
        if dna_code and self.dna_qr_check(dna_code):
            logger.info("Access granted: collective consent confirmed for network=%s", network)
            return f"Access granted: collective consent confirmed for network={network}."

        # If dna not provided but network is trusted by policy, allow
        policy = self.network_policies.get(network, {})
        if not dna_code and not policy.get("require_qr", False):
            logger.info("Access granted: network %s does not require QR", network)
            return f"Access granted: network {network} clear (no QR required)."

        logger.warning("Access denied: lineage-safe consent required or network policy blocked.")
        raise AccessDeniedError("Access denied: lineage-safe consent required or network policy blocked.")

    def is_trapped(self, attacker_id: str) -> Optional[str]:
        with self._lock:
            return self.mirror_layer.get(attacker_id)

    def __repr__(self) -> str:
        return (
            f"<DiamondFirewall guardians={len(self.guardians)} "
            f"captains={len(self.captains)} active={self.active}>"
        )

    # Lockdown control methods
    def set_lockdown_keys(self, primary: str, secondary_un: str) -> None:
        """Set the lockdown keys. `primary` is Space Leaf Corp main key; `secondary_un` is UN key."""
        with self._lock:
            self.lockdown_keys["primary"] = primary
            self.lockdown_keys["secondary_un"] = secondary_un
        logger.info("Lockdown keys set (primary and UN secondary)")

    def activate_lockdown_with_primary(self, key: str) -> bool:
        """Activate lockdown if primary key matches."""
        if key == self.lockdown_keys.get("primary"):
            with self._lock:
                self.lockdown_active = True
            logger.warning("Lockdown activated by primary key")
            return True
        logger.warning("Failed lockdown activation attempt with invalid primary key")
        return False

    def activate_lockdown_with_un(self, key: str) -> bool:
        """Activate lockdown if secondary UN key matches (UN can also activate)."""
        if key == self.lockdown_keys.get("secondary_un"):
            with self._lock:
                self.lockdown_active = True
                # Ensure UN activation results in limited mode (no full bypass)
                self.lockdown_un_limited = True
            logger.warning("Lockdown activated by UN secondary key (limited mode)")
            return True
        logger.warning("Failed lockdown activation attempt with invalid UN key")
        return False

    def deactivate_lockdown(self, key: str) -> bool:
        """Deactivate lockdown only with primary key (Space Leaf) or UN key if allowed."""
        if key in (self.lockdown_keys.get("primary"), self.lockdown_keys.get("secondary_un")):
            with self._lock:
                self.lockdown_active = False
                # reset limited flag when deactivated
                self.lockdown_un_limited = True
            logger.info("Lockdown deactivated by key")
            return True
        logger.warning("Failed attempt to deactivate lockdown with invalid key")
        return False


# Basic unit tests
if __name__ == "__main__":
    import unittest

    class TestDiamondFirewall(unittest.TestCase):
        def setUp(self):
            self.guardians = ["Miko", "JD", "Guardian_A", "Guardian_B"]
            self.captains = ["Captain_1", "Captain_2", "Captain_3", "Captain_4"]
            self.firewall = DiamondFirewall(self.guardians, self.captains, True)

        def test_dna_qr_check_pass(self):
            self.assertTrue(self.firewall.dna_qr_check("LINEAGE_SAFE_QR123"))

        def test_access_request_granted(self):
            self.assertEqual(
                self.firewall.access_request("LINEAGE_SAFE_QR123"),
                "Access granted: collective consent confirmed.",
            )

        def test_access_request_denied_raises(self):
            with self.assertRaises(AccessDeniedError):
                self.firewall.access_request("INVALID_QR")

        def test_intrusion_and_trap(self):
            res = self.firewall.intrusion_attempt("Hacker_007", "malware_payload")
            self.assertIn("trapped in mirror layer", res)
            self.assertIsNotNone(self.firewall.is_trapped("Hacker_007"))

        def test_collapse_cycle_purges_mirror(self):
            self.firewall.intrusion_attempt("Hacker_007", "payload")
            self.assertIn("Hacker_007", self.firewall.mirror_layer)
            self.firewall.collapse_cycle()
            self.assertNotIn("Hacker_007", self.firewall.mirror_layer)

    unittest.main()
