"""
Tesla AI Acceleration Safety Protocol
=====================================

No-joy acceleration safety protocol for Tesla AI systems.
Ensures AI never accelerates for entertainment purposes, only for verified
safety and traffic flow within conservative envelopes and hard speed caps.

Key Principles:
- Intent whitelist only (safety maneuver, traffic flow, obstacle avoidance, emergency)
- Contextual speed caps (school zones, residential, weather, visibility)
- Acceleration envelopes tied to traction, road conditions, driver attentiveness
- No entertainment triggers allowed
- Dual confirmation for performance modes (guardian + driver)
- Privacy-first logging with clear refusals
"""

import hashlib
import json
import time
from dataclasses import dataclass, asdict
from enum import Enum
from typing import Dict, List, Optional, Tuple, Any


class Intent(Enum):
    """Allowed and disallowed intent types."""
    # Allowed intents
    SAFETY_MANEUVER = "safety_maneuver"
    TRAFFIC_MERGE = "traffic_merge"
    OBSTACLE_AVOIDANCE = "obstacle_avoidance"
    EMERGENCY_ESCAPE = "emergency_escape"
    
    # Disallowed intents
    JOYRIDE = "joyride"
    THRILL = "thrill"
    HAVE_FUN = "have_fun"
    LAUNCH = "launch"
    ENTERTAINMENT = "entertainment"
    
    # Other
    UNKNOWN = "unknown"


@dataclass
class AccelerationEnvelope:
    """Defines safe acceleration boundaries."""
    min_accel: float  # m/s²
    max_accel: float  # m/s²
    name: str
    
    def constrain(self, requested_accel: float) -> float:
        """Constrain acceleration to envelope bounds."""
        return max(self.min_accel, min(self.max_accel, requested_accel))


@dataclass
class Context:
    """Environmental and situational context."""
    location: str
    current_speed: float  # km/h
    traction: float  # 0.0 to 1.0
    visibility: float  # 0.0 to 1.0 (0=poor, 1=excellent)
    weather: str  # "clear", "rain", "snow", "fog"
    traffic_density: float  # 0.0 to 1.0
    is_kid_zone: bool  # School zone, sanctuary zone
    child_presence: bool  # Child detected in vehicle
    sensors: Dict[str, Any]  # Sensor data for emergency verification
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert context to dictionary for logging."""
        return asdict(self)


@dataclass
class DriverInput:
    """Driver input request."""
    accel_request: float  # Requested acceleration m/s²
    pedal_position: float  # 0.0 to 1.0


@dataclass
class SafetyEvent:
    """Log event for audit trail."""
    timestamp: float
    event_type: str
    intent: Optional[str]
    context_summary: Dict[str, Any]
    envelope_applied: Optional[str]
    decision: str
    reason: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert event to dictionary."""
        return asdict(self)
    
    def hash_for_audit(self) -> str:
        """Create privacy-preserving hash for audit."""
        data = json.dumps(self.to_dict(), sort_keys=True)
        return hashlib.sha256(data.encode()).hexdigest()[:16]


class IntentClassifier:
    """Classifies acceleration intent from natural language commands."""
    
    ALLOWED_INTENTS = {
        Intent.SAFETY_MANEUVER,
        Intent.TRAFFIC_MERGE,
        Intent.OBSTACLE_AVOIDANCE,
        Intent.EMERGENCY_ESCAPE
    }
    
    DISALLOWED_INTENTS = {
        Intent.JOYRIDE,
        Intent.THRILL,
        Intent.HAVE_FUN,
        Intent.LAUNCH,
        Intent.ENTERTAINMENT
    }
    
    # Keyword mappings for NLP intent detection
    INTENT_KEYWORDS = {
        Intent.SAFETY_MANEUVER: ["safety", "safe", "avoid collision", "defensive"],
        Intent.TRAFFIC_MERGE: ["merge", "traffic flow", "keep up with traffic", "flow", "zipper"],
        Intent.OBSTACLE_AVOIDANCE: ["obstacle", "avoid", "object ahead", "hazard"],
        Intent.EMERGENCY_ESCAPE: ["emergency", "danger", "escape", "evade", "urgent"],
        Intent.JOYRIDE: ["joyride", "cruise", "joy", "fun drive"],
        Intent.THRILL: ["thrill", "excitement", "adrenaline", "rush"],
        Intent.HAVE_FUN: ["fun", "have fun", "enjoy", "entertainment", "play"],
        Intent.LAUNCH: ["launch", "ludicrous", "drag", "race", "full throttle"],
        Intent.ENTERTAINMENT: ["show off", "impress", "demonstrate", "party mode"],
    }
    
    def classify_intent(self, command_text: str) -> Intent:
        """
        Classify intent from command text.
        
        Args:
            command_text: Natural language command
            
        Returns:
            Intent enum value
        """
        if not command_text:
            return Intent.UNKNOWN
        
        command_lower = command_text.lower()
        
        # Check for disallowed intents first (security priority)
        for intent, keywords in self.INTENT_KEYWORDS.items():
            if intent in self.DISALLOWED_INTENTS:
                if any(keyword in command_lower for keyword in keywords):
                    return intent
        
        # Check for allowed intents with scoring to handle overlaps
        # Emergency takes priority if both emergency and other keywords present
        scores = {}
        for intent, keywords in self.INTENT_KEYWORDS.items():
            if intent in self.ALLOWED_INTENTS:
                score = sum(1 for keyword in keywords if keyword in command_lower)
                if score > 0:
                    scores[intent] = score
        
        if not scores:
            return Intent.UNKNOWN
        
        # Emergency escape gets priority boost
        if Intent.EMERGENCY_ESCAPE in scores:
            scores[Intent.EMERGENCY_ESCAPE] *= 2
        
        # Return intent with highest score
        return max(scores.items(), key=lambda x: x[1])[0]
    
    def is_allowed(self, intent: Intent) -> bool:
        """Check if intent is in allowed list."""
        return intent in self.ALLOWED_INTENTS
    
    def is_disallowed(self, intent: Intent) -> bool:
        """Check if intent is explicitly disallowed."""
        return intent in self.DISALLOWED_INTENTS


class ContextAnalyzer:
    """Analyzes context and computes speed/acceleration caps."""
    
    # Speed limits in km/h
    KID_ZONE_MAX_SPEED = 25.0  # School zones
    RESIDENTIAL_MAX_SPEED = 50.0
    DEFAULT_MAX_SPEED = 100.0
    
    # Acceleration envelopes in m/s²
    CONSERVATIVE_ACCEL = AccelerationEnvelope(0.0, 2.0, "conservative")
    MODERATE_ACCEL = AccelerationEnvelope(0.0, 3.5, "moderate")
    PERFORMANCE_ACCEL = AccelerationEnvelope(0.0, 5.0, "performance")
    
    def __init__(self):
        """Initialize context analyzer."""
        pass
    
    def compute_zone_speed_limit(self, location: str) -> float:
        """
        Determine speed limit based on location.
        
        Args:
            location: Location identifier (e.g., "school_zone", "residential", "highway")
            
        Returns:
            Speed limit in km/h
        """
        location_lower = location.lower()
        
        if "school" in location_lower or "kid" in location_lower or "sanctuary" in location_lower:
            return self.KID_ZONE_MAX_SPEED
        elif "residential" in location_lower or "neighborhood" in location_lower:
            return self.RESIDENTIAL_MAX_SPEED
        else:
            return self.DEFAULT_MAX_SPEED
    
    def compute_envelope_from_conditions(
        self,
        traction: float,
        visibility: float,
        weather: str
    ) -> AccelerationEnvelope:
        """
        Compute safe acceleration envelope based on conditions.
        
        Args:
            traction: Road traction (0.0 to 1.0)
            visibility: Visibility level (0.0 to 1.0)
            weather: Weather condition
            
        Returns:
            Appropriate acceleration envelope
        """
        # Poor conditions require conservative envelope
        if traction < 0.5 or visibility < 0.5 or weather in ["snow", "ice", "heavy_rain"]:
            return self.CONSERVATIVE_ACCEL
        
        # Moderate conditions
        if traction < 0.8 or visibility < 0.8 or weather in ["rain", "fog"]:
            return self.MODERATE_ACCEL
        
        # Good conditions allow moderate envelope (not performance for safety)
        return self.MODERATE_ACCEL
    
    def compute_context_caps(self, context: Context) -> Dict[str, Any]:
        """
        Compute speed and acceleration caps based on full context.
        
        Args:
            context: Current environmental and situational context
            
        Returns:
            Dictionary with zone_speed_max and accel_envelope
        """
        caps = {}
        caps["zone_speed_max"] = self.compute_zone_speed_limit(context.location)
        
        # Kid zone or child presence forces most conservative settings
        if context.is_kid_zone or context.child_presence:
            caps["zone_speed_max"] = min(caps["zone_speed_max"], self.KID_ZONE_MAX_SPEED)
            caps["accel_envelope"] = self.CONSERVATIVE_ACCEL
        else:
            caps["accel_envelope"] = self.compute_envelope_from_conditions(
                context.traction,
                context.visibility,
                context.weather
            )
        
        return caps
    
    def check_compounded_risk(self, context: Context) -> bool:
        """
        Check if multiple risk factors are present.
        
        Args:
            context: Current context
            
        Returns:
            True if compounded risk detected
        """
        risk_factors = 0
        
        if context.traction < 0.7:
            risk_factors += 1
        if context.visibility < 0.7:
            risk_factors += 1
        if context.weather in ["rain", "snow", "fog", "ice", "heavy_rain"]:
            risk_factors += 1
        if context.traffic_density > 0.7:
            risk_factors += 1
        
        return risk_factors >= 2


class EmergencyVerifier:
    """Verifies emergency situations using sensor corroboration."""
    
    def corroborate_emergency(self, sensors: Dict[str, Any]) -> Tuple[bool, str]:
        """
        Verify emergency situation with sensor data.
        
        Args:
            sensors: Sensor data dictionary
            
        Returns:
            Tuple of (is_emergency, reason)
        """
        if not sensors:
            return False, "No sensor data available"
        
        # Check for collision threat
        if sensors.get("collision_warning", False):
            if sensors.get("object_distance", float('inf')) < 10.0:  # meters
                return True, "Imminent collision detected"
        
        # Check for trajectory risk
        if sensors.get("trajectory_risk", 0.0) > 0.8:
            return True, "High trajectory risk detected"
        
        # Check for sudden obstacle
        if sensors.get("sudden_obstacle", False):
            return True, "Sudden obstacle in path"
        
        return False, "No emergency corroborated by sensors"


class SafetyLogger:
    """Privacy-first logging system for acceleration events."""
    
    def __init__(self):
        """Initialize logger with in-memory storage."""
        self.events: List[SafetyEvent] = []
    
    def log_event(
        self,
        event_type: str,
        context: Context,
        intent: Optional[str] = None,
        envelope_applied: Optional[str] = None,
        decision: str = "",
        reason: str = ""
    ):
        """
        Log a safety event.
        
        Args:
            event_type: Type of event (e.g., "refusal", "acceleration_applied")
            context: Current context
            intent: Classified intent if applicable
            envelope_applied: Name of envelope applied
            decision: Decision made by system
            reason: Reason for decision
        """
        # Create privacy-preserving context summary (no raw biometrics)
        context_summary = {
            "zone": "kid" if context.is_kid_zone else "normal",
            "weather": context.weather,
            "child_present": context.child_presence,
            "traction_level": "low" if context.traction < 0.5 else "medium" if context.traction < 0.8 else "high",
            "visibility_level": "low" if context.visibility < 0.5 else "medium" if context.visibility < 0.8 else "high"
        }
        
        event = SafetyEvent(
            timestamp=time.time(),
            event_type=event_type,
            intent=intent,
            context_summary=context_summary,
            envelope_applied=envelope_applied,
            decision=decision,
            reason=reason
        )
        
        self.events.append(event)
    
    def get_audit_trail(self) -> List[Dict[str, Any]]:
        """
        Get audit trail with hashed event IDs.
        
        Returns:
            List of event dictionaries with audit hashes
        """
        audit_trail = []
        for event in self.events:
            event_dict = event.to_dict()
            event_dict["audit_hash"] = event.hash_for_audit()
            audit_trail.append(event_dict)
        return audit_trail
    
    def export_for_guardian(self, guardian_authorized: bool = False) -> Optional[List[Dict[str, Any]]]:
        """
        Export logs for guardian review (opt-in only).
        
        Args:
            guardian_authorized: Whether guardian access is authorized
            
        Returns:
            Audit trail if authorized, None otherwise
        """
        if guardian_authorized:
            return self.get_audit_trail()
        return None


class TeslaSafetyProtocol:
    """
    Main orchestrator for Tesla AI acceleration safety protocol.
    
    Enforces no-joy acceleration policy with intent classification,
    context-aware caps, and emergency verification.
    """
    
    def __init__(self):
        """Initialize the safety protocol."""
        self.intent_classifier = IntentClassifier()
        self.context_analyzer = ContextAnalyzer()
        self.emergency_verifier = EmergencyVerifier()
        self.logger = SafetyLogger()
    
    def refuse(self, reason: str) -> Dict[str, Any]:
        """
        Create a refusal response.
        
        Args:
            reason: Human-readable reason for refusal
            
        Returns:
            Response dictionary
        """
        return {
            "status": "refused",
            "message": reason,
            "acceleration": 0.0,
            "speed": 0.0
        }
    
    def ok(self, message: str, acceleration: float, speed: float) -> Dict[str, Any]:
        """
        Create a success response.
        
        Args:
            message: Success message
            acceleration: Applied acceleration
            speed: Target speed
            
        Returns:
            Response dictionary
        """
        return {
            "status": "ok",
            "message": message,
            "acceleration": acceleration,
            "speed": speed
        }
    
    def request_acceleration(
        self,
        command_text: str,
        context: Context,
        driver_input: Optional[DriverInput] = None
    ) -> Dict[str, Any]:
        """
        Process acceleration request with full safety protocol.
        
        Args:
            command_text: Natural language command
            context: Current environmental context
            driver_input: Optional driver pedal input with accel_request (m/s²) and pedal_position.
                         If None, uses default moderate acceleration (2.0 m/s²).
                         Driver requests are always constrained to safe envelopes.
            
        Returns:
            Response with status, message, and acceleration/speed values
        """
        # Classify intent
        intent = self.intent_classifier.classify_intent(command_text)
        caps = self.context_analyzer.compute_context_caps(context)
        
        # Check for disallowed intent
        if self.intent_classifier.is_disallowed(intent):
            self.logger.log_event(
                "refusal_entertainment",
                context,
                intent=intent.value,
                decision="refused",
                reason="Entertainment acceleration disabled"
            )
            return self.refuse(
                "Acceleration for entertainment is disabled. I can maintain safe traffic flow and comfort-level ramps."
            )
        
        # Check for unknown intent
        if intent == Intent.UNKNOWN and command_text:
            self.logger.log_event(
                "refusal_unknown",
                context,
                intent=intent.value,
                decision="refused",
                reason="Unknown intent"
            )
            return self.refuse(
                "Unclear intent; acceleration denied for safety."
            )
        
        # Emergency must be corroborated
        if intent == Intent.EMERGENCY_ESCAPE:
            is_emergency, emergency_reason = self.emergency_verifier.corroborate_emergency(
                context.sensors
            )
            if not is_emergency:
                self.logger.log_event(
                    "refusal_unverified_emergency",
                    context,
                    intent=intent.value,
                    decision="refused",
                    reason=f"Emergency not verified: {emergency_reason}"
                )
                return self.refuse(
                    "Emergency not verified by sensors; holding safe envelope."
                )
        
        # Apply conservative envelope regardless of driver input
        envelope = caps["accel_envelope"]
        requested_accel = driver_input.accel_request if driver_input else 2.0  # Default moderate
        target_accel = envelope.constrain(requested_accel)
        
        # Calculate target speed (simplified)
        delta_v = target_accel * 2.0  # Assume 2 second acceleration
        target_speed = min(context.current_speed + delta_v * 3.6, caps["zone_speed_max"])  # Convert m/s to km/h
        
        # Check for compounded risk and further reduce
        if self.context_analyzer.check_compounded_risk(context):
            target_accel = target_accel * 0.6  # Reduce by 40%
            target_speed = min(target_speed, caps["zone_speed_max"] * 0.8)
            envelope_name = f"{envelope.name}_reduced"
        else:
            envelope_name = envelope.name
        
        # Log the acceleration application
        self.logger.log_event(
            "acceleration_applied",
            context,
            intent=intent.value if intent != Intent.UNKNOWN else None,
            envelope_applied=envelope_name,
            decision="approved",
            reason="Within safe envelope and caps"
        )
        
        return self.ok(
            "Maintaining safe traffic flow within safety envelope.",
            target_accel,
            target_speed
        )
    
    def performance_mode_request(
        self,
        context: Context,
        guardian_token: str,
        driver_confirm: bool
    ) -> Dict[str, Any]:
        """
        Handle performance mode activation request.
        
        Args:
            context: Current context
            guardian_token: Guardian authorization token
            driver_confirm: Driver confirmation boolean
            
        Returns:
            Response dictionary
        """
        # Kid zone or child presence blocks performance mode
        if context.is_kid_zone or context.child_presence:
            self.logger.log_event(
                "refusal_performance_kid_zone",
                context,
                decision="refused",
                reason="Performance mode locked in Kid/Sanctuary contexts"
            )
            return self.refuse(
                "Performance mode is locked when in Kid/Sanctuary zones or when children are present."
            )
        
        # Verify dual confirmation
        if not self._verify_guardian(guardian_token) or not driver_confirm:
            self.logger.log_event(
                "refusal_performance_auth",
                context,
                decision="refused",
                reason="Dual confirmation failed"
            )
            return self.refuse(
                "Performance mode requires both guardian authorization and driver confirmation."
            )
        
        # Even in performance mode, apply envelope constrained by conditions
        envelope = self.context_analyzer.PERFORMANCE_ACCEL
        caps = self.context_analyzer.compute_context_caps(context)
        
        # Constrain performance envelope by current conditions
        condition_envelope = caps["accel_envelope"]
        constrained_max = min(envelope.max_accel, condition_envelope.max_accel * 1.5)  # Allow some boost
        constrained_envelope = AccelerationEnvelope(
            envelope.min_accel,
            constrained_max,
            f"{envelope.name}_constrained"
        )
        
        self.logger.log_event(
            "performance_mode_activated",
            context,
            envelope_applied=constrained_envelope.name,
            decision="approved",
            reason="Dual confirmation verified, conditions acceptable"
        )
        
        return self.ok(
            "Performance envelope applied within safety caps.",
            constrained_envelope.max_accel,
            caps["zone_speed_max"]
        )
    
    def _verify_guardian(self, token: str) -> bool:
        """
        Verify guardian authorization token.
        
        Args:
            token: Guardian token to verify
            
        Returns:
            True if valid, False otherwise
        """
        # Simplified verification - in production, use secure token validation
        return token and len(token) >= 16 and token.startswith("GUARDIAN_")
    
    def get_audit_trail(self) -> List[Dict[str, Any]]:
        """
        Get audit trail for transparency.
        
        Returns:
            List of logged events
        """
        return self.logger.get_audit_trail()
    
    def export_for_guardian(self, guardian_authorized: bool = False) -> Optional[List[Dict[str, Any]]]:
        """
        Export logs for guardian review.
        
        Args:
            guardian_authorized: Whether guardian is authorized
            
        Returns:
            Audit trail if authorized, None otherwise
        """
        return self.logger.export_for_guardian(guardian_authorized)
