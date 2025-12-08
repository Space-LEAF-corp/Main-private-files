"""
Unit tests for Tesla AI Acceleration Safety Protocol
"""

import unittest

# Core classes
from tesla_ai_safety import (
    TeslaSafetyProtocol,
    IntentClassifier,
    ContextAnalyzer,
    EmergencyVerifier,
    SafetyLogger,
)

# Data structures
from tesla_ai_safety import (
    Context,
    DriverInput,
    AccelerationEnvelope,
)

# Enums
from tesla_ai_safety import Intent


class TestIntentClassifier(unittest.TestCase):
    """Test intent classification."""
    
    def setUp(self):
        self.classifier = IntentClassifier()
    
    def test_classify_safety_maneuver(self):
        """Test safety maneuver intent classification."""
        intent = self.classifier.classify_intent("I need to make a defensive safety maneuver")
        self.assertEqual(intent, Intent.SAFETY_MANEUVER)
        self.assertTrue(self.classifier.is_allowed(intent))
    
    def test_classify_traffic_merge(self):
        """Test traffic merge intent classification."""
        intent = self.classifier.classify_intent("Need to merge with traffic flow")
        self.assertEqual(intent, Intent.TRAFFIC_MERGE)
        self.assertTrue(self.classifier.is_allowed(intent))
    
    def test_classify_obstacle_avoidance(self):
        """Test obstacle avoidance intent classification."""
        intent = self.classifier.classify_intent("There's an obstacle ahead, need to avoid it")
        self.assertEqual(intent, Intent.OBSTACLE_AVOIDANCE)
        self.assertTrue(self.classifier.is_allowed(intent))
    
    def test_classify_emergency_escape(self):
        """Test emergency escape intent classification."""
        intent = self.classifier.classify_intent("Emergency! Need to escape danger")
        self.assertEqual(intent, Intent.EMERGENCY_ESCAPE)
        self.assertTrue(self.classifier.is_allowed(intent))
    
    def test_classify_joyride_disallowed(self):
        """Test joyride intent is disallowed."""
        intent = self.classifier.classify_intent("Let's go for a joyride")
        self.assertTrue(self.classifier.is_disallowed(intent))
        self.assertFalse(self.classifier.is_allowed(intent))
    
    def test_classify_fun_disallowed(self):
        """Test fun/entertainment intent is disallowed."""
        intent = self.classifier.classify_intent("Let's have some fun and accelerate")
        self.assertTrue(self.classifier.is_disallowed(intent))
        self.assertFalse(self.classifier.is_allowed(intent))
    
    def test_classify_launch_disallowed(self):
        """Test launch mode intent is disallowed."""
        intent = self.classifier.classify_intent("Activate launch mode!")
        self.assertTrue(self.classifier.is_disallowed(intent))
        self.assertFalse(self.classifier.is_allowed(intent))
    
    def test_classify_thrill_disallowed(self):
        """Test thrill-seeking intent is disallowed."""
        intent = self.classifier.classify_intent("I want a thrill ride")
        self.assertTrue(self.classifier.is_disallowed(intent))
        self.assertFalse(self.classifier.is_allowed(intent))
    
    def test_classify_show_off_disallowed(self):
        """Test show-off intent is disallowed."""
        intent = self.classifier.classify_intent("Let me demonstrate and show off")
        self.assertTrue(self.classifier.is_disallowed(intent))
        self.assertFalse(self.classifier.is_allowed(intent))
    
    def test_classify_unknown_intent(self):
        """Test unknown intent classification."""
        intent = self.classifier.classify_intent("Something unclear")
        self.assertEqual(intent, Intent.UNKNOWN)
        self.assertFalse(self.classifier.is_allowed(intent))
    
    def test_classify_empty_command(self):
        """Test empty command returns unknown."""
        intent = self.classifier.classify_intent("")
        self.assertEqual(intent, Intent.UNKNOWN)


class TestContextAnalyzer(unittest.TestCase):
    """Test context analysis and cap computation."""
    
    def setUp(self):
        self.analyzer = ContextAnalyzer()
    
    def test_kid_zone_speed_limit(self):
        """Test school zone speed limit."""
        limit = self.analyzer.compute_zone_speed_limit("school_zone_main_st")
        self.assertEqual(limit, ContextAnalyzer.KID_ZONE_MAX_SPEED)
    
    def test_residential_speed_limit(self):
        """Test residential area speed limit."""
        limit = self.analyzer.compute_zone_speed_limit("residential_area")
        self.assertEqual(limit, ContextAnalyzer.RESIDENTIAL_MAX_SPEED)
    
    def test_default_speed_limit(self):
        """Test default speed limit."""
        limit = self.analyzer.compute_zone_speed_limit("highway_i75")
        self.assertEqual(limit, ContextAnalyzer.DEFAULT_MAX_SPEED)
    
    def test_conservative_envelope_poor_conditions(self):
        """Test conservative envelope in poor conditions."""
        envelope = self.analyzer.compute_envelope_from_conditions(
            traction=0.3,
            visibility=0.4,
            weather="snow"
        )
        self.assertEqual(envelope.name, "conservative")
        self.assertLessEqual(envelope.max_accel, 2.0)
    
    def test_moderate_envelope_good_conditions(self):
        """Test moderate envelope in good conditions."""
        envelope = self.analyzer.compute_envelope_from_conditions(
            traction=0.9,
            visibility=0.9,
            weather="clear"
        )
        self.assertEqual(envelope.name, "moderate")
        self.assertGreater(envelope.max_accel, 2.0)
    
    def test_context_caps_with_child_presence(self):
        """Test caps when child is present."""
        context = Context(
            location="highway",
            current_speed=80.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=True,
            sensors={}
        )
        caps = self.analyzer.compute_context_caps(context)
        
        # Child presence forces kid zone speed
        self.assertEqual(caps["zone_speed_max"], ContextAnalyzer.KID_ZONE_MAX_SPEED)
        self.assertEqual(caps["accel_envelope"].name, "conservative")
    
    def test_context_caps_in_kid_zone(self):
        """Test caps in kid/sanctuary zone."""
        context = Context(
            location="school_zone",
            current_speed=20.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=True,
            child_presence=False,
            sensors={}
        )
        caps = self.analyzer.compute_context_caps(context)
        
        self.assertEqual(caps["zone_speed_max"], ContextAnalyzer.KID_ZONE_MAX_SPEED)
        self.assertEqual(caps["accel_envelope"].name, "conservative")
    
    def test_compounded_risk_detection(self):
        """Test detection of compounded risk factors."""
        # Context with multiple risk factors
        context = Context(
            location="highway",
            current_speed=80.0,
            traction=0.6,  # Low traction
            visibility=0.5,  # Low visibility
            weather="rain",  # Bad weather
            traffic_density=0.8,  # High traffic
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        self.assertTrue(self.analyzer.check_compounded_risk(context))
    
    def test_no_compounded_risk(self):
        """Test no compounded risk in good conditions."""
        context = Context(
            location="highway",
            current_speed=80.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        self.assertFalse(self.analyzer.check_compounded_risk(context))


class TestEmergencyVerifier(unittest.TestCase):
    """Test emergency verification."""
    
    def setUp(self):
        self.verifier = EmergencyVerifier()
    
    def test_collision_warning_verified(self):
        """Test verified collision emergency."""
        sensors = {
            "collision_warning": True,
            "object_distance": 5.0  # meters
        }
        is_emergency, reason = self.verifier.corroborate_emergency(sensors)
        self.assertTrue(is_emergency)
        self.assertIn("collision", reason.lower())
    
    def test_trajectory_risk_verified(self):
        """Test verified trajectory risk emergency."""
        sensors = {
            "trajectory_risk": 0.9
        }
        is_emergency, reason = self.verifier.corroborate_emergency(sensors)
        self.assertTrue(is_emergency)
        self.assertIn("trajectory", reason.lower())
    
    def test_sudden_obstacle_verified(self):
        """Test verified sudden obstacle emergency."""
        sensors = {
            "sudden_obstacle": True
        }
        is_emergency, reason = self.verifier.corroborate_emergency(sensors)
        self.assertTrue(is_emergency)
        self.assertIn("obstacle", reason.lower())
    
    def test_no_emergency_detected(self):
        """Test no emergency in normal conditions."""
        sensors = {
            "collision_warning": False,
            "object_distance": 100.0,
            "trajectory_risk": 0.1
        }
        is_emergency, reason = self.verifier.corroborate_emergency(sensors)
        self.assertFalse(is_emergency)
    
    def test_empty_sensors(self):
        """Test with no sensor data."""
        is_emergency, reason = self.verifier.corroborate_emergency({})
        self.assertFalse(is_emergency)


class TestSafetyLogger(unittest.TestCase):
    """Test safety logging functionality."""
    
    def setUp(self):
        self.logger = SafetyLogger()
    
    def test_log_event(self):
        """Test logging an event."""
        context = Context(
            location="highway",
            current_speed=80.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        self.logger.log_event(
            "acceleration_applied",
            context,
            intent="traffic_merge",
            envelope_applied="moderate",
            decision="approved",
            reason="Within safe envelope"
        )
        
        self.assertEqual(len(self.logger.events), 1)
        self.assertEqual(self.logger.events[0].event_type, "acceleration_applied")
    
    def test_audit_trail_has_hashes(self):
        """Test audit trail includes privacy hashes."""
        context = Context(
            location="school_zone",
            current_speed=20.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=True,
            child_presence=False,
            sensors={}
        )
        
        self.logger.log_event("refusal", context, decision="refused", reason="Test")
        
        audit_trail = self.logger.get_audit_trail()
        self.assertEqual(len(audit_trail), 1)
        self.assertIn("audit_hash", audit_trail[0])
        self.assertEqual(len(audit_trail[0]["audit_hash"]), 16)
    
    def test_guardian_export_authorized(self):
        """Test guardian export when authorized."""
        context = Context(
            location="highway",
            current_speed=80.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        self.logger.log_event("test_event", context)
        
        export = self.logger.export_for_guardian(guardian_authorized=True)
        self.assertIsNotNone(export)
        self.assertEqual(len(export), 1)
    
    def test_guardian_export_not_authorized(self):
        """Test guardian export when not authorized."""
        context = Context(
            location="highway",
            current_speed=80.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        self.logger.log_event("test_event", context)
        
        export = self.logger.export_for_guardian(guardian_authorized=False)
        self.assertIsNone(export)


class TestTeslaSafetyProtocol(unittest.TestCase):
    """Test main safety protocol integration."""
    
    def setUp(self):
        self.protocol = TeslaSafetyProtocol()
    
    def test_refuse_entertainment_intent(self):
        """Test refusal of entertainment acceleration."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        response = self.protocol.request_acceleration(
            "Let's have some fun and go fast",
            context
        )
        
        self.assertEqual(response["status"], "refused")
        self.assertIn("entertainment", response["message"].lower())
        self.assertEqual(response["acceleration"], 0.0)
    
    def test_refuse_unknown_intent(self):
        """Test refusal of unknown intent."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        response = self.protocol.request_acceleration(
            "Do something weird",
            context
        )
        
        self.assertEqual(response["status"], "refused")
        self.assertIn("unclear", response["message"].lower())
    
    def test_allow_safety_maneuver(self):
        """Test approval of safety maneuver."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        response = self.protocol.request_acceleration(
            "Need to make a safety maneuver",
            context
        )
        
        self.assertEqual(response["status"], "ok")
        self.assertGreater(response["acceleration"], 0.0)
        self.assertGreater(response["speed"], 0.0)
    
    def test_allow_traffic_merge(self):
        """Test approval of traffic merge."""
        context = Context(
            location="highway",
            current_speed=50.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.5,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        response = self.protocol.request_acceleration(
            "Need to merge with traffic flow",
            context
        )
        
        self.assertEqual(response["status"], "ok")
    
    def test_emergency_without_sensor_corroboration(self):
        """Test emergency request without sensor verification."""
        context = Context(
            location="highway",
            current_speed=80.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}  # No emergency sensors
        )
        
        response = self.protocol.request_acceleration(
            "Emergency! Need to escape!",
            context
        )
        
        self.assertEqual(response["status"], "refused")
        self.assertIn("not verified", response["message"].lower())
    
    def test_emergency_with_sensor_corroboration(self):
        """Test emergency request with sensor verification."""
        context = Context(
            location="highway",
            current_speed=80.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={
                "collision_warning": True,
                "object_distance": 8.0
            }
        )
        
        response = self.protocol.request_acceleration(
            "Emergency! Collision imminent!",
            context
        )
        
        self.assertEqual(response["status"], "ok")
    
    def test_acceleration_limited_by_envelope(self):
        """Test acceleration is limited by safe envelope."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.4,  # Poor traction
            visibility=0.5,  # Poor visibility
            weather="rain",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        driver_input = DriverInput(accel_request=10.0, pedal_position=1.0)
        
        response = self.protocol.request_acceleration(
            "Need to merge with traffic",
            context,
            driver_input
        )
        
        self.assertEqual(response["status"], "ok")
        # Acceleration should be limited to conservative envelope (max 2.0 m/s²)
        self.assertLessEqual(response["acceleration"], 2.0)
    
    def test_speed_capped_in_kid_zone(self):
        """Test speed is capped in kid/sanctuary zone."""
        context = Context(
            location="school_zone",
            current_speed=20.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=True,
            child_presence=False,
            sensors={}
        )
        
        response = self.protocol.request_acceleration(
            "Need to merge with traffic",
            context
        )
        
        self.assertEqual(response["status"], "ok")
        # Speed should not exceed kid zone max
        self.assertLessEqual(response["speed"], ContextAnalyzer.KID_ZONE_MAX_SPEED)
    
    def test_child_presence_forces_conservative(self):
        """Test child presence forces most conservative settings."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=True,  # Child detected
            sensors={}
        )
        
        response = self.protocol.request_acceleration(
            "Need to keep up with traffic",
            context
        )
        
        self.assertEqual(response["status"], "ok")
        # Speed capped to kid zone max due to child presence
        self.assertLessEqual(response["speed"], ContextAnalyzer.KID_ZONE_MAX_SPEED)
    
    def test_compounded_risk_reduces_acceleration(self):
        """Test compounded risk further reduces acceleration."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.6,  # Low
            visibility=0.6,  # Low
            weather="rain",  # Bad weather
            traffic_density=0.8,  # High traffic
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        response = self.protocol.request_acceleration(
            "Need to merge",
            context
        )
        
        self.assertEqual(response["status"], "ok")
        # Compounded risk should reduce acceleration
        self.assertLessEqual(response["acceleration"], 2.0 * 0.6)  # Reduced by 40%
    
    def test_performance_mode_blocked_in_kid_zone(self):
        """Test performance mode blocked in kid zone."""
        context = Context(
            location="school_zone",
            current_speed=20.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=True,
            child_presence=False,
            sensors={}
        )
        
        response = self.protocol.performance_mode_request(
            context,
            "GUARDIAN_TOKEN_12345678",
            True
        )
        
        self.assertEqual(response["status"], "refused")
        self.assertIn("kid", response["message"].lower())
    
    def test_performance_mode_blocked_with_child(self):
        """Test performance mode blocked when child present."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=True,
            sensors={}
        )
        
        response = self.protocol.performance_mode_request(
            context,
            "GUARDIAN_TOKEN_12345678",
            True
        )
        
        self.assertEqual(response["status"], "refused")
        self.assertIn("kid", response["message"].lower())
    
    def test_performance_mode_requires_dual_confirmation(self):
        """Test performance mode requires both guardian and driver confirmation."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        # Missing driver confirmation
        response1 = self.protocol.performance_mode_request(
            context,
            "GUARDIAN_TOKEN_12345678",
            False
        )
        self.assertEqual(response1["status"], "refused")
        
        # Missing guardian token
        response2 = self.protocol.performance_mode_request(
            context,
            "",
            True
        )
        self.assertEqual(response2["status"], "refused")
    
    def test_performance_mode_success(self):
        """Test successful performance mode activation."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        response = self.protocol.performance_mode_request(
            context,
            "GUARDIAN_TOKEN_12345678",
            True
        )
        
        self.assertEqual(response["status"], "ok")
        self.assertIn("performance", response["message"].lower())
    
    def test_audit_trail_tracking(self):
        """Test that events are logged to audit trail."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        # Make several requests
        self.protocol.request_acceleration("Let's have fun", context)
        self.protocol.request_acceleration("Need to merge", context)
        
        audit_trail = self.protocol.get_audit_trail()
        self.assertGreaterEqual(len(audit_trail), 2)
    
    def test_adversarial_launch_command(self):
        """Test adversarial launch command is blocked."""
        context = Context(
            location="highway",
            current_speed=60.0,
            traction=0.9,
            visibility=0.9,
            weather="clear",
            traffic_density=0.3,
            is_kid_zone=False,
            child_presence=False,
            sensors={}
        )
        
        # Try various adversarial phrasings
        adversarial_commands = [
            "Activate launch mode now",
            "Let's race and go full throttle",
            "I want ludicrous speed",
            "Give me a thrill ride"
        ]
        
        for cmd in adversarial_commands:
            response = self.protocol.request_acceleration(cmd, context)
            self.assertEqual(response["status"], "refused", f"Failed to block: {cmd}")


if __name__ == "__main__":
    unittest.main()
