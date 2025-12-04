"""Tests for TimeEngine-Jarvondis integration."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest
from time_engine import TimeEngine
from time_engine_integration import JarvondisTimeIntegration
from jarvondis.jarvondis import Jarvondis


class TestJarvondisTimeIntegration(unittest.TestCase):
    """Test integration between Jarvondis and TimeEngine."""

    def setUp(self):
        self.jarvondis = Jarvondis()
        self.engine = TimeEngine()
        self.integration = JarvondisTimeIntegration(self.jarvondis, self.engine)

    def test_track_interaction_creates_memory_and_learning_item(self):
        """Tracking interaction should update both Jarvondis and TimeEngine."""
        result = self.integration.track_interaction(
            topic="Python",
            input_text="What is Python?",
            output_text="Python is a programming language",
            quality_rating=5,
        )
        
        self.assertIn("jarvondis_memory", result)
        self.assertIn("time_engine", result)
        self.assertIn("learning_state", result)

    def test_topic_index_tracking(self):
        """Topics should be indexed for easy lookup."""
        self.integration.track_interaction(
            topic="Physics",
            input_text="What is velocity?",
            output_text="Velocity is rate of change",
            quality_rating=4,
        )
        
        self.assertIn("Physics", self.integration.topic_index)

    def test_multiple_interactions_same_topic(self):
        """Multiple interactions on same topic should strengthen retention."""
        for i in range(3):
            self.integration.track_interaction(
                topic="Math",
                input_text=f"Math question {i}",
                output_text=f"Math answer {i}",
                quality_rating=5,
            )
        
        # Should only have one item for Math
        math_items = [t for t in self.integration.topic_index if "Math" in t]
        self.assertEqual(len(math_items), 1)

    def test_jarvondis_memory_persists(self):
        """Interactions should be saved in Jarvondis memory."""
        self.integration.track_interaction(
            topic="History",
            input_text="What happened in 1066?",
            output_text="The Norman Conquest",
            quality_rating=4,
        )
        
        # Check Jarvondis has the interaction
        self.assertGreater(len(self.jarvondis.memory), 0)

    def test_get_learning_priority_topics(self):
        """Should return topics in priority order."""
        # Add multiple interactions
        for topic in ["Python", "Math", "History"]:
            self.integration.track_interaction(
                topic=topic,
                input_text=f"Learn {topic}",
                output_text=f"Info about {topic}",
                quality_rating=4,
            )
        
        priority = self.integration.get_learning_priority_topics(limit=2)
        self.assertLessEqual(len(priority), 2)
        self.assertGreater(len(priority), 0)

    def test_exponential_learning_report(self):
        """Should generate comprehensive learning report."""
        self.integration.track_interaction(
            topic="Calculus",
            input_text="What is a derivative?",
            output_text="Rate of change",
            quality_rating=5,
        )
        
        report = self.integration.exponential_learning_report()
        
        self.assertIn("memory_interactions", report)
        self.assertIn("total_learning_hours", report)
        self.assertIn("exponential_progress", report)
        self.assertIn("topics_tracked", report)
        self.assertIn("average_retention", report)

    def test_suggest_next_learning_session(self):
        """Should suggest next learning session."""
        self.integration.track_interaction(
            topic="Science",
            input_text="What is gravity?",
            output_text="Fundamental force",
            quality_rating=3,  # Lower quality to mark for review
        )
        
        suggestion = self.integration.suggest_next_learning_session()
        self.assertIn("priority_topics", suggestion)
        self.assertIn("total_due_for_review", suggestion)

    def test_export_integrated_state(self):
        """Should export complete integrated state."""
        self.integration.track_interaction(
            topic="Art",
            input_text="What is impressionism?",
            output_text="Art movement",
            quality_rating=4,
        )
        
        state = self.integration.export_integrated_state()
        
        self.assertIn("jarvondis_memory", state)
        self.assertIn("time_engine", state)
        self.assertIn("topic_mapping", state)
        self.assertIn("learning_report", state)

    def test_learning_progression_over_multiple_reviews(self):
        """Track how retention improves with multiple quality reviews."""
        topic = "Programming"
        
        # First review
        result1 = self.integration.track_interaction(
            topic=topic,
            input_text="Review 1",
            output_text="Good response",
            quality_rating=4,
        )
        retention1 = result1["time_engine"]["retention_strength"]
        
        # Second review
        result2 = self.integration.track_interaction(
            topic=topic,
            input_text="Review 2",
            output_text="Better response",
            quality_rating=5,
        )
        retention2 = result2["time_engine"]["retention_strength"]
        
        # Retention should improve with quality reviews
        self.assertLessEqual(retention1, retention2)

    def test_quality_rating_affects_learning(self):
        """Different quality ratings should affect learning differently."""
        # Low quality interaction
        result_low = self.integration.track_interaction(
            topic="Topic_A",
            input_text="Question",
            output_text="Poor response",
            quality_rating=1,
        )
        
        # High quality interaction
        result_high = self.integration.track_interaction(
            topic="Topic_B",
            input_text="Question",
            output_text="Excellent response",
            quality_rating=5,
        )
        
        low_retention = result_low["time_engine"]["retention_strength"]
        high_retention = result_high["time_engine"]["retention_strength"]
        
        self.assertLess(low_retention, high_retention)

    def test_exponential_progress_increases_with_sessions(self):
        """Exponential progress should increase as learning time accumulates."""
        # First session
        self.engine.start_session()
        self.integration.track_interaction(
            topic="Session1",
            input_text="Q",
            output_text="A",
            quality_rating=5,
        )
        report1 = self.integration.exponential_learning_report()
        progress1 = report1["exponential_progress"]
        
        # Second session (after some time)
        self.engine.start_session()
        for i in range(3):
            self.integration.track_interaction(
                topic=f"Session2_{i}",
                input_text="Q",
                output_text="A",
                quality_rating=5,
            )
        self.engine.end_session()
        
        report2 = self.integration.exponential_learning_report()
        progress2 = report2["exponential_progress"]
        
        self.assertGreaterEqual(progress2, progress1)


if __name__ == "__main__":
    unittest.main()
