"""Tests for TimeEngine: exponential learning and spaced repetition."""
import time
import unittest
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from time_engine import (
    LearningItem,
    EbbinghausDecay,
    SM2Algorithm,
    ExponentialLearningCurve,
    TimeEngine,
)


class TestEbbinghausDecay(unittest.TestCase):
    """Test Ebbinghaus forgetting curve calculations."""

    def test_retention_at_zero_time(self):
        """Retention should equal strength at time 0."""
        strength = 0.8
        retention = EbbinghausDecay.retention_at_time(0, strength)
        self.assertAlmostEqual(retention, strength, places=5)

    def test_retention_decreases_over_time(self):
        """Retention should decrease as time increases."""
        strength = 0.8
        retention_1h = EbbinghausDecay.retention_at_time(1, strength)
        retention_24h = EbbinghausDecay.retention_at_time(24, strength)
        self.assertGreater(retention_1h, retention_24h)

    def test_retention_with_strong_item(self):
        """Strong items should retain better over time."""
        retention_weak_24h = EbbinghausDecay.retention_at_time(24, 0.3)
        retention_strong_24h = EbbinghausDecay.retention_at_time(24, 0.9)
        self.assertLess(retention_weak_24h, retention_strong_24h)

    def test_next_review_time_calculation(self):
        """Next review time should increase with strength."""
        strength = 0.8
        review_time = EbbinghausDecay.next_review_time(strength)
        self.assertGreater(review_time, 0)
        self.assertIsInstance(review_time, float)

    def test_next_review_time_stronger_items(self):
        """Stronger items should have longer review times."""
        time_weak = EbbinghausDecay.next_review_time(0.3)
        time_strong = EbbinghausDecay.next_review_time(0.9)
        self.assertLess(time_weak, time_strong)


class TestSM2Algorithm(unittest.TestCase):
    """Test SuperMemo 2 spaced repetition algorithm."""

    def test_failed_review_restarts_interval(self):
        """Failed review (quality < 3) should restart interval."""
        interval, ease = SM2Algorithm.next_interval(7.0, 2.5, 2)
        self.assertEqual(interval, 1.0)

    def test_successful_review_increases_interval(self):
        """Successful review (quality >= 4) should increase interval."""
        interval, ease = SM2Algorithm.next_interval(1.0, 2.5, 4)
        self.assertGreater(interval, 1.0)

    def test_ease_factor_minimum(self):
        """Ease factor should never go below 1.3."""
        _, ease = SM2Algorithm.next_interval(7.0, 1.3, 2)
        self.assertGreaterEqual(ease, 1.3)

    def test_perfect_review_increases_ease(self):
        """Perfect review (quality 5) should increase ease factor."""
        _, ease_good = SM2Algorithm.next_interval(7.0, 2.5, 4)
        _, ease_perfect = SM2Algorithm.next_interval(7.0, 2.5, 5)
        self.assertLess(ease_good, ease_perfect)

    def test_three_review_sequence(self):
        """First three reviews should follow SM-2 sequence: 1 → 3 → 7 days."""
        interval = 1.0
        ease = 2.5

        # Second review after 3 days
        interval, ease = SM2Algorithm.next_interval(interval, ease, 4)
        self.assertEqual(interval, 3.0)

        # Third review after 7 days
        interval, ease = SM2Algorithm.next_interval(interval, ease, 4)
        self.assertEqual(interval, 7.0)


class TestExponentialLearningCurve(unittest.TestCase):
    """Test exponential learning curve calculations."""

    def test_learning_progress_at_zero(self):
        """Progress should be 0 at time 0."""
        progress = ExponentialLearningCurve.learning_progress(0)
        self.assertAlmostEqual(progress, 0.0, places=5)

    def test_learning_progress_increases(self):
        """Progress should increase with time."""
        progress_1h = ExponentialLearningCurve.learning_progress(1)
        progress_10h = ExponentialLearningCurve.learning_progress(10)
        self.assertLess(progress_1h, progress_10h)

    def test_learning_progress_asymptote(self):
        """Progress should approach asymptote as time → ∞."""
        progress_large = ExponentialLearningCurve.learning_progress(1000)
        self.assertLess(progress_large, 1.0)
        self.assertGreater(progress_large, 0.99)

    def test_higher_rate_faster_learning(self):
        """Higher rate should lead to faster learning."""
        progress_slow = ExponentialLearningCurve.learning_progress(5, rate=0.05)
        progress_fast = ExponentialLearningCurve.learning_progress(5, rate=0.2)
        self.assertLess(progress_slow, progress_fast)

    def test_accelerated_growth(self):
        """Acceleration should modify progress."""
        base = 0.5
        accelerated = ExponentialLearningCurve.accelerated_growth(base)
        self.assertNotEqual(base, accelerated)
        self.assertIsInstance(accelerated, float)


class TestLearningItem(unittest.TestCase):
    """Test LearningItem dataclass."""

    def test_learning_item_creation(self):
        """Should create learning item with defaults."""
        item = LearningItem("id1", "content")
        self.assertEqual(item.item_id, "id1")
        self.assertEqual(item.content, "content")
        self.assertEqual(item.review_count, 0)
        self.assertEqual(item.retention_strength, 0.0)
        self.assertEqual(item.ease_factor, 2.5)

    def test_learning_item_timestamps(self):
        """Should track created_at and last_review."""
        item = LearningItem("id1", "content")
        self.assertGreater(item.created_at, 0)
        self.assertGreater(item.last_review, 0)
        self.assertAlmostEqual(item.created_at, item.last_review, delta=0.1)


class TestTimeEngine(unittest.TestCase):
    """Test TimeEngine main class."""

    def setUp(self):
        self.engine = TimeEngine()

    def test_add_item(self):
        """Should add learning items."""
        item = self.engine.add_item("id1", "Learn Python")
        self.assertIn("id1", self.engine.items)
        self.assertEqual(item.content, "Learn Python")

    def test_session_tracking(self):
        """Should track learning sessions."""
        self.engine.start_session()
        self.assertIsNotNone(self.engine.session_start)
        
        time.sleep(0.1)  # Sleep 100ms
        result = self.engine.end_session()
        
        self.assertGreater(result["session_duration_hours"], 0)
        self.assertGreater(result["total_learning_time_hours"], 0)

    def test_review_item_quality_affects_retention(self):
        """Quality ratings should affect retention strength."""
        self.engine.add_item("id1", "content")
        
        # Low quality review
        result_low = self.engine.review_item("id1", 2)
        low_retention = result_low["retention_strength"]
        
        # High quality review
        result_high = self.engine.review_item("id1", 5)
        high_retention = result_high["retention_strength"]
        
        self.assertLess(low_retention, high_retention)

    def test_review_count_increments(self):
        """Review count should increment."""
        self.engine.add_item("id1", "content")
        
        self.engine.review_item("id1", 4)
        self.assertEqual(self.engine.items["id1"].review_count, 1)
        
        self.engine.review_item("id1", 4)
        self.assertEqual(self.engine.items["id1"].review_count, 2)

    def test_get_next_reviews(self):
        """Should return items due for review."""
        self.engine.add_item("id1", "content1")
        self.engine.add_item("id2", "content2")
        
        # Set past last review for id1
        self.engine.items["id1"].last_review = time.time() - 3 * 24 * 3600
        self.engine.items["id1"].interval_days = 1.0
        
        next_reviews = self.engine.get_next_reviews()
        self.assertGreater(len(next_reviews), 0)
        self.assertEqual(next_reviews[0].item_id, "id1")

    def test_learning_analytics(self):
        """Should generate learning analytics."""
        self.engine.add_item("id1", "Learn A")
        self.engine.add_item("id2", "Learn B")
        
        self.engine.review_item("id1", 4)
        self.engine.review_item("id2", 5)
        
        analytics = self.engine.learning_analytics()
        self.assertEqual(analytics["total_items"], 2)
        self.assertEqual(analytics["total_reviews"], 2)
        self.assertGreater(analytics["average_retention"], 0)

    def test_export_state(self):
        """Should export complete learning state."""
        self.engine.add_item("id1", "Learn Python")
        self.engine.review_item("id1", 4)
        
        state = self.engine.export_state()
        
        self.assertIn("items", state)
        self.assertIn("id1", state["items"])
        self.assertIn("total_learning_time", state)
        self.assertIn("review_history", state)
        self.assertIn("analytics", state)

    def test_multiple_items_workflow(self):
        """Test complete workflow with multiple items."""
        # Add items
        self.engine.add_item("python", "Learn Python")
        self.engine.add_item("math", "Learn Calculus")
        self.engine.add_item("history", "Learn WW2")
        
        # Start session
        self.engine.start_session()
        
        # Review items
        for item_id in ["python", "math", "history"]:
            self.engine.review_item(item_id, 4)
        
        # End session
        result = self.engine.end_session()
        
        # Check analytics
        analytics = self.engine.learning_analytics()
        self.assertEqual(analytics["total_items"], 3)
        self.assertEqual(analytics["total_reviews"], 3)
        self.assertGreater(analytics["exponential_progress"], 0)

    def test_review_nonexistent_item(self):
        """Reviewing nonexistent item should return error."""
        result = self.engine.review_item("nonexistent", 4)
        self.assertIn("error", result)

    def test_end_session_without_start(self):
        """Ending session without start should return error."""
        result = self.engine.end_session()
        self.assertIn("error", result)

    def test_learning_acceleration_workflow(self):
        """Test learning with acceleration."""
        self.engine.add_item("accelerated", "Learning with boost")
        
        # Simulate multiple quality reviews
        for _ in range(5):
            self.engine.review_item("accelerated", 5)
        
        analytics = self.engine.learning_analytics()
        self.assertEqual(analytics["total_reviews"], 5)
        
        item = self.engine.items["accelerated"]
        self.assertGreater(item.ease_factor, 2.5)  # Ease increased
        self.assertGreater(item.retention_strength, 0)


if __name__ == "__main__":
    unittest.main()
