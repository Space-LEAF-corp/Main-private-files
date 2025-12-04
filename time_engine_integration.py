"""Integration between TimeEngine and Jarvondis for intelligent memory growth.

This module enables Jarvondis to learn exponentially by tracking which topics
grow in strength over time through spaced repetition.
"""
from __future__ import annotations

import time
from typing import Dict, Optional
from jarvondis.jarvondis import Jarvondis
from time_engine import TimeEngine


class JarvondisTimeIntegration:
    """Bridges Jarvondis memory system with TimeEngine learning curves."""

    def __init__(self, jarvondis: Jarvondis, time_engine: TimeEngine):
        self.jarvondis = jarvondis
        self.engine = time_engine
        self.topic_index: Dict[str, str] = {}  # topic -> item_id mapping

    def track_interaction(
        self,
        topic: str,
        input_text: str,
        output_text: str,
        quality_rating: int = 4,
    ) -> Dict:
        """Track a Jarvondis interaction and register with TimeEngine.
        
        Args:
            topic: Learning topic (e.g., "Python", "Mathematics")
            input_text: User input/question
            output_text: Jarvondis response
            quality_rating: How well the response was (0-5)
        
        Returns:
            Integration result with learning metrics
        """
        # Save to Jarvondis memory
        self.jarvondis.remember(input_text, output_text, topic)
        
        # Track in TimeEngine
        item_id = f"{topic}_{int(time.time())}"
        
        if topic not in self.topic_index:
            # First time seeing this topic
            self.engine.add_item(item_id, f"Topic: {topic}")
            self.topic_index[topic] = item_id
        else:
            item_id = self.topic_index[topic]
        
        # Review the item with quality rating
        review_result = self.engine.review_item(item_id, quality_rating)
        
        return {
            "jarvondis_memory": {
                "input": input_text,
                "output": output_text,
                "topic": topic,
            },
            "time_engine": review_result,
            "learning_state": self.engine.learning_analytics(),
        }

    def get_learning_priority_topics(self, limit: int = 5) -> list:
        """Get topics most needing review, sorted by priority.
        
        Returns:
            List of (topic, learning_item) tuples
        """
        due_items = self.engine.get_next_reviews(limit * 2)
        
        result = []
        for item in due_items[:limit]:
            # Find topic from item_id
            for topic, iid in self.topic_index.items():
                if iid == item.item_id:
                    result.append((topic, item))
                    break
        
        return result

    def exponential_learning_report(self) -> Dict:
        """Generate exponential learning analysis."""
        analytics = self.engine.learning_analytics()
        
        # Get memory stats
        memory_count = len(self.jarvondis.memory)
        
        return {
            "memory_interactions": memory_count,
            "total_learning_hours": analytics["total_learning_hours"],
            "exponential_progress": analytics["exponential_progress"],
            "topics_tracked": len(self.topic_index),
            "items_requiring_review": analytics["items_due_for_review"],
            "average_retention": analytics["average_retention"],
            "average_ease_factor": analytics["average_ease"],
            "engine_state": self.engine.export_state(),
        }

    def suggest_next_learning_session(self) -> Dict:
        """Suggest what to focus on in next learning session."""
        due_topics = self.get_learning_priority_topics(limit=3)
        
        if not due_topics:
            return {
                "suggestion": "Add new topics to learn",
                "action": "Call track_interaction() with new topics",
            }
        
        return {
            "priority_topics": [(t, i.content) for t, i in due_topics],
            "total_due_for_review": len(self.engine.get_next_reviews()),
            "estimated_session_minutes": len(due_topics) * 5,
            "learning_progress": self.engine.learning_analytics()["exponential_progress"],
        }

    def export_integrated_state(self) -> Dict:
        """Export complete integrated learning state."""
        return {
            "jarvondis_memory": {
                "total_interactions": len(self.jarvondis.memory),
                "interactions": self.jarvondis.memory,
            },
            "time_engine": self.engine.export_state(),
            "topic_mapping": self.topic_index,
            "learning_report": self.exponential_learning_report(),
        }
