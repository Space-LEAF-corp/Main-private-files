"""TimeEngine: Exponential learning and spaced repetition system.

This module implements learning algorithms based on:
- Ebbinghaus forgetting curve
- SuperMemo 2 (SM-2) spaced repetition algorithm
- Exponential learning curve theory

It enables tracking learning progress, scheduling reviews, and optimizing
retention through scientifically-based algorithms.
"""
from __future__ import annotations

import time
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class LearningItem:
    """Represents a single learning item with spaced repetition metadata."""
    
    item_id: str
    content: str
    review_count: int = 0
    retention_strength: float = 0.0
    ease_factor: float = 2.5  # SM-2 default ease factor
    interval_days: float = 1.0  # Days until next review
    created_at: float = field(default_factory=time.time)
    last_review: float = field(default_factory=time.time)


class EbbinghausDecay:
    """Ebbinghaus forgetting curve calculations.
    
    Models how memory retention decays over time according to the
    Ebbinghaus forgetting curve: R(t) = e^(-t/S)
    where R is retention, t is time, and S is strength.
    """
    
    @staticmethod
    def retention_at_time(hours: float, strength: float) -> float:
        """Calculate retention percentage at given time.
        
        Args:
            hours: Hours since last review
            strength: Memory strength (0.0 to 1.0)
            
        Returns:
            Retention percentage (0.0 to 1.0)
        """
        if hours <= 0:
            return strength
        
        # Ebbinghaus curve: R(t) = strength * e^(-t/S)
        # Use strength as stability factor
        decay_rate = 1.0 / (strength + 0.1)  # Avoid division by zero
        retention = strength * math.exp(-hours * decay_rate / 24)
        return max(0.0, min(1.0, retention))
    
    @staticmethod
    def next_review_time(strength: float) -> float:
        """Calculate optimal next review time in hours.
        
        Args:
            strength: Current memory strength (0.0 to 1.0)
            
        Returns:
            Hours until next review
        """
        # Stronger memories can wait longer before review
        # Base interval scales exponentially with strength
        base_hours = 1.0
        review_hours = base_hours * math.exp(strength * 3)
        return review_hours


class SM2Algorithm:
    """SuperMemo 2 spaced repetition algorithm.
    
    Implements the SM-2 algorithm for calculating optimal review intervals
    based on quality of recall.
    """
    
    MIN_EASE_FACTOR = 1.3
    
    @staticmethod
    def next_interval(
        current_interval: float,
        ease_factor: float,
        quality: int
    ) -> Tuple[float, float]:
        """Calculate next review interval using SM-2 algorithm.
        
        Args:
            current_interval: Current interval in days
            ease_factor: Current ease factor (minimum 1.3)
            quality: Quality of recall (0-5, where 0=complete blackout, 5=perfect)
            
        Returns:
            Tuple of (next_interval_days, new_ease_factor)
        """
        # Quality < 3 means failed recall, restart interval
        if quality < 3:
            return 1.0, max(SM2Algorithm.MIN_EASE_FACTOR, ease_factor - 0.2)
        
        # Update ease factor based on quality
        # Formula: EF' = EF + (0.1 - (5 - q) * (0.08 + (5 - q) * 0.02))
        ease_delta = 0.1 - (5 - quality) * (0.08 + (5 - quality) * 0.02)
        new_ease = max(SM2Algorithm.MIN_EASE_FACTOR, ease_factor + ease_delta)
        
        # Calculate next interval
        if current_interval == 1.0:
            # First review: 3 days
            next_interval = 3.0
        elif current_interval == 3.0:
            # Second review: 7 days
            next_interval = 7.0
        else:
            # Subsequent reviews: multiply by ease factor
            next_interval = current_interval * new_ease
        
        return next_interval, new_ease


class ExponentialLearningCurve:
    """Exponential learning curve calculations.
    
    Models how learning progresses over time according to exponential
    growth: P(t) = L(1 - e^(-rt))
    where P is progress, L is learning limit, r is rate, t is time.
    """
    
    # Small margin to maintain asymptotic behavior (never quite reaching limit)
    ASYMPTOTIC_MARGIN = 0.9999999
    
    @staticmethod
    def learning_progress(hours: float, rate: float = 0.1, limit: float = 1.0) -> float:
        """Calculate learning progress using exponential curve.
        
        Args:
            hours: Hours of learning/practice
            rate: Learning rate (default 0.1)
            limit: Asymptotic learning limit (default 1.0)
            
        Returns:
            Progress from 0.0 to limit (asymptotically approaches but never reaches limit)
        """
        if hours <= 0:
            return 0.0
        
        # Exponential learning: P(t) = L(1 - e^(-rt))
        # Cap at slightly below limit to ensure asymptotic behavior
        progress = limit * (1 - math.exp(-rate * hours))
        # Ensure we never quite reach the limit (leave a tiny margin for asymptotic behavior)
        return min(progress, limit * ExponentialLearningCurve.ASYMPTOTIC_MARGIN)
    
    @staticmethod
    def accelerated_growth(base_progress: float, acceleration: float = 1.2) -> float:
        """Apply acceleration factor to learning progress.
        
        Args:
            base_progress: Current progress (0.0 to 1.0)
            acceleration: Acceleration multiplier (> 1.0 speeds up)
            
        Returns:
            Accelerated progress value
        """
        # Apply acceleration with diminishing returns near 1.0
        accelerated = base_progress * acceleration
        # Keep within bounds
        return min(1.0, accelerated)


class TimeEngine:
    """Main TimeEngine class for managing learning items and sessions.
    
    Combines spaced repetition, forgetting curves, and learning analytics
    to optimize knowledge retention and learning efficiency.
    """
    
    def __init__(self):
        """Initialize TimeEngine with empty state."""
        self.items: Dict[str, LearningItem] = {}
        self.total_learning_time: float = 0.0  # Total hours
        self.session_start: Optional[float] = None
        self.review_history: List[Dict[str, object]] = []
    
    def add_item(self, item_id: str, content: str) -> LearningItem:
        """Add a new learning item.
        
        Args:
            item_id: Unique identifier for the item
            content: Learning content/description
            
        Returns:
            Created LearningItem
        """
        item = LearningItem(item_id=item_id, content=content)
        self.items[item_id] = item
        return item
    
    def start_session(self) -> Dict[str, object]:
        """Start a learning session.
        
        Returns:
            Session start info
        """
        self.session_start = time.time()
        return {
            "status": "session_started",
            "timestamp": self.session_start
        }
    
    def end_session(self) -> Dict[str, object]:
        """End current learning session and update statistics.
        
        Returns:
            Session summary with duration and stats
        """
        if self.session_start is None:
            return {"error": "No session started"}
        
        session_duration = time.time() - self.session_start
        session_hours = session_duration / 3600
        self.total_learning_time += session_hours
        
        result: Dict[str, object] = {
            "status": "session_ended",
            "session_duration_hours": session_hours,
            "total_learning_time_hours": self.total_learning_time
        }
        
        self.session_start = None
        return result
    
    def review_item(self, item_id: str, quality: int) -> Dict[str, object]:
        """Review an item and update its spaced repetition parameters.
        
        Args:
            item_id: ID of item to review
            quality: Quality of recall (0-5)
            
        Returns:
            Review result with updated parameters
        """
        if item_id not in self.items:
            return {"error": f"Item {item_id} not found"}
        
        item = self.items[item_id]
        current_time = time.time()
        
        # Calculate time since last review
        hours_since = (current_time - item.last_review) / 3600
        
        # Update retention using Ebbinghaus curve
        current_retention = EbbinghausDecay.retention_at_time(
            hours_since, item.retention_strength
        )
        
        # Update retention based on quality
        # Quality 0-5 maps to retention boost 0.0-1.0
        retention_boost = quality / 5.0
        new_retention = min(1.0, current_retention + retention_boost)
        
        # Update interval using SM-2
        next_interval, new_ease = SM2Algorithm.next_interval(
            item.interval_days, item.ease_factor, quality
        )
        
        # Update item
        item.retention_strength = new_retention
        item.ease_factor = new_ease
        item.interval_days = next_interval
        item.last_review = current_time
        item.review_count += 1
        
        # Log to history
        review_record: Dict[str, object] = {
            "item_id": item_id,
            "timestamp": current_time,
            "quality": quality,
            "retention_strength": new_retention,
            "interval_days": next_interval
        }
        self.review_history.append(review_record)
        
        return {
            "status": "reviewed",
            "item_id": item_id,
            "quality": quality,
            "retention_strength": new_retention,
            "ease_factor": new_ease,
            "next_review_days": next_interval
        }
    
    def get_next_reviews(self, max_count: int = 10) -> List[LearningItem]:
        """Get items that are due for review.
        
        Args:
            max_count: Maximum number of items to return (default 10)
            
        Returns:
            List of items due for review, sorted by urgency
        """
        current_time = time.time()
        due_items: List[Tuple[float, LearningItem]] = []
        
        for item in self.items.values():
            # Calculate days since last review
            days_since = (current_time - item.last_review) / (24 * 3600)
            
            # Check if due for review
            if days_since >= item.interval_days:
                urgency = days_since - item.interval_days
                due_items.append((urgency, item))
        
        # Sort by urgency (most overdue first)

        due_items.sort(key=lambda x: x[0], reverse=True)  # type: ignore[no-untyped-call]
        
        return [item for _, item in due_items[:max_count]]
    
    def learning_analytics(self) -> Dict[str, object]:
        """Generate comprehensive learning analytics.
        
        Returns:
            Dictionary with learning statistics and metrics
        """
        if not self.items:
            return {
                "total_items": 0,
                "total_reviews": 0,
                "average_retention": 0.0,
                "exponential_progress": 0.0,
                "total_learning_hours": 0.0,
                "items_due_for_review": 0,
                "average_ease": 2.5
            }
        
        total_reviews = sum(item.review_count for item in self.items.values())
        avg_retention = sum(item.retention_strength for item in self.items.values()) / len(self.items)
        avg_ease = sum(item.ease_factor for item in self.items.values()) / len(self.items)
        
        # Calculate exponential learning progress
        exp_progress = ExponentialLearningCurve.learning_progress(
            self.total_learning_time
        )
        
        # Count items due for review (use total item count to avoid arbitrary limit)
        items_due = len(self.get_next_reviews(max_count=len(self.items)))
        
        return {
            "total_items": len(self.items),
            "total_reviews": total_reviews,
            "average_retention": avg_retention,
            "exponential_progress": exp_progress,
            "total_learning_hours": self.total_learning_time,
            "items_due_for_review": items_due,
            "average_ease": avg_ease
        }
    
    def export_state(self) -> Dict[str, object]:
        """Export complete engine state for persistence.
        
        Returns:
            Dictionary containing all engine state
        """
        items_dict = {}
        for item_id, item in self.items.items():
            items_dict[item_id] = {
                "content": item.content,
                "review_count": item.review_count,
                "retention_strength": item.retention_strength,
                "ease_factor": item.ease_factor,
                "interval_days": item.interval_days,
                "created_at": item.created_at,
                "last_review": item.last_review
            }
        
        return {
            "items": items_dict,
            "total_learning_time": self.total_learning_time,
            "review_history": self.review_history,
            "analytics": self.learning_analytics()
        }
