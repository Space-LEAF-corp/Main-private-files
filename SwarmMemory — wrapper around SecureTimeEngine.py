class SwarmMemory:
    """
    Swarm Memory v1.0
    - Stores events as TimeEngine items
    - Computes priority from severity * event weight
    - Exposes retraining batches for ML / control
    - Tracks per-node trust with decay
    """

    def __init__(self, engine: SecureTimeEngine, config: Optional[SwarmMemoryConfig] = None):
        self.engine = engine
        self.config = config or SwarmMemoryConfig()
        self.trust_scores: Dict[str, float] = {}   # node_id -> trust [0..1]
        self._items_count = 0

    # -----------------------------
    # Internal helpers
    # -----------------------------

    def _event_weight(self, event_type: str) -> float:
        return EVENT_WEIGHTS.get(event_type, EVENT_WEIGHTS["OTHER"])

    def _compute_priority(self, event: SwarmEvent) -> float:
        return event.severity * self._event_weight(event.event_type)

    def _ensure_capacity(self):
        # simple cap: if too many items, export/rotate externally
        if self._items_count > self.config.max_items:
            # in v1.0 we just stop adding; later: implement pruning
            return False
        return True

    # -----------------------------
    # Public API
    # -----------------------------

    def record_event(self, event: SwarmEvent):
        """Ingest a swarm event into TimeEngine if it matters enough."""
        priority = self._compute_priority(event)
        if priority < self.config.min_priority:
            return  # ignore trivial noise

        if not self._ensure_capacity():
            return

        item = {
            "timestamp": event.timestamp,
            "node_id": event.node_id,
            "event_type": event.event_type,
            "severity": event.severity,
            "priority": priority,
            "context": event.context,
        }

        self.engine.add_item(item)
        self._items_count += 1

        # Update trust immediately for trust-related events
        if event.event_type == "TRUST_ANOMALY":
            self._apply_trust_penalty(event.node_id, event.severity)

    def get_retraining_batch(self) -> List[Dict]:
        """
        Use TimeEngine's review schedule as 'what the swarm should think about now'.
        Returns a list of items to feed into ML retraining / control tuning.
        """
        reviews = self.engine.get_next_reviews()
        return reviews.get("next", [])

    # -----------------------------
    # Trust model
    # -----------------------------

    def _apply_trust_penalty(self, node_id: str, severity: float):
        """Lower trust based on anomaly severity."""
        current = self.trust_scores.get(node_id, 1.0)
        penalty = severity * 0.5  # up to -0.5 per severe anomaly
        self.trust_scores[node_id] = max(0.0, current - penalty)

    def _apply_trust_reward(self, node_id: str, quality: float):
        """Raise trust for good behavior / clean intervals."""
        current = self.trust_scores.get(node_id, 0.5)
        reward = quality * 0.2
        self.trust_scores[node_id] = min(1.0, current + reward)

    def decay_trust(self, now: Optional[float] = None):
        """
        Periodic trust decay: called e.g. once per day.
        Nodes drift toward 0.5 unless reinforced.
        """
        if now is None:
            now = time.time()

        # simple exponential decay toward 0.5
        half_life = self.config.trust_decay_half_life
        if half_life <= 0:
            return

        # decay factor per call; in v1.0 we assume this is called ~once per half-life
        decay_factor = 0.5

        for node_id, score in list(self.trust_scores.items()):
            # move score toward 0.5 by decay_factor
            delta = (0.5 - score) * decay_factor
            self.trust_scores[node_id] = max(0.0, min(1.0, score + delta))

    # -----------------------------
    # Export / introspection
    # -----------------------------

    def export_state(self) -> Dict:
        base = self.engine.export_state()
        base["trust_scores"] = self.trust_scores
        base["items_count"] = self._items_count
        return base
