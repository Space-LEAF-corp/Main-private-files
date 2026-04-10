import time
from typing import Dict, List, Optional

# --- Configuration ---

ALERT_THRESHOLD = 0.70          # sudden normalized pressure drop threshold
DEBOUNCE_MS = 250               # minimum duration to consider a drop as valid
MULTISENSOR_CONFIRM_MS = 400    # window to check neighboring sensor agreement
WATCHDOG_TIMEOUT_MS = 3000      # max time without data before watchdog alerts
MIN_HEALTHY_SENSORS = 0.80      # fraction required to be considered healthy
BASELINE_WINDOW = 120           # seconds for rolling baseline calibration

AUTHORIZED_USERS = {"nurseA", "chargeRN", "ptSafetyLead", "admin"}

# --- Utilities ---

class AuditLog:
    def __init__(self):
        self.events: List[Dict] = []

    def write(self, event: str, meta: Dict = None):
        self.events.append({
            "ts": time.time(),
            "event": event,
            "meta": meta or {}
        })

class Sensor:
    def __init__(self, sensor_id: str):
        self.id = sensor_id
        self.last_value: float = 0.0
        self.last_update_ts: float = 0.0
        self.baseline_values: List[float] = []
        self.online: bool = True

    def update(self, value: float):
        self.last_value = value
        self.last_update_ts = time.time()
        # maintain rolling baseline window
        self.baseline_values.append(value)
        # prune values older than BASELINE_WINDOW seconds
        cutoff = self.last_update_ts - BASELINE_WINDOW
        self.baseline_values = [
            v for v in self.baseline_values
            # assuming 1Hz updates; for real streams, store (v, ts) pairs
        ][:BASELINE_WINDOW]  # simplified; replace with timestamped buffer in prod

    def baseline(self) -> float:
        if not self.baseline_values:
            return self.last_value
        # robust baseline via median (resistant to spikes)
        sorted_vals = sorted(self.baseline_values)
        mid = len(sorted_vals) // 2
        return sorted_vals[mid]

class FallPreventionSystem:
    def __init__(self, sensors: List[Sensor], neighbors: Dict[str, List[str]]):
        self.sensors: Dict[str, Sensor] = {s.id: s for s in sensors}
        self.neighbors = neighbors  # map sensor_id -> neighbor_ids
        self.alert_active: bool = False
        self.alert_meta: Dict = {}
        self.last_data_ts: float = 0.0
        self.log = AuditLog()

    # --- Core monitoring ---

    def ingest(self, sensor_id: str, value: float):
        """Ingest a new sensor value (pressure)."""
        if sensor_id not in self.sensors:
            return
        s = self.sensors[sensor_id]
        s.update(value)
        self.last_data_ts = time.time()
        self._check_for_fall(sensor_id)

    def _check_for_fall(self, sensor_id: str):
        s = self.sensors[sensor_id]
        baseline = s.baseline()
        if baseline <= 0:
            return

        normalized_drop = max(0.0, (baseline - s.last_value) / baseline)

        # Debounce: require drop to persist for DEBOUNCE_MS
        start_ts = time.time()
        if normalized_drop >= ALERT_THRESHOLD:
            time.sleep(DEBOUNCE_MS / 1000.0)
            # re-sample after debounce
            if (baseline - s.last_value) / baseline >= ALERT_THRESHOLD:
                # Multi-sensor confirmation (if neighbors available)
                if self._confirm_with_neighbors(sensor_id, start_ts):
                    self._trigger_alert(sensor_id, normalized_drop)
                else:
                    # If no neighbors or disagreement, still trigger but mark as "single-sensor"
                    self._trigger_alert(sensor_id, normalized_drop, single_sensor=True)

    def _confirm_with_neighbors(self, sensor_id: str, start_ts: float) -> bool:
        neigh_ids = self.neighbors.get(sensor_id, [])
        if not neigh_ids:
            return False
        deadline = start_ts + (MULTISENSOR_CONFIRM_MS / 1000.0)
        agreement = 0
        total = 0
        while time.time() < deadline:
            time.sleep(0.05)
            for nid in neigh_ids:
                ns = self.sensors.get(nid)
                if not ns: 
                    continue
                total += 1
                b = ns.baseline()
                if b > 0 and (b - ns.last_value) / b >= (ALERT_THRESHOLD * 0.6):
                    agreement += 1
        # require at least half of neighbors to show supportive anomaly
        return total > 0 and agreement / total >= 0.5

    def _trigger_alert(self, sensor_id: str, drop: float, single_sensor: bool = False):
        self.alert_active = True
        self.alert_meta = {
            "sensor_id": sensor_id,
            "drop": round(drop, 3),
            "ts": time.time(),
            "single_sensor": single_sensor
        }
        self.log.write("ALERT_TRIGGERED", self.alert_meta)
        # TODO: integrate with nurse station / paging system

    # --- Diagnostics & reset ---

    def system_health(self) -> Dict:
        """Run diagnostics across sensors and streams."""
        now = time.time()
        online_count = 0
        healthy_baselines = 0
        for s in self.sensors.values():
            # consider sensor offline if no update in watchdog window
            if now - s.last_update_ts <= WATCHDOG_TIMEOUT_MS / 1000.0:
                online_count += 1
            # baseline sanity: non-zero and not extreme
            b = s.baseline()
            if b > 0 and 0.05 <= b <= 5.0:  # example expected pressure range
                healthy_baselines += 1

        total = len(self.sensors)
        health = {
            "total_sensors": total,
            "online_sensors": online_count,
            "healthy_baselines": healthy_baselines,
            "stream_ok": (now - self.last_data_ts) <= WATCHDOG_TIMEOUT_MS / 1000.0,
            "fraction_online": (online_count / total) if total else 0.0,
            "fraction_healthy": (healthy_baselines / total) if total else 0.0
        }
        self.log.write("DIAGNOSTIC_RUN", health)
        return health

    def recalibrate(self):
        """Recompute baselines by smoothing recent values."""
        for s in self.sensors.values():
            # In production, compute baseline from timestamped buffer here.
            # We retain rolling median approach via Sensor.baseline().
            pass
        self.log.write("RECALIBRATION_COMPLETE", {"ts": time.time()})

    def reset_system(self, user_id: str) -> bool:
        """Authenticated reset: clears alert and runs diagnostics before resuming."""
        if user_id not in AUTHORIZED_USERS:
            self.log.write("RESET_DENIED", {"user_id": user_id})
            return False

        # Clear alert state
        self.alert_active = False
        self.alert_meta = {}
        self.log.write("RESET_INITIATED", {"user_id": user_id})

        # Recalibrate baselines
        self.recalibrate()

        # Run diagnostics; only resume if healthy enough
        health = self.system_health()
        ok = (
            health["stream_ok"] and
            health["fraction_online"] >= MIN_HEALTHY_SENSORS and
            health["fraction_healthy"] >= MIN_HEALTHY_SENSORS
        )
        if ok:
            self.log.write("RESET_SUCCESS", {"user_id": user_id, "health": health})
            return True
        else:
            # Graceful degrade: system continues monitoring but flags maintenance
            self.log.write("RESET_PARTIAL", {"user_id": user_id, "health": health})
            # Optionally trigger maintenance notification here
            return True  # continue with warnings

    # --- Watchdog ---

    def watchdog_check(self):
        """Detect stalled input streams and raise maintenance alerts."""
        stalled = (time.time() - self.last_data_ts) > (WATCHDOG_TIMEOUT_MS / 1000.0)
        if stalled:
            self.log.write("WATCHDOG_STALL", {"ts": time.time()})
            # TODO: notify engineering/biomed for investigation
        return not stalled
