# uknptts_system_kids_v1.py
# Kid-facing conceptual model for the UK/UN Planetary Transportation Train System (UKNPTTS)
# Focus: story, observation, learning. No real control, just safe simulation snapshots.

from enum import Enum
import random
import math

# -----------------------------
# Friendly enums
# -----------------------------

class World(Enum):
    EARTH = "Earth"
    MOON = "Moon"
    MARS = "Mars"

class TrackPieceState(Enum):
    RESTING = 0
    WIGGLING_STRAIGHT = 1
    LOCKED_STRAIGHT = 2
    HAZARD_PAUSE = 3

class RailType(Enum):
    MAGNET_FLOAT = 0

# -----------------------------
# Friendly track pieces
# -----------------------------

class TrackPiece:
    """
    One little piece of magic track that can wiggle itself straight.
    Kids see it as a 'smart rail tile'.
    """
    def __init__(self, piece_id, position, target_position):
        self.piece_id = piece_id
        self.position = position
        self.target_position = target_position
        self.state = TrackPieceState.RESTING
        self.locked = False
        self.hazard_flag = False
        self.rail_type = RailType.MAGNET_FLOAT

    def alignment_error(self):
        dx = self.target_position[0] - self.position[0]
        dy = self.target_position[1] - self.position[1]
        dz = self.target_position[2] - self.position[2]
        error = math.sqrt(dx*dx + dy*dy + dz*dz)
        return error, (dx, dy, dz)

    def check_for_hazard(self):
        # Tiny chance of a pretend hazard so kids can see safety in action
        self.hazard_flag = (random.random() < 0.001)
        return self.hazard_flag

    def wiggle_toward_target(self, delta):
        step_scale = 0.1
        self.position = (
            self.position[0] + delta[0] * step_scale,
            self.position[1] + delta[1] * step_scale,
            self.position[2] + delta[2] * step_scale,
        )

    def lock_straight(self):
        self.locked = True
        self.state = TrackPieceState.LOCKED_STRAIGHT

    def unlock(self):
        self.locked = False
        self.state = TrackPieceState.RESTING

    def update(self, max_allowed_error):
        if self.check_for_hazard():
            self.state = TrackPieceState.HAZARD_PAUSE
            self.locked = True
            return

        if self.state == TrackPieceState.HAZARD_PAUSE:
            return

        error, delta = self.alignment_error()
        if error <= max_allowed_error:
            if not self.locked:
                self.lock_straight()
        else:
            self.state = TrackPieceState.WIGGLING_STRAIGHT
            self.unlock()
            self.wiggle_toward_target(delta)


class SmartTrack:
    """
    A whole line of smart track pieces that kids can watch as they wiggle into a straight line.
    """
    def __init__(self, max_allowed_error=0.002):
        self.pieces = []
        self.max_allowed_error = max_allowed_error
        self.global_pause = False

    def add_piece(self, piece):
        self.pieces.append(piece)

    def update_all(self):
        for p in self.pieces:
            p.update(self.max_allowed_error)

    def check_safety(self):
        self.global_pause = any(p.state == TrackPieceState.HAZARD_PAUSE for p in self.pieces)

    def step(self):
        self.update_all()
        self.check_safety()
        return {
            "global_pause": self.global_pause,
            "pieces": [
                {
                    "id": p.piece_id,
                    "state": p.state.name,
                    "locked": p.locked,
                    "hazard": p.hazard_flag,
                }
                for p in self.pieces
            ],
        }

# -----------------------------
# Friendly tunneling robot
# -----------------------------

class TunnelBuddy:
    """
    A friendly tunneling robot that digs safe tunnels and gently places smart track pieces.
    Kids only watch its progress.
    """
    def __init__(self, name="TunnelBuddy-1", world=World.MOON):
        self.name = name
        self.world = world
        self.progress_m = 0.0

    def dig_step(self, meters=1.0):
        self.progress_m += meters

    def place_track_piece(self, smart_track):
        new_id = len(smart_track.pieces)
        target_pos = (self.progress_m, 0.0, 0.0)
        start_pos = (
            target_pos[0] + random.uniform(-0.05, 0.05),
            target_pos[1] + random.uniform(-0.05, 0.05),
            target_pos[2] + random.uniform(-0.05, 0.05),
        )
        piece = TrackPiece(piece_id=new_id, position=start_pos, target_position=target_pos)
        smart_track.add_piece(piece)

    def step(self, smart_track):
        self.dig_step(meters=1.0)
        self.place_track_piece(smart_track)

# -----------------------------
# Friendly worker train
# -----------------------------

class CarKind(Enum):
    BREAK_ROOM = 0
    BUNK_WOMEN = 1
    BUNK_MEN = 2
    WORKER_RIDE = 3

class CozyCar:
    """
    A cozy car where workers can ride, rest, or relax.
    Kids learn that grown-ups need safe, calm places too.
    """
    def __init__(self, car_id, kind, capacity):
        self.car_id = car_id
        self.kind = kind
        self.capacity = capacity
        self.occupants = 0

    def status(self):
        return {
            "car_id": self.car_id,
            "kind": self.kind.name,
            "capacity": self.capacity,
            "occupants": self.occupants,
        }


class CozyTrain:
    """
    The worker train:
    - Women’s bunk car
    - Men’s bunk car
    - Shared break room
    - Worker ride car
    """
    def __init__(self):
        self.cars = []
        self._build_default()

    def _build_default(self):
        self.cars.append(CozyCar(0, CarKind.BUNK_WOMEN, capacity=20))
        self.cars.append(CozyCar(1, CarKind.BUNK_MEN, capacity=20))
        self.cars.append(CozyCar(2, CarKind.BREAK_ROOM, capacity=30))
        self.cars.append(CozyCar(3, CarKind.WORKER_RIDE, capacity=40))

    def status(self):
        return [c.status() for c in self.cars]

# -----------------------------
# Kid observation mission
# -----------------------------

class BuildStage(Enum):
    STAGE_1_START = 1
    STAGE_2_TUNNEL = 2
    STAGE_3_MORE_TUNNELS = 3
    STAGE_4_COZY_TRAIN_READY = 4
    STAGE_5_FULL_SYSTEM_TEST = 5

class KidsObservationMission:
    """
    Kids watch the whole 5-stage build:
    - Robots only on the Moon/Mars
    - No one is in danger
    - Kids see progress snapshots and learn how big projects grow safely
    """
    def __init__(self):
        self.stage = BuildStage.STAGE_1_START
        self.smart_track = SmartTrack()
        self.tunnel_buddy = TunnelBuddy()
        self.cozy_train = CozyTrain()

    def advance_stage_if_ready(self):
        if self.stage == BuildStage.STAGE_1_START and len(self.smart_track.pieces) > 0:
            self.stage = BuildStage.STAGE_2_TUNNEL
        elif self.stage == BuildStage.STAGE_2_TUNNEL and self.tunnel_buddy.progress_m > 1000:
            self.stage = BuildStage.STAGE_3_MORE_TUNNELS
        elif self.stage == BuildStage.STAGE_3_MORE_TUNNELS and self.tunnel_buddy.progress_m > 5000:
            self.stage = BuildStage.STAGE_4_COZY_TRAIN_READY
        elif self.stage == BuildStage.STAGE_4_COZY_TRAIN_READY:
            self.stage = BuildStage.STAGE_5_FULL_SYSTEM_TEST

    def simulate_step(self):
        self.tunnel_buddy.step(self.smart_track)
        track_status = self.smart_track.step()
        self.advance_stage_if_ready()

        snapshot = {
            "stage": self.stage.name,
            "tunnel_progress_m": self.tunnel_buddy.progress_m,
            "track_global_pause": track_status["global_pause"],
            "num_track_pieces": len(self.smart_track.pieces),
            "cozy_train": self.cozy_train.status(),
        }
        return snapshot


if __name__ == "__main__":
    mission = KidsObservationMission()
    for step in range(2000):
        snap = mission.simulate_step()
        if step % 200 == 0:
            print(f"Step {step}")
            print(f"  Stage: {snap['stage']}")
            print(f"  Tunnel progress (m): {snap['tunnel_progress_m']:.1f}")
            print(f"  Track pieces: {snap['num_track_pieces']}")
            print(f"  Track global pause: {snap['track_global_pause']}")
            print(f"  Cozy train: {snap['cozy_train']}")
            print("-" * 40)
