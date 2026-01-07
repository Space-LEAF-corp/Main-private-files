from dataclasses import dataclass, asdict
from typing import List, Dict
import random

# ----- Core Data Models -----

@dataclass
class SymbolEvent:
    t: float  # timestamp in seconds
    symbol: str
    meaning: str

@dataclass
class SubtitleEvent:
    t: float
    text_en: str

@dataclass
class HapticEvent:
    t: float
    pattern: str  # e.g. "short-short-long"
    intensity: float  # 0.0 - 1.0

@dataclass
class AudioEvent:
    t: float
    cue_type: str  # e.g. "narration", "motif", "transition"
    content: str   # text for narration, label for motif

@dataclass
class NarrativePacket:
    seed: str
    duration: float
    symbols: List[SymbolEvent]
    subtitles: List[SubtitleEvent]
    haptics: List[HapticEvent]
    audio: List[AudioEvent]


# ----- Symbolic Engine -----

class SymbolicEngine:
    def __init__(self):
        # Simple starter lexicon – can be expanded
        self.lexicon = {
            "threshold": "A moment of choice or crossing",
            "river": "Flow of time or emotion",
            "star": "Guidance or distant goal",
            "branch": "Diverging paths or options",
            "helix": "Growth, learning, or lineage",
        }

    def generate_sequence(self, duration: float) -> List[SymbolEvent]:
        events = []
        t = 0.0
        step = duration / 5  # 5 symbolic beats for now
        for _ in range(5):
            symbol = random.choice(list(self.lexicon.keys()))
            events.append(SymbolEvent(
                t=round(t, 2),
                symbol=symbol,
                meaning=self.lexicon[symbol]
            ))
            t += step
        return events


# ----- Accessibility Engines -----

class SubtitleEngine:
    def generate(self, symbols: List[SymbolEvent]) -> List[SubtitleEvent]:
        subtitles = []
        for ev in symbols:
            text = f"{ev.meaning}."
            subtitles.append(SubtitleEvent(t=ev.t, text_en=text))
        return subtitles


class HapticEngine:
    def generate(self, symbols: List[SymbolEvent]) -> List[HapticEvent]:
        pattern_map = {
            "threshold": "short-long",
            "river": "wave-wave",
            "star": "short-short-short",
            "branch": "short-short-long",
            "helix": "ramp-up",
        }
        haptics = []
        for ev in symbols:
            pattern = pattern_map.get(ev.symbol, "short")
            intensity = 0.6 if ev.symbol in ["threshold", "branch"] else 0.4
            haptics.append(HapticEvent(
                t=ev.t,
                pattern=pattern,
                intensity=intensity
            ))
        return haptics


class AudioEngine:
    def generate(self, seed: str, symbols: List[SymbolEvent]) -> List[AudioEvent]:
        audio = []
        # Opening narration
        audio.append(AudioEvent(
            t=0.0,
            cue_type="narration",
            content=f"Here we go. This is the story of {seed} as it unfolds."
        ))
        # Symbol-linked cues
        for ev in symbols:
            audio.append(AudioEvent(
                t=ev.t,
                cue_type="narration",
                content=f"{ev.meaning}."
            ))
        # Closing line
        audio.append(AudioEvent(
            t=symbols[-1].t + 1.0,
            cue_type="narration",
            content="The moment settles. The lesson remains."
        ))
        return audio


# ----- Generator Orchestrator -----

class NarrativeGenerator:
    def __init__(self):
        self.symbolic = SymbolicEngine()
        self.subtitles_engine = SubtitleEngine()
        self.haptics_engine = HapticEngine()
        self.audio_engine = AudioEngine()

    def generate_packet(self, seed: str, duration: float = 30.0) -> Dict:
        symbols = self.symbolic.generate_sequence(duration)
        subtitles = self.subtitles_engine.generate(symbols)
        haptics = self.haptics_engine.generate(symbols)
        audio = self.audio_engine.generate(seed, symbols)

        packet = NarrativePacket(
            seed=seed,
            duration=duration,
            symbols=symbols,
            subtitles=subtitles,
            haptics=haptics,
            audio=audio
        )
        return asdict(packet)


# ----- Example usage -----

if __name__ == "__main__":
    gen = NarrativeGenerator()
    packet = gen.generate_packet(seed="Quantum Particle Intelligence Generator", duration=45.0)
    # In a real system, this packet would be handed off to:
    # - a video editor
    # - a TTS engine (RP voice)
    # - a haptics driver
    # - a subtitle renderer
    print(packet)
