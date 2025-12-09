"""Refactored Jarvondis module.

Provides: Jarvondis, Personality, ErebusSync (placeholder)

This module aims to be a small, testable rework of the original
`Jarvondis 2.0` script with safer save/load and a simple API usable
from a CLI or tests.
"""
from __future__ import annotations

import csv
import json
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime
from typing import List, Dict, Optional

try:
    import pandas as pd
except Exception:  # pragma: no cover - fallback behavior exercised in tests if needed
    pd = None


class ErebusSync:
    """Placeholder bridge to whatever backend/service Jarvondis queries.

    Replace or subclass for real integrations.
    """

    def query(self, input_str: str) -> str:
        return f"Echoing back: {input_str}"


@dataclass
class Personality:
    tone: str = "neutral"

    def stylize(self, response: str) -> str:
        if self.tone == "witty":
            return f"{response} 😉"
        if self.tone == "formal":
            return f"{response}. I hope that satisfies your query."
        if self.tone == "mythic":
            return f"⚔️ {response} — inscribed in the Captain’s Log."
        if self.tone == "playful":
            return f"{response} 🎮✨"
        return response


class Jarvondis:
    def __init__(self, personality: Optional[Personality] = None, memory_file: str = "jarvondis_memory.csv", memory_format: str = "csv"):
        self.eremus_sync = ErebusSync()
        self.personality = personality or Personality(tone="witty")
        self.memory_file = memory_file
        self.memory_format = memory_format.lower()
        self._memory: List[Dict[str, str]] = []
        self._load_memory()

    def initialize(self) -> None:
        print("🟢 Jarvondis online. Ready to assist.")

    def learn(self, input_str: str) -> str:
        try:
            response = self.eremus_sync.query(input_str)
            styled = self.personality.stylize(response)
            entry = {
                "timestamp": datetime.now().isoformat(),
                "input": input_str,
                "response": styled,
                "tone": self.personality.tone,
            }
            self._memory.append(entry)
            return styled
        except Exception as e:
            return f"⚠️ Error processing input: {e}"

    def respond(self, input_str: str) -> str:
        return self.learn(input_str)
    
    def remember(self, input_str: str, response: str, topic: Optional[str] = None) -> None:
        """Store a custom memory entry with input, response, and optional topic.
        
        Args:
            input_str: The input text
            response: The response text
            topic: Optional topic/category for the memory
        """
        entry = {
            "timestamp": datetime.now().isoformat(),
            "input": input_str,
            "response": response,
            "tone": self.personality.tone,
        }
        if topic:
            entry["topic"] = topic
        self._memory.append(entry)

    # Persistence
    def save_memory(self, filename: Optional[str] = None) -> None:
        filename = filename or self.memory_file
        if self.memory_format == "json" or filename.endswith(".json"):
            self._save_json(filename)
        else:
            self._save_csv(filename)

    def _save_json(self, filename: str) -> None:
        tmp = tempfile.NamedTemporaryFile("w", delete=False, dir=os.path.dirname(os.path.abspath(filename)) or None)
        try:
            with open(tmp.name, "w", encoding="utf-8") as f:
                json.dump(self._memory, f, ensure_ascii=False, indent=2)
            os.replace(tmp.name, filename)
        finally:
            if os.path.exists(tmp.name):
                try:
                    os.remove(tmp.name)
                except Exception:
                    pass

    def _save_csv(self, filename: str) -> None:
        # atomic write
        tmp = tempfile.NamedTemporaryFile("w", delete=False, newline="", dir=os.path.dirname(os.path.abspath(filename)) or None)
        try:
            fieldnames = ["timestamp", "input", "response", "tone"]
            with open(tmp.name, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for row in self._memory:
                    writer.writerow({k: row.get(k, "") for k in fieldnames})
            os.replace(tmp.name, filename)
        finally:
            if os.path.exists(tmp.name):
                try:
                    os.remove(tmp.name)
                except Exception:
                    pass

    def _load_memory(self, filename: Optional[str] = None) -> None:
        filename = filename or self.memory_file
        if not os.path.exists(filename):
            self._memory = []
            return
        try:
            if self.memory_format == "json" or filename.endswith(".json"):
                with open(filename, "r", encoding="utf-8") as f:
                    self._memory = json.load(f)
            else:
                # prefer pandas if available for robust reading
                if pd is not None:
                    df = pd.read_csv(filename)
                    self._memory = df.fillna("").to_dict(orient="records")
                else:
                    with open(filename, "r", encoding="utf-8", newline="") as f:
                        reader = csv.DictReader(f)
                        self._memory = [dict(row) for row in reader]
        except Exception:
            # on error, reset memory but don't raise to keep interactive loop robust
            self._memory = []

    @property
    def memory(self) -> List[Dict[str, str]]:
        return list(self._memory)


__all__ = ["Jarvondis", "Personality", "ErebusSync"]
