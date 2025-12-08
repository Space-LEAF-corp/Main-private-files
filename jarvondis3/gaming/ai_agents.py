#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from jarvondis3.core.sovereign import Jarvondis3AgentHarness

class GameAgent:
    """
    Example game agent constrained by Jarvondis harness.
    """
    def __init__(self, harness: Jarvondis3AgentHarness, name: str):
        self.harness = harness
        self.name = name

    def act(self, command: str):
        """
        Agent attempts an action via sovereign interface.
        Non-sensitive commands permitted under emergency-mode for operators;
        sensitive commands require admin signature (which agents cannot provide).
        """
        return self.harness.execute(command)
