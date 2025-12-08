#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from jarvondis3.core.sovereign import milestone_banner

def level_up_ritual(level: int) -> str:
    return f"{milestone_banner('COMMAND_EXECUTED')} Player reached level {level}!"

def boss_defeated_ritual(boss_name: str) -> str:
    return f"{milestone_banner('DIAGNOSTICS_COMPLETE')} Boss {boss_name} defeated under sovereign watch."
