"""Simple CLI runner for the refactored Jarvondis module.

Usage:
  python cli.py run --memory-file jarvondis_memory.csv --format csv --tone witty
"""
from __future__ import annotations

import argparse
import sys

from jarvondis.jarvondis import Jarvondis, Personality


def interactive_loop(j: Jarvondis) -> None:
    j.initialize()
    try:
        while True:
            user_input = input("You: ")
            if not user_input:
                continue
            if user_input.lower() in ("exit", "quit"):
                print("Jarvondis: Shutting down. 💤")
                j.save_memory()
                break
            print("Jarvondis:", j.respond(user_input))
    except (KeyboardInterrupt, EOFError):
        print("\nJarvondis: Input closed. Shutting down.")
        j.save_memory()


def main(argv=None):
    parser = argparse.ArgumentParser(prog="jarvondis-cli")
    sub = parser.add_subparsers(dest="cmd")

    runp = sub.add_parser("run", help="Run interactive Jarvondis")
    runp.add_argument("--memory-file", default="jarvondis_memory.csv")
    runp.add_argument("--format", default="csv", choices=("csv", "json"))
    runp.add_argument("--tone", default="witty", help="Personality tone (witty, formal, mythic, playful, neutral)")

    args = parser.parse_args(argv)
    if args.cmd == "run":
        p = Personality(tone=args.tone)
        j = Jarvondis(personality=p, memory_file=args.memory_file, memory_format=args.format)
        interactive_loop(j)
    else:
        parser.print_help()


if __name__ == "__main__":
    main(sys.argv[1:])
