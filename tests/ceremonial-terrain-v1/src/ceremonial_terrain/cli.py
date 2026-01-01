import argparse
from .config import SimulationConfig
from .simulation import run_365_day_simulation, validate_forward_backward

def main():
    parser = argparse.ArgumentParser(description="Ceremonial Terrain v1 – Map Generator Prototype")
    subparsers = parser.add_subparsers(dest="command", required=True)

    sim_parser = subparsers.add_parser("run-sim", help="Run a 365-day simulation for a seed")
    sim_parser.add_argument("--seed", type=int, default=42)

    val_parser = subparsers.add_parser("validate-sim", help="Validate determinism for a seed")
    val_parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    if args.command == "run-sim":
        config = SimulationConfig(seed=args.seed)
        snapshots = run_365_day_simulation(config)
        print(f"Simulation complete. Stored {len(snapshots)} snapshots.")
        last_day = max(snapshots.keys())
        last_state = snapshots[last_day]
        print(f"Last day: {last_day}, volcano events: {len(last_state.volcano_events)}")

    elif args.command == "validate-sim":
        config = SimulationConfig(seed=args.seed)
        ok = validate_forward_backward(config)
        if ok:
            print("Validation PASSED: forward and replay states match.")
        else:
            print("Validation FAILED: mismatch between forward and replay states.")
