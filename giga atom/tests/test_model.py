# tests/test_model.py
from giga_atom.simulations import compute_excess_and_lambda, baseline_counts # pyright: ignore[reportAttributeAccessIssue, reportUnknownVariableType]

# Type hints for imported functions (for type checkers)
from typing import Callable, List, Any, Protocol

# Protocol for the expected state object
class StateLike(Protocol):
    excess_fraction: float
    lambda_daily: float

compute_excess_and_lambda: Callable[[Any, float], List[StateLike]]
baseline_counts: Callable[[], list] # pyright: ignore[reportMissingTypeArgument]

def test_baseline_no_excess():
    counts = baseline_counts() # pyright: ignore[reportUnknownVariableType]
    states: List[StateLike] = compute_excess_and_lambda(counts, 3.0) # pyright: ignore[reportUnknownVariableType]
    from typing import cast
    assert all(cast(StateLike, s).excess_fraction == 0.0 for s in states)  # type: ignore[attr-defined]
    assert all(cast(StateLike, s).lambda_daily == 0.0 for s in states)  # type: ignore[attr-defined]
