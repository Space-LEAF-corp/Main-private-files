Jarvondis (refactor)
=====================

This workspace contains a refactored `Jarvondis` module and a small CLI runner.

Quick start
-----------

- Install dependencies (optional, recommended):

```bash
python -m pip install -r requirements.txt
```

- Run interactively:

```bash
python cli.py run --memory-file jarvondis_memory.csv --format csv --tone witty
```

Files added
- `jarvondis/jarvondis.py` — refactored class and persistence
- `cli.py` — small interactive CLI
- `tests/test_jarvondis.py` — unit tests (uses built-in unittest)
- `requirements.txt` — pandas, numpy

Notes
- The `ErebusSync` class is a placeholder. Replace with your real integration.
- Memory can be saved as CSV (default) or JSON using `--format json`.
# Main-private-files
Big or small
taking a step into becoming a professional developer and creating a new type of product for secure and impregnable purposes.