"""
utils.py — Plotting and CSV export utilities
Space Leaf Corp — Internal Use Only
"""

import pandas as pd  # type: ignore
import matplotlib.pyplot as plt

def export_csv(df: pd.DataFrame, filename: str):
    df.to_csv(filename, index=False) # pyright: ignore[reportUnknownMemberType]

def plot_series(
    t: list[float] | pd.Series,
    series_dict: dict[str, pd.Series | list[float]],
    title: str,
    xlabel: str = "Days",
    ylabel: str = "Value",
    filename: str | None = None
) -> None:
    plt.figure(figsize=(10,5))
    for label, arr in series_dict.items():
        plt.plot(t, arr, label=str(label))
    plt.title(title) # pyright: ignore[reportUnknownMemberType]
    plt.xlabel(xlabel) # pyright: ignore[reportUnknownMemberType]
    plt.ylabel(ylabel)
    plt.grid(visible=True)
    legend = plt.legend()
    if filename:
        plt.savefig(filename, dpi=150)
    plt.show()
