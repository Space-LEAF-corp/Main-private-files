"""
utils.py — Plotting and CSV export utilities
Space Leaf Corp — Internal Use Only
"""

import pandas as pd  # type: ignore
import matplotlib.pyplot as plt # pyright: ignore[reportMissingModuleSource]

def export_csv(df: pd.DataFrame, filename: str):
    df.to_csv(filename, index=False) # pyright: ignore[reportUnknownMemberType]

def plot_series(
    t: list[float] | pd.Series, # pyright: ignore[reportGeneralTypeIssues]
    series_dict: dict[str, pd.Series | list[float]], # pyright: ignore[reportGeneralTypeIssues]
    title: str,
    xlabel: str = "Days",
    ylabel: str = "Value",
    filename: str | None = None # pyright: ignore[reportGeneralTypeIssues]
) -> None:
    plt.figure(figsize=(10,5)) # pyright: ignore[reportUnknownMemberType]
    for label, arr in series_dict.items():
        plt.plot(t, arr, label=str(label)) # pyright: ignore[reportUnknownMemberType]
    plt.title(title) # pyright: ignore[reportUnknownMemberType]
    plt.xlabel(xlabel) # pyright: ignore[reportUnknownMemberType]
    plt.ylabel(ylabel) # pyright: ignore[reportUnknownMemberType]
    plt.grid(visible=True) # pyright: ignore[reportUnknownMemberType]
    plt.legend() # pyright: ignore[reportUnknownMemberType]
    if filename:
        plt.savefig(filename, dpi=150) # pyright: ignore[reportUnknownMemberType]
    plt.show() # pyright: ignore[reportUnknownMemberType]
