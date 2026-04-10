# Predictive model schematic: Weather vs Climate (10-year outlook for Ocala, FL)
# Dependencies: numpy, pandas, matplotlib, seaborn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(42)

# -----------------------------
# 1) Configuration (assumptions)
# -----------------------------
years = 10
start_year = 2026
months = years * 12
date_index = pd.date_range(f"{start_year}-01-01", periods=months, freq="MS")

# Baseline climatology assumptions for Ocala-like conditions (illustrative)
# Temperature: baseline monthly mean (°F) with a seasonal cycle roughly 55–92°F
T_baseline_mean_annual = 72.0  # average annual mean
T_season_amp = 15.0            # seasonal amplitude (peak-to-mean)
T_phase = -np.pi/2             # phase to align warmest months to Jul-Aug
T_noise_sigma0 = 5.0           # weather noise std (°F)
T_trend_per_year = 0.4         # warming trend (°F per year)
T_var_inflate_per_year = 0.02  # variance inflation per year (extremes)

# Precipitation: baseline monthly mean (inches), seasonal wet summer pattern
P_baseline_mean_annual = 4.5   # monthly mean precipitation (inches)
P_season_amp = 1.5             # seasonal amplitude
P_phase = 0.0
P_noise_sigma0 = 1.2           # noise std (inches)
P_trend_per_year = 0.05        # change in monthly mean (inches per year)
P_var_inflate_per_year = 0.03  # variance inflation per year

# Mild AR(1) for “weather” autocorrelation
phi_T = 0.3
phi_P = 0.25

# -----------------------------
# 2) Build seasonal cycles
# -----------------------------
m = np.arange(months)
month_in_year = (m % 12)

# Temperature seasonal signal
T_season = T_season_amp * np.sin(2 * np.pi * month_in_year / 12 + T_phase)
T_mean0 = T_baseline_mean_annual + T_season

# Precipitation seasonal signal
P_season = P_season_amp * np.sin(2 * np.pi * month_in_year / 12 + P_phase)
P_mean0 = P_baseline_mean_annual + P_season
P_mean0 = np.clip(P_mean0, 0.1, None)  # no negative precip

# -----------------------------
# 3) Climate trend + variance drift
# -----------------------------
t_years = m / 12.0

T_mean_trend = T_mean0 + T_trend_per_year * t_years
P_mean_trend = P_mean0 + P_trend_per_year * t_years

T_sigma = T_noise_sigma0 * np.sqrt(1 + T_var_inflate_per_year * t_years)
P_sigma = P_noise_sigma0 * np.sqrt(1 + P_var_inflate_per_year * t_years)

# -----------------------------
# 4) Weather realizations (AR(1))
# -----------------------------
def ar1_noise(sigma_series, phi):
    eps = np.zeros_like(sigma_series)
    for i in range(len(sigma_series)):
        innov = np.random.normal(0, sigma_series[i])
        eps[i] = (phi * eps[i-1] + innov) if i > 0 else innov
    return eps

T_weather = ar1_noise(T_sigma, phi_T)
P_weather = ar1_noise(P_sigma, phi_P)

T_obs = T_mean_trend + T_weather
P_obs = np.clip(P_mean_trend + P_weather, 0.0, None)

# -----------------------------
# 5) Uncertainty bands (climate view)
# -----------------------------
# Assume normal uncertainty around climate mean due to internal variability
T_upper = T_mean_trend + 1.64 * T_sigma  # ~95th percentile
T_lower = T_mean_trend - 1.64 * T_sigma  # ~5th percentile

P_upper = np.clip(P_mean_trend + 1.64 * P_sigma, 0.0, None)
P_lower = np.clip(P_mean_trend - 1.64 * P_sigma, 0.0, None)

# -----------------------------
# 6) Summaries
# -----------------------------
df = pd.DataFrame({
    "date": date_index,
    "month": month_in_year,
    "year": start_year + (m // 12),
    "T_mean_trend": T_mean_trend,
    "T_obs": T_obs,
    "T_lower": T_lower,
    "T_upper": T_upper,
    "P_mean_trend": P_mean_trend,
    "P_obs": P_obs,
    "P_lower": P_lower,
    "P_upper": P_upper
})

annual = df.groupby("year").agg({
    "T_mean_trend": "mean",
    "T_obs": "mean",
    "P_mean_trend": "sum",
    "P_obs": "sum"
}).rename(columns={
    "T_mean_trend": "AnnualMeanTemp_F",
    "T_obs": "AnnualMeanTempObserved_F",
    "P_mean_trend": "AnnualTotalPrecip_in",
    "P_obs": "AnnualTotalPrecipObserved_in"
}).reset_index()

warming_per_decade = (annual["AnnualMeanTemp_F"].iloc[-1] - annual["AnnualMeanTemp_F"].iloc[0])
precip_change_per_decade = (annual["AnnualTotalPrecip_in"].iloc[-1] - annual["AnnualTotalPrecip_in"].iloc[0])

# -----------------------------
# 7) Plots
# -----------------------------
sns.set(style="whitegrid", context="talk")

fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# Temperature
axes[0].plot(df["date"], df["T_obs"], color="tab:red", alpha=0.6, label="Weather realization (daily-like monthly)")
axes[0].plot(df["date"], df["T_mean_trend"], color="darkred", linewidth=2, label="Climate mean trend")
axes[0].fill_between(df["date"], df["T_lower"], df["T_upper"], color="lightcoral", alpha=0.3, label="Uncertainty band (5–95%)")
axes[0].set_ylabel("Temperature (°F)")
axes[0].set_title("Temperature: Weather vs Climate (Ocala, FL, next 10 years)")
axes[0].legend(loc="upper left")

# Precipitation
axes[1].plot(df["date"], df["P_obs"], color="tabblue", alpha=0.6, label="Weather realization")
axes[1].plot(df["date"], df["P_mean_trend"], color="navy", linewidth=2, label="Climate mean trend")
axes[1].fill_between(df["date"], df["P_lower"], df["P_upper"], color="lightblue", alpha=0.3, label="Uncertainty band (5–95%)")
axes[1].set_ylabel("Monthly precipitation (in)")
axes[1].set_title("Precipitation: Weather vs Climate (Ocala, FL, next 10 years)")
axes[1].legend(loc="upper left")

plt.xlabel("Date")
plt.tight_layout()
plt.show()

# -----------------------------
# 8) Print concise metrics
# -----------------------------
print(f"Warming over 10 years (mean climatological): {warming_per_decade:.2f} °F")
print(f"Change in annual precipitation over 10 years (climatological): {precip_change_per_decade:.2f} inches")
print("\nAnnual summary (first and last years):")
print(annual.head(1))
print(annual.tail(1))
