import numpy as np
from scipy.constants import G, pi

# Constants
mu_sun = 1.327e20      # m^3/s^2
AU = 1.496e11          # meters
mu_earth = 3.986e14    # m^3/s^2

def hohmann_dv(r1, r2, mu=mu_sun):
    """General Hohmann transfer Δv for circular→circular orbits."""
    if r1 <= 0 or r2 <= 0:
        raise ValueError("Orbital radii must be positive.")

    a_trans = 0.5 * (r1 + r2)

    v1 = np.sqrt(mu / r1)
    v_trans1 = np.sqrt(mu * (2/r1 - 1/a_trans))
    dv1 = abs(v_trans1 - v1)

    v2 = np.sqrt(mu / r2)
    v_trans2 = np.sqrt(mu * (2/r2 - 1/a_trans))
    dv2 = abs(v_trans2 - v2)

    return dv1 + dv2


def nea_rendezvous_dv_approx(a_nea_au=1.05, e_nea=0.08, inc_deg=5.0):
    """Approximate Δv from LEO to NEA rendezvous (km/s)."""
    r_earth = AU
    r_nea_eq = a_nea_au * AU * np.sqrt(1 - e_nea**2)

    dv_hoh = hohmann_dv(r_earth, r_nea_eq) / 1000  # km/s

    # Plane change at NEA orbital speed
    v_nea = np.sqrt(mu_sun / r_nea_eq) / 1000
    plane_change_dv = 2 * v_nea * np.sin(np.deg2rad(inc_deg) / 2)

    dv_total = dv_hoh + plane_change_dv + 0.5  # capture fudge

    print(f"Approx Hohmann dv: {dv_hoh:.2f} km/s")
    print(f"Plane change: {plane_change_dv:.2f} km/s")
    print(f"Total est rendezvous dv: {dv_total:.2f} km/s")

    return dv_total


def swarm_harvest_sim(num_bots=10, bot_mass_kg=100, dv_available_km_s=5.0,
                      isp_ion=3000, harvest_rate_kg_day=5.0):
    """Swarm propellant + harvest throughput model."""
    ve = isp_ion * 9.81 / 1000  # km/s

    # Rocket equation (no thrust fraction inside exponent)
    mass_ratio = np.exp(dv_available_km_s / ve)
    propellant_per_bot = bot_mass_kg * (mass_ratio - 1)
    total_prop = num_bots * propellant_per_bot

    total_harvest_kg = num_bots * harvest_rate_kg_day * 365

    print(f"\nSwarm Sim ({num_bots} bots):")
    print(f"  Propellant per bot: {propellant_per_bot:.0f} kg")
    print(f"  Total propellant: {total_prop:.0f} kg")
    print(f"  Annual harvest: {total_harvest_kg/1000:.1f} tons")
    print(f"  Go/No-Go: Viable if dv < {dv_available_km_s} km/s and on-site refuel available.")

    return {
        "propellant_per_bot": propellant_per_bot,
        "total_propellant": total_prop,
        "annual_harvest_kg": total_harvest_kg
    }


# Example run
print("Delta-V Calc for low-dv NEA")
dv_example = nea_rendezvous_dv_approx(a_nea_au=1.02, e_nea=0.05, inc_deg=2.0)
swarm_harvest_sim(num_bots=20, dv_available_km_s=dv_example)