def generate_atomic_fingerprint(sat_id, launch_ts, orbit_class, mission_type):
    seed = hash(f"{sat_id}-{launch_ts}-{orbit_class}-{mission_type}")
    rng = PRNG(seed)

    # Example: pick elements and roles from curated sets
    element = pick_element(rng)
    role = pick_role(rng)
    state_level = rng.randint(1, 9)
    position = generate_position_pattern(rng)

    return {
        "element": element,
        "role": role,
        "state_level": state_level,
        "position_label": position
    }