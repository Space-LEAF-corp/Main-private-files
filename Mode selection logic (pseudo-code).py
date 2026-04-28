S = stability_margin()          # CoG vs support polygon
C = min_joint_capacity_ratio()  # worst joint capacity fraction

if (S < S_low) or (C < C_low):
    mode = BRACED

elif (S < S_mid) or (C < C_mid):
    mode = LIMPING

elif (S < S_high) or (C < C_high):
    mode = GUARDED

else:
    mode = AGILE
