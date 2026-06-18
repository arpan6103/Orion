import numpy as np
import matplotlib.pyplot as plt

# Load the CSV
data = np.genfromtxt("src/porkchop.csv", delimiter=",", names=True)

dep = data["departure_day"]
arr = data["arrival_day"]
dv  = data["total_dv"]

# Build a regular grid for contouring
dep_unique = np.unique(dep)
arr_unique = np.unique(arr)
DEP, ARR = np.meshgrid(dep_unique, arr_unique)

# Fill grid with delta-V, NaN where no data
DV = np.full(DEP.shape, np.nan)
dep_idx = {v: i for i, v in enumerate(dep_unique)}
arr_idx = {v: i for i, v in enumerate(arr_unique)}
for d, a, v in zip(dep, arr, dv):
    DV[arr_idx[a], dep_idx[d]] = v

# Cap delta-V for cleaner contours (ignore absurdly expensive transfers)
DV_capped = np.clip(DV, None, 15)

plt.figure(figsize=(10, 8))
levels = np.linspace(np.nanmin(DV_capped), 15, 30)
cs = plt.contourf(DEP, ARR, DV_capped, levels=levels, cmap="jet")
plt.colorbar(cs, label="Total Δv (km/s)")

# Mark the minimum
imin = np.nanargmin(dv)
plt.plot(dep[imin], arr[imin], "w*", markersize=20, markeredgecolor="black")
plt.annotate(f"  min Δv = {dv[imin]:.2f} km/s",
             (dep[imin], arr[imin]), color="black", fontsize=11)

plt.xlabel("Departure day (from epoch)")
plt.ylabel("Arrival day (from epoch)")
plt.title("Earth → Mars Porkchop Plot")
plt.tight_layout()
plt.savefig("porkchop.png", dpi=150)
plt.show()
print(f"Minimum Δv = {dv[imin]:.3f} km/s at departure={dep[imin]:.0f}, arrival={arr[imin]:.0f}")