import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================
# INPUT
# =========================

qzv_file = "table.qzv"

# =========================
# Extract sample depth
# =========================

z = zipfile.ZipFile(qzv_file)

sample_file = [f for f in z.namelist() if "sample-frequency-detail.csv" in f][0]

df = pd.read_csv(z.open(sample_file))

sample_depth = df.iloc[:,1]

print("===== Sample depth summary =====")
print(sample_depth.describe())

# =========================
# Rarefaction depth candidates
# =========================

depth_range = np.arange(1000, sample_depth.max(), 500)

retention = []

for d in depth_range:
    retained = (sample_depth >= d).sum()
    retention.append(retained)

retention_df = pd.DataFrame({
    "depth": depth_range,
    "samples_retained": retention
})

print("\n===== Retention table =====")
print(retention_df)

# =========================
# recommended rarefaction depth
# =========================

total_samples = len(sample_depth)

retention_df["retention_rate"] = retention_df["samples_retained"] / total_samples

recommended = retention_df[
    retention_df["retention_rate"] >= 0.9
]

print("\n===== Recommended depths (>=90% samples retained) =====")
print(recommended)

# =========================
# Plot retention curve
# =========================

plt.figure(figsize=(6,4))

plt.plot(retention_df["depth"], retention_df["samples_retained"])

plt.xlabel("Rarefaction depth")
plt.ylabel("Samples retained")

plt.title("Sample Retention Curve")

plt.grid()

plt.savefig("rarefaction_retention_curve.png", dpi=300)

print("\nRetention curve saved as:")
print("rarefaction_retention_curve.png")

# =========================
# Percentile reference
# =========================

print("\n===== Depth percentiles =====")

for p in [5,10,25,50]:
    val = int(sample_depth.quantile(p/100))
    print(f"{p}% percentile depth:", val)
