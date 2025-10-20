import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import yaml

# -------- Load Config --------
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
METRIC = config["metric"]
THRESHOLD = config["outlier_z_threshold"]

RESULTS_DIR = config["results_dir"]

PLOT_DIR = os.path.join(RESULTS_DIR, "qc_plots")
os.makedirs(PLOT_DIR, exist_ok=True)

# ---------- Filenames with METRIC -----------
EXPR_PATH = os.path.join(RESULTS_DIR, f"expr_topimmune__{CANCER_TYPE}__{METRIC}.parquet")
CLIN_PATH = os.path.join(RESULTS_DIR, f"clinical__{CANCER_TYPE}.tsv")

OUT_EXPR_CLEAN = os.path.join(RESULTS_DIR, f"expr_topimmune_cleaned__{CANCER_TYPE}__{METRIC}.parquet")
OUT_CLIN_CLEAN = os.path.join(RESULTS_DIR, f"clinical_cleaned__{CANCER_TYPE}__{METRIC}.tsv")
OUTLIER_REPORT = os.path.join(RESULTS_DIR, f"outliers__{CANCER_TYPE}__{METRIC}.tsv")
PCA_PLOT_PATH = os.path.join(PLOT_DIR, f"pca_outlier_detection__{CANCER_TYPE}__{METRIC}.pdf")

# ----------------- LOAD DATA -----------------
print(f"Loading expression: {EXPR_PATH}")
expr = pd.read_parquet(EXPR_PATH)
clinical = pd.read_csv(CLIN_PATH, sep="\t")

# ----------------- PCA + Z-SCORE -----------------
print("Running PCA on expression data...")
expr_T = expr.T  # shape: samples x genes
expr_scaled = StandardScaler().fit_transform(expr_T)
pca = PCA(n_components=2)
pc = pca.fit_transform(expr_scaled)

pc_df = pd.DataFrame(pc, columns=["PC1", "PC2"], index=expr_T.index)

# Z-score to detect outliers
z_scores = np.abs((pc_df - pc_df.mean()) / pc_df.std())
outliers = pc_df[(z_scores["PC1"] > THRESHOLD) | (z_scores["PC2"] > THRESHOLD)]
outlier_ids = outliers.index.tolist()

# Save outlier report
outliers.to_csv(OUTLIER_REPORT, sep="\t")
print(f"Outliers detected: {len(outlier_ids)} → saved to {OUTLIER_REPORT}")

# ----------------- PLOT PCA -----------------
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pc_df, x="PC1", y="PC2", alpha=0.7, label="Samples")
if len(outlier_ids) > 0:
    sns.scatterplot(data=pc_df.loc[outlier_ids], x="PC1", y="PC2", color="red", label="Outliers")
plt.title(f"PCA Outlier Detection ({CANCER_TYPE} | {METRIC})")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend()
plt.tight_layout()

# --- Dual export: PDF + PNG ---
base_path = PCA_PLOT_PATH.replace(".pdf", "")
plt.savefig(base_path + ".pdf", dpi=300)
plt.savefig(base_path + ".png", dpi=300)
plt.close()

print(f"PCA plot saved to: {base_path}.pdf and .png")

# ----------------- FILTER EXPRESSION + CLINICAL -----------------
print("Filtering out outlier samples...")
expr_cleaned = expr.drop(columns=outlier_ids, errors="ignore")
clinical_cleaned = clinical[~clinical["sample"].isin(outlier_ids)].copy()

# Save cleaned datasets
expr_cleaned.to_parquet(OUT_EXPR_CLEAN)
clinical_cleaned.to_csv(OUT_CLIN_CLEAN, sep="\t", index=False)

print(f"Cleaned expression matrix saved to: {OUT_EXPR_CLEAN}")
print(f"Cleaned clinical data saved to: {OUT_CLIN_CLEAN}")