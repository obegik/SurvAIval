import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sksurv.util import Surv
from sksurv.nonparametric import kaplan_meier_estimator
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import yaml
import mygene


# -------- Load config --------
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
METRIC = config["metric"]
RESULTS_DIR = config["results_dir"]



EXPR_PATH = os.path.join(RESULTS_DIR, f"expr_topimmune_cleaned__{CANCER_TYPE}__{METRIC}.parquet")
CLIN_PATH = os.path.join(RESULTS_DIR, f"clinical_cleaned__{CANCER_TYPE}__{METRIC}.tsv")
KM_PLOT_PATH = os.path.join(RESULTS_DIR, f"km_PC1_split__{CANCER_TYPE}__{METRIC}.pdf")

# -------- Load data --------
print("Loading data...")
expr = pd.read_parquet(EXPR_PATH)
clinical = pd.read_csv(CLIN_PATH, sep="\t")

# Transpose expression: samples as rows
expr_T = expr.T
expr_T.index.name = "sample"
expr_T.reset_index(inplace=True)

# Merge with clinical
merged = expr_T.merge(clinical, on="sample", how="inner")
print(f"Merged samples: {merged.shape[0]}")

# -------- PCA --------
gene_cols = expr.index.tolist()
X = merged[gene_cols]
X_scaled = StandardScaler().fit_transform(X)
pc1 = PCA(n_components=1).fit_transform(X_scaled).flatten()

# Split by PC1 median
pc1_split = pd.Series(pc1, index=merged["sample"])
median_score = pc1_split.median()
risk_group = np.where(pc1_split > median_score, "High PC1", "Low PC1")

# -------- KM Plot --------
times = merged["OS.time"]
events = merged["OS"] == 1

# Ensure no missing survival data


fig, ax = plt.subplots(figsize=(8, 6))

for label, color in zip(["Low PC1", "High PC1"], ["#1f77b4", "#ff7f0e"]):
    # build logical mask for this group, excluding missing values
    mask = (risk_group == label) & ~pd.isna(events) & ~pd.isna(times)
    
    # subset for survival times and events
    t, s = kaplan_meier_estimator(events[mask], times[mask])
    
    n_samples = mask.sum()
    label_with_n = f"{label} (n={n_samples})"
    
    ax.step(t, s, where="post", label=label_with_n, color=color)
    
    # confidence region
    s_lower = np.clip(s * 0.95, 0, 1)
    s_upper = np.clip(s * 1.05, 0, 1)
    ax.fill_between(t, s_lower, s_upper, step="post", alpha=0.2, color=color)

# -------- Log-rank test --------
result = logrank_test(
    times[risk_group == "High PC1"], times[risk_group == "Low PC1"],
    event_observed_A=events[risk_group == "High PC1"],
    event_observed_B=events[risk_group == "Low PC1"]
)
pval = result.p_value


ax.set_title
# -------- Finalize plot --------
pval_str = f"{pval:.4f}".rstrip("0").rstrip(".") if pval >= 0.0001 else "<0.0001"
ax.set_title(f"KM Stratification by PC1 (p = {pval_str})\n{CANCER_TYPE}, {METRIC}")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Survival Probability")
ax.legend()
plt.tight_layout()

# --- Dual export (PDF + PNG) ---
base_path = KM_PLOT_PATH.replace(".pdf", "")
plt.savefig(base_path + ".pdf", dpi=300)
plt.savefig(base_path + ".png", dpi=300)
plt.close()

print(f"KM plot saved to: {base_path}.pdf and .png")





# -------- Assumptions --------
# - expr = genes x samples (DataFrame)
# - expr_T = samples x genes (already transposed)
# - RESULTS_DIR, CANCER_TYPE, METRIC defined
# - gene_cols = list(expr.index)

# -------- PCA Computation --------
print(" Running PCA on scaled expression...")
X_scaled = StandardScaler().fit_transform(expr_T[gene_cols])
pca = PCA(n_components=1)
pca.fit(X_scaled)

# Get PC1 loadings
pc1_loadings = pd.Series(pca.components_[0], index=gene_cols)
pc1_sorted = pc1_loadings.sort_values(key=np.abs, ascending=False)

# Save raw PC1 loadings
pc1_sorted.to_csv(
    os.path.join(RESULTS_DIR, f"PC1_gene_loadings__{CANCER_TYPE}__{METRIC}.tsv"),
    sep="\t"
)

# -------- Gene Annotation --------
print("Annotating gene symbols with MyGene.info...")
mg = mygene.MyGeneInfo()
query_result = mg.querymany(
    pc1_sorted.index.tolist(),
    scopes="ensembl.gene",
    fields="symbol,name",
    species="human"
)

annot_df = pd.DataFrame(query_result).set_index("query")
pc1_annotated = pc1_sorted.rename("PC1_loading").to_frame().join(
    annot_df[["symbol", "name"]], how="left"
)

# Save annotated loadings
pc1_annotated.to_csv(
    os.path.join(RESULTS_DIR, f"PC1_gene_loadings_annotated__{CANCER_TYPE}__{METRIC}.tsv"),
    sep="\t"
)
print(" Annotated PC1 gene loadings saved.")



# -------- Plot Top PC1 Contributors with Gene Symbols --------
top_n = 15

# Select top positive and negative loadings
top_positive = pc1_annotated.sort_values("PC1_loading", ascending=False).head(top_n)
top_negative = pc1_annotated.sort_values("PC1_loading", ascending=True).head(top_n)

# Combine and sort properly by PC1_loading value
top_contributors = pd.concat([top_negative, top_positive]).sort_values("PC1_loading")

# Fallback to Ensembl ID if symbol missing
top_contributors["label"] = top_contributors["symbol"].combine_first(top_contributors.index.to_series())

# -------- Export Top PC1 Contributor Genes --------
top_contributors_export = top_contributors.reset_index()[["index", "symbol", "name", "PC1_loading"]]
top_contributors_export.columns = ["Ensembl_ID", "Gene_Symbol", "Gene_Name", "PC1_Loading"]

# Save file
export_path = os.path.join(RESULTS_DIR, f"top_PC1_contributor_genes__{CANCER_TYPE}__{METRIC}.tsv")
top_contributors_export.to_csv(export_path, sep="\t", index=False)
print(f"Top contributor gene table exported to: {export_path}")

# Define colors: orange for positive, blue for negative
colors = ["#ff7f0e" if x > 0 else "#1f77b4" for x in top_contributors["PC1_loading"]]

# Plot horizontal bar chart
plt.figure(figsize=(10, 6))
plt.barh(
    y=top_contributors["label"],
    width=top_contributors["PC1_loading"],
    color=colors
)
plt.title("Top PC1 Gene Contributors")
plt.xlabel("PC1 Loading")
plt.tight_layout()


# --- Dual export (PDF + PNG) ---
output_path = os.path.join(RESULTS_DIR, f"PC1_top_gene_contributors_annotated__{CANCER_TYPE}__{METRIC}")
plt.savefig(output_path + ".pdf", dpi=300)
plt.savefig(output_path + ".png", dpi=300)
plt.close()

print(f"Sorted PC1 contributor plot saved to: {output_path}.pdf and .png")


