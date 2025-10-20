import os
import pandas as pd
import yaml

# ----------- LOAD CONFIG -----------
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
METRIC = config["metric"]
DATA_DIR = config["data_dir"]
RESULTS_DIR = config["results_dir"]


# ---------- PATHS ----------
EXPR_PATH = os.path.join(RESULTS_DIR, f"expr__{CANCER_TYPE}.parquet")
TOP_GENES_PATH = os.path.join(RESULTS_DIR, f"top_immune_genes__{CANCER_TYPE}__{METRIC}.tsv")

OUTPUT_EXPR_PATH = os.path.join(RESULTS_DIR, f"expr_topimmune__{CANCER_TYPE}__{METRIC}.parquet")
OUTPUT_ANNOT_PATH = os.path.join(RESULTS_DIR, f"expr_topimmune_annotated__{CANCER_TYPE}__{METRIC}.parquet")

# ---------- LOAD DATA ----------
print(f"Loading expression matrix from: {EXPR_PATH}")
expr = pd.read_parquet(EXPR_PATH)
expr.index = expr.index.str.replace(r"\.\d+$", "", regex=True)  # Remove Ensembl version suffix

print(f"Loading immune genes from: {TOP_GENES_PATH}")
top_immune_genes = pd.read_csv(TOP_GENES_PATH, sep="\t")

# ---------- KEEP ONLY RELEVANT COLUMNS ----------
cols_to_keep = ["ensembl_id", "Symbol", METRIC]
top_immune_genes = top_immune_genes[cols_to_keep].dropna(subset=["ensembl_id"])

# ---------- SUBSET GENES ----------
ensembl_ids = top_immune_genes["ensembl_id"].unique()
expr_subset = expr.loc[expr.index.isin(ensembl_ids)].copy()

print(f"Matched immune genes in expression matrix: {expr_subset.shape[0]} / {len(ensembl_ids)}")

# ---------- SAVE MATRIX ONLY ----------
expr_subset.to_parquet(OUTPUT_EXPR_PATH)
print(f"Saved subsetted expression matrix to: {OUTPUT_EXPR_PATH}")

# ---------- SAVE ANNOTATED VERSION ----------
expr_annotated = top_immune_genes.set_index("ensembl_id").join(expr_subset, how="inner")
expr_annotated.to_parquet(OUTPUT_ANNOT_PATH)
print(f"Saved annotated expression matrix to: {OUTPUT_ANNOT_PATH}")