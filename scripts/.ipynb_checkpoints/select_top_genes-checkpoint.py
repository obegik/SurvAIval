import os
import pandas as pd
import mygene
import yaml

# ----------- LOAD CONFIG -----------
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
CANCER_NAME = config["cancer_name"]
COHORT = "TCGA"
TOP_N = config["top_n"]
METRIC = config["metric"]

DATA_DIR = config["data_dir"]
RESULTS_DIR = config["results_dir"]
SCORES_FILE = os.path.join(DATA_DIR, "all_scores.tsv")
OUTPUT_FILE = os.path.join(RESULTS_DIR, f"top_immune_genes__{CANCER_TYPE}__{METRIC}.tsv")


# ----------- LOAD IMMUNE SCORES -----------
print("Loading immune scores...")
scores = pd.read_csv(SCORES_FILE, sep="\t", low_memory=False)

print(scores)
# Filter to specific cancer, cohort, and RNASeq platform
scores_filtered = scores[
    (scores["Cancer"] == CANCER_NAME) &
    (scores["Cohort"] == COHORT) &
    (scores["Platform"] == "RNASeq")
].copy()

# Drop genes with missing values for selected metric
scores_filtered = scores_filtered.dropna(subset=[METRIC])


print(scores_filtered)

# ----------- SELECT TOP IMMUNE GENES -----------
print(f"Selecting top immune genes using metric: {METRIC}")
if METRIC == "CTL Cor":
    top_half = scores_filtered.sort_values(METRIC, ascending=False).head(TOP_N // 2)
    bottom_half = scores_filtered.sort_values(METRIC, ascending=True).head(TOP_N // 2)
    top_immune_genes = pd.concat([top_half, bottom_half], axis=0)
else:
    top_immune_genes = scores_filtered.sort_values(METRIC, ascending=False).head(TOP_N)

# ----------- MAP SYMBOL → ENSEMBL ID -----------
print("Mapping gene symbols to Ensembl IDs using MyGene...")
mg = mygene.MyGeneInfo()
query = top_immune_genes["Symbol"].dropna().unique().tolist()
gene_info = mg.querymany(query, scopes="symbol", fields="ensembl.gene", species="human")

gene_map = pd.DataFrame(gene_info)

if "notfound" in gene_map.columns:
    gene_map = gene_map[gene_map["notfound"] != True]


# Handle nested 'ensembl' fields
def extract_ensg(x):
    if isinstance(x, list) and isinstance(x[0], dict):
        return x[0].get("gene")
    elif isinstance(x, dict):
        return x.get("gene")
    return None

gene_map["ensembl_id"] = gene_map["ensembl"].apply(extract_ensg)
gene_map = gene_map[["query", "ensembl_id"]].rename(columns={"query": "Symbol"})

# Merge into main table
top_immune_genes = top_immune_genes.merge(gene_map, on="Symbol", how="left")
top_immune_genes = top_immune_genes.dropna(subset=["ensembl_id"])

# ----------- SAVE OUTPUT -----------
os.makedirs(RESULTS_DIR, exist_ok=True)
top_immune_genes.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"Saved top {TOP_N} immune genes based on {METRIC} for {CANCER_TYPE}:")
print(f"{OUTPUT_FILE}")