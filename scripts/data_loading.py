import os
import pandas as pd
import urllib.request
import yaml


# ---------- LOAD CONFIG ----------
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
DATA_DIR = config["data_dir"]
RESULTS_DIR = config["results_dir"]
#EXPR_URL = config["expression_url"]
#CLIN_URL = config["clinical_url"]

#EXPR_FILE = os.path.basename(EXPR_URL)
#CLIN_FILE = os.path.basename(CLIN_URL)
#EXPR_PATH = os.path.join(DATA_DIR, EXPR_FILE)
#CLIN_PATH = os.path.join(DATA_DIR, CLIN_FILE)

# ---------- FIXED FILENAMES ----------
EXPR_PATH = os.path.join(DATA_DIR, "tcga_RSEM_gene_tpm")
CLIN_PATH = os.path.join(DATA_DIR, "Survival_SupplementalTable_S1_20171025_xena_sp")

if not os.path.exists(EXPR_PATH):
    raise FileNotFoundError(f"Expression file not found: {EXPR_PATH}")
if not os.path.exists(CLIN_PATH):
    raise FileNotFoundError(f"Clinical file not found: {CLIN_PATH}")

print(f"Using local expression file: {EXPR_PATH}")
print(f"Using local clinical file: {CLIN_PATH}")


# ---------- CREATE DIRS ----------
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ---------- DOWNLOAD IF MISSING ----------
#def download_if_missing(url, path):
#    if not os.path.exists(path):
#        print(f"Downloading {os.path.basename(path)}...")
#        urllib.request.urlretrieve(url, path)
#    else:
#        print(f"{os.path.basename(path)} already exists.")

#download_if_missing(EXPR_URL, EXPR_PATH)
#download_if_missing(CLIN_URL, CLIN_PATH)

# ---------- LOAD DATA ----------
print("Loading expression matrix...")
expr = pd.read_csv(EXPR_PATH, sep="\t", compression="infer", index_col=0)

print("Loading clinical metadata...")
clinical = pd.read_csv(CLIN_PATH, sep="\t")

print(f"Expression shape: {expr.shape}")
print(f"Clinical columns: {clinical.columns.tolist()}")

# ---------- FILTER TO SELECTED CANCER TYPE ----------
# Filter clinical table
clinical_cancer = clinical[clinical["cancer type abbreviation"] == CANCER_TYPE].copy()

# Extract sample barcodes (first 15 characters)
barcodes = clinical_cancer["sample"].str.slice(0, 15)

# Match expression data to those barcodes
expr_filtered = expr.loc[:, expr.columns.str.slice(0, 15).isin(barcodes)].copy()

# Match clinical data to expression samples
clinical_matched = clinical_cancer[clinical_cancer["sample"].str.slice(0, 15).isin(expr_filtered.columns.str.slice(0, 15))].copy()

# ---------- FINAL SHAPES ----------
print(f"Filtered expression shape ({CANCER_TYPE}): {expr_filtered.shape}")
print(f"Matched clinical shape ({CANCER_TYPE}): {clinical_matched.shape}")

# ---------- EXPORT RESULTS ----------
expr_filtered.to_parquet(f"{RESULTS_DIR}/expr__{CANCER_TYPE}.parquet", engine="fastparquet")
clinical_matched.to_csv(f"{RESULTS_DIR}/clinical__{CANCER_TYPE}.tsv", sep="\t", index=False)


