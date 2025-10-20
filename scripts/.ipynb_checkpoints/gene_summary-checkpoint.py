############################################################
# gene_summary.py — AI-generated summaries for top genes
# ----------------------------------------------------------
# Dynamically selects model based on available hardware.
# Saves TSV file compatible with semantic clustering + report.
#
# Author: Oguzhan Begik, PhD
############################################################

import os
import yaml
import torch
import psutil
import pandas as pd
from tqdm import tqdm
from transformers import pipeline

# === Load config ===
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
METRIC = config["metric"]
RESULTS_DIR = config["results_dir"]

in_path = os.path.join(RESULTS_DIR, f"top_PC1_contributor_genes__{CANCER_TYPE}__{METRIC}.tsv")
out_path = os.path.join(RESULTS_DIR, f"PC1_gene_summaries_HF__{CANCER_TYPE}__{METRIC}.tsv")

# === Load genes ===
df = pd.read_csv(in_path, sep="\t")
genes = df["Gene_Symbol"].dropna().unique().tolist()
print(f"🧬 Generating summaries for {len(genes)} genes...")

# === Detect system resources ===
ram_gb = psutil.virtual_memory().total / (1024**3)
device = "cuda" if torch.cuda.is_available() else ("mps" if torch.backends.mps.is_available() else "cpu")
print(f"🖥️ Detected device: {device}, RAM ≈ {ram_gb:.1f} GB")

# === Choose model intelligently ===
if device != "cpu" and ram_gb > 20:
    model_name = "mistralai/Mistral-7B-Instruct-v0.2"
    task = "text-generation"
    print("🚀 Using Mistral-7B for deep summarization...")
elif ram_gb > 12:
    model_name = "google/flan-t5-base"
    task = "text2text-generation"
    print("⚡ Using Flan-T5-Base for balanced summarization...")
else:
    model_name = "google/flan-t5-small"
    task = "text2text-generation"
    print("💨 Using Flan-T5-Small for lightweight summarization...")

# === Initialize model ===
summarizer = pipeline(task, model=model_name, device_map="auto")

# === Generate summaries ===
results = []
for gene in tqdm(genes, desc="Summarizing genes"):
    prompt = f"Summarize the biological role of {gene} in immune biology and cancer."
    try:
        output = summarizer(prompt, max_new_tokens=128, temperature=0.5)[0]
        summary = (
            output["generated_text"]
            if "generated_text" in output
            else output.get("summary_text", "")
        )
    except Exception as e:
        summary = f"Error: {e}"
    results.append({"Gene_Symbol": gene, "Summary": summary.strip()})

# === Save output ===
out_df = pd.DataFrame(results)
out_df.to_csv(out_path, sep="\t", index=False)
print(f"✅ Summaries saved to: {out_path}")