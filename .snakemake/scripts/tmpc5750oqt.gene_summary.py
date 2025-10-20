
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/opt/miniconda3/envs/pipeline_az/lib/python3.10/site-packages', '/Users/boguzhan/Library/Caches/snakemake/snakemake/source-cache/runtime-cache/tmpfg1l0dk3/file/Users/boguzhan/Desktop/tic_pipeline/scripts', '/Users/boguzhan/Desktop/tic_pipeline/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\x8e\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c:results/top_PC1_contributor_genes__COAD__T_Dysfunction.tsv\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x0ccontributors\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x18)}\x94\x8c\x05_name\x94h\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c6results/PC1_gene_summaries_HF__COAD__T_Dysfunction.tsv\x94a}\x94(h\x0c}\x94\x8c\x07summary\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh)h&ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c0/var/folders/55/00_v6_jj6_vcp7817htf4_zh0000gn/T\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bhZK\x01h\\K\x01h^hWub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x0bcancer_type\x94\x8c\x04COAD\x94\x8c\x0bcancer_name\x94\x8c\nColorectal\x94\x8c\x06metric\x94\x8c\rT_Dysfunction\x94\x8c\x08data_dir\x94\x8c\x04data\x94\x8c\x0bresults_dir\x94\x8c\x07results\x94\x8c\x0bscripts_dir\x94\x8c\x07scripts\x94\x8c\x0bscores_file\x94\x8c\x13data/all_scores.tsv\x94\x8c\x05top_n\x94M\xf4\x01\x8c\x13outlier_z_threshold\x94K\x03\x8c\x0bserpapi_key\x94\x8c@f992519b303a0d170272cbf4fae5a16699d5eb2a1a349ad992ffb6da6f58ce72\x94u\x8c\x04rule\x94\x8c\x0cgene_summary\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c,/Users/boguzhan/Desktop/tic_pipeline/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = True; __real_file__ = __file__; __file__ = '/Users/boguzhan/Desktop/tic_pipeline/scripts/gene_summary.py';
######## snakemake preamble end #########
############################################################
# 06_gene_summary.py
# ----------------------------------------------------------
# Summarize top PC1 genes using Hugging Face LLM API
# - Reads top_PC1_contributor_genes__{CANCER_TYPE}__{METRIC}.tsv
# - Generates short 2–3 sentence summaries per gene
# - Outputs to PC1_gene_summaries_HF__{CANCER_TYPE}__{METRIC}.tsv
#
# Author: Oguzhan Begik, PhD
############################################################

import os
import yaml
import pandas as pd
from tqdm import tqdm
from transformers import pipeline

# -------- Load config --------
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
METRIC = config["metric"]
RESULTS_DIR = config["results_dir"]

# -------- Paths --------
input_path = os.path.join(RESULTS_DIR, f"top_PC1_contributor_genes__{CANCER_TYPE}__{METRIC}.tsv")
output_path = os.path.join(RESULTS_DIR, f"PC1_gene_summaries_HF__{CANCER_TYPE}__{METRIC}.tsv")

# -------- Load input --------
genes_df = pd.read_csv(input_path, sep="\t")
genes = genes_df["Gene_Symbol"].dropna().unique().tolist()

print(f"🧬 Generating Hugging Face summaries for {len(genes)} genes...")

# -------- Initialize HF model --------
summarizer = pipeline("text-generation", model="mistralai/Mistral-7B-Instruct-v0.2", device_map="auto")

# -------- Generate summaries --------
summaries = []
for gene in tqdm(genes, desc="Summarizing genes"):
    prompt = (
        f"Summarize the biological function of the human gene {gene} "
        f"in the context of immune biology and cancer. "
        f"Keep it under 3 sentences and highlight relevance to tumor immunity."
    )
    try:
        response = summarizer(prompt, max_new_tokens=100, temperature=0.3)[0]["generated_text"]
        clean_response = response.replace(prompt, "").strip()
    except Exception as e:
        clean_response = f"Error generating summary: {e}"

    summaries.append({"Gene_Symbol": gene, "Summary": clean_response})

# -------- Save output --------
out_df = pd.DataFrame(summaries)
out_df.to_csv(output_path, sep="\t", index=False)

print(f"✅ Summaries saved to: {output_path}")