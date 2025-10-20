############################################################
# generate_report.py
# ----------------------------------------------------------
# Final report generator for the TCGA immune survival pipeline
#
# Includes:
# - Executive summary
# - Top PC1 genes table
# - PCA outlier, KM survival, semantic clustering plots
# - Gene summaries (HuggingFace)
#
# Compatible with pipeline scripts:
#   data_loading.py
#   subset_expression.py
#   select_top_genes.py
#   qc_outlier_removal.py
#   pc1_survival.py
#   sementic_clustering.py
#
# Author: Oguzhan Begik, PhD
############################################################

import os
import yaml
from datetime import datetime
import pandas as pd
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors

# -------- Load config --------
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
METRIC = config["metric"]
RESULTS_DIR = config.get("results_dir", "results")

# -------- File paths --------
TOP_GENES = os.path.join(RESULTS_DIR, f"top_PC1_contributor_genes__{CANCER_TYPE}__{METRIC}.tsv")
OUTLIERS = os.path.join(RESULTS_DIR, f"outliers__{CANCER_TYPE}__{METRIC}.tsv")
SUMMARY_PATH = os.path.join(RESULTS_DIR, f"PC1_gene_summaries_HF__{CANCER_TYPE}__{METRIC}.tsv")

PCA_PLOT = os.path.join(RESULTS_DIR, f"qc_plots/pca_outlier_detection__{CANCER_TYPE}__{METRIC}.png")
KM_PLOT = os.path.join(RESULTS_DIR, f"km_PC1_split__{CANCER_TYPE}__{METRIC}.png")
PC1_PLOT = os.path.join(RESULTS_DIR, f"PC1_top_gene_contributors_annotated__{CANCER_TYPE}__{METRIC}.png")
SEMANTIC_PLOT = os.path.join(RESULTS_DIR, f"semantic_clusters_stable__{CANCER_TYPE}__{METRIC}.png")

OUTPUT_PATH = os.path.join(RESULTS_DIR, f"final_report__{CANCER_TYPE}__{METRIC}.pdf")

# -------- Load metrics --------
num_genes = 0
num_outliers = 0
km_pval = "N/A"

if os.path.exists(TOP_GENES):
    df_genes = pd.read_csv(TOP_GENES, sep="\t")
    num_genes = len(df_genes)
else:
    df_genes = pd.DataFrame()

if os.path.exists(OUTLIERS):
    df_out = pd.read_csv(OUTLIERS, sep="\t")
    num_outliers = len(df_out)

if os.path.exists(KM_PLOT.replace(".png", ".txt")):
    with open(KM_PLOT.replace(".png", ".txt")) as f:
        km_pval = f.read().strip()

# -------- ReportLab styles --------
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="MainTitle", alignment=TA_CENTER, fontSize=18,
                          textColor="#660099", leading=22, spaceAfter=15))
styles.add(ParagraphStyle(name="CustomSubTitle", alignment=TA_CENTER, fontSize=12,
                          textColor=colors.HexColor("#555555"), spaceAfter=10))
styles.add(ParagraphStyle(name="Body", alignment=TA_LEFT, fontSize=10, leading=14))
styles.add(ParagraphStyle(name="Section", fontSize=14, leading=18,
                          textColor="#003366", spaceBefore=10, spaceAfter=6))

doc = SimpleDocTemplate(
    OUTPUT_PATH, pagesize=A4,
    rightMargin=40, leftMargin=40, topMargin=50, bottomMargin=40
)
story = []

# ==========================================================
# COVER PAGE
# ==========================================================
story.append(Spacer(1, 40))
story.append(Paragraph("Immuno-Cancer Project", styles["MainTitle"]))
story.append(Paragraph("Automated TCGA Immune Survival Report", styles["CustomSubTitle"]))
story.append(Spacer(1, 20))
story.append(Paragraph(f"<b>Cancer Type:</b> {CANCER_TYPE}", styles["Body"]))
story.append(Paragraph(f"<b>Metric:</b> {METRIC}", styles["Body"]))
story.append(Paragraph(f"<b>Date:</b> {datetime.today().strftime('%Y-%m-%d %H:%M')}", styles["Body"]))
story.append(Paragraph(f"<b>Genes analyzed:</b> {num_genes}", styles["Body"]))
story.append(Paragraph(f"<b>Outliers removed:</b> {num_outliers}", styles["Body"]))
story.append(Paragraph(f"<b>KM p-value:</b> {km_pval}", styles["Body"]))
story.append(Spacer(1, 60))
story.append(Paragraph("<b>Prepared by:</b>", styles["Body"]))
story.append(Paragraph("Oguzhan Begik, PhD", styles["Body"]))
story.append(Paragraph("Postdoctoral Scientist", styles["Body"]))
story.append(PageBreak())

# ==========================================================
# EXECUTIVE HIGHLIGHTS
# ==========================================================
story.append(Paragraph("Executive Highlights", styles["Section"]))
story.append(Paragraph(
    f"This report summarizes the outputs from the automated TCGA Immune Survival pipeline "
    f"for {CANCER_TYPE}, using immune dysfunction metric <b>{METRIC}</b>. "
    f"The workflow includes PCA-based outlier removal, PC1-driven survival stratification, "
    f"and semantic clustering of top immune-associated gene summaries.", styles["Body"]
))
story.append(Spacer(1, 10))

if not df_genes.empty:
    df_genes_sorted = df_genes.sort_values(by="PC1_Loading", key=abs, ascending=False)
    table_data = [["Gene", "PC1 Loading"]] + df_genes_sorted.head(30)[["Gene_Symbol", "PC1_Loading"]].values.tolist()
    table = Table(table_data, hAlign="LEFT", colWidths=[180, 150])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E6E6FA")),
        ("ALIGN", (1, 1), (-1, -1), "RIGHT"),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.whitesmoke, colors.lightgrey]),
    ]))
    story.append(table)
    story.append(Spacer(1, 15))

story.append(PageBreak())

# ==========================================================
# VISUALIZATION PANELS
# ==========================================================
plots = [
    ("PCA Outlier Detection (QC)", PCA_PLOT),
    ("PC1-Based Kaplan–Meier Survival", KM_PLOT),
    ("Top PC1 Gene Contributors", PC1_PLOT),
    ("Semantic Clustering of AI Gene Summaries", SEMANTIC_PLOT)
]

for title, img_path in plots:
    story.append(Paragraph(title, styles["Section"]))
    if os.path.exists(img_path):
        try:
            story.append(Image(img_path, width=400, height=300))
        except Exception as e:
            story.append(Paragraph(f"⚠️ Could not load {os.path.basename(img_path)}: {e}", styles["Body"]))
    else:
        story.append(Paragraph(f"⚠️ Image not found: {os.path.basename(img_path)}", styles["Body"]))
    story.append(PageBreak())

# ==========================================================
# GENE SUMMARIES
# ==========================================================
story.append(Paragraph("Gene Summary Insights", styles["Section"]))
if os.path.exists(SUMMARY_PATH):
    df_sum = pd.read_csv(SUMMARY_PATH, sep="\t")
    for _, row in df_sum.iterrows():
        story.append(Paragraph(f"<b>{row['Gene_Symbol']}</b>: {row['Summary']}", styles["Body"]))
        story.append(Spacer(1, 5))
else:
    story.append(Paragraph("Gene summary file not found.", styles["Body"]))

# ==========================================================
# FOOTER
# ==========================================================
story.append(Spacer(1, 40))
story.append(Paragraph(
    "<i>Generated automatically using Snakemake & Python — "
    "Reproducible TCGA Immune Survival Pipeline (v1.0)</i>", styles["CustomSubTitle"]
))

doc.build(story)
print(f"✅ Final integrated report saved to: {OUTPUT_PATH}")