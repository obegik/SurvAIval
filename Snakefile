############################################################
# Snakefile — TCGA Immune Survival Pipeline
# Author: Oguzhan Begik, PhD
# ----------------------------------------------------------
# Full pipeline integrating all modular scripts:
#   data_loading.py
#   subset_expression.py
#   select_top_genes.py
#   qc_outlier_removal.py
#   pc1_survival.py
#   sementic_clustering.py
#   generate_report.py
#
# Output: results/final_report__{CANCER_TYPE}__{METRIC}.pdf
############################################################

configfile: "config/config.yaml"

CANCER_TYPE = config["cancer_type"]
METRIC = config["metric"]
RESULTS_DIR = config["results_dir"]

# ----------------------------------------------------------
rule all:
    input:
        f"{RESULTS_DIR}/final_report__{CANCER_TYPE}__{METRIC}.pdf"

# ----------------------------------------------------------
rule data_loading:
    output:
        expr=f"{RESULTS_DIR}/expr__{CANCER_TYPE}.parquet",
        clin=f"{RESULTS_DIR}/clinical__{CANCER_TYPE}.tsv"
    script:
        "scripts/data_loading.py"

# ----------------------------------------------------------
rule select_top_genes:
    input:
        scores="data/all_scores.tsv"
    output:
        top_genes=f"{RESULTS_DIR}/top_immune_genes__{CANCER_TYPE}__{METRIC}.tsv"
    script:
        "scripts/select_top_genes.py"

# ----------------------------------------------------------
rule subset_expression:
    input:
        expr=f"{RESULTS_DIR}/expr__{CANCER_TYPE}.parquet",
        top_genes=f"{RESULTS_DIR}/top_immune_genes__{CANCER_TYPE}__{METRIC}.tsv"
    output:
        matrix=f"{RESULTS_DIR}/expr_topimmune__{CANCER_TYPE}__{METRIC}.parquet",
        annotated=f"{RESULTS_DIR}/expr_topimmune_annotated__{CANCER_TYPE}__{METRIC}.parquet"
    script:
        "scripts/subset_expression.py"

# ----------------------------------------------------------
rule remove_outliers:
    input:
        matrix=f"{RESULTS_DIR}/expr_topimmune__{CANCER_TYPE}__{METRIC}.parquet",
        clin=f"{RESULTS_DIR}/clinical__{CANCER_TYPE}.tsv"
    output:
        expr_clean=f"{RESULTS_DIR}/expr_topimmune_cleaned__{CANCER_TYPE}__{METRIC}.parquet",
        clin_clean=f"{RESULTS_DIR}/clinical_cleaned__{CANCER_TYPE}__{METRIC}.tsv",
        outliers=f"{RESULTS_DIR}/outliers__{CANCER_TYPE}__{METRIC}.tsv",
        pca_plot=f"{RESULTS_DIR}/qc_plots/pca_outlier_detection__{CANCER_TYPE}__{METRIC}.png"
    script:
        "scripts/qc_outlier_removal.py"

# ----------------------------------------------------------
rule pca_km:
    input:
        expr=f"{RESULTS_DIR}/expr_topimmune_cleaned__{CANCER_TYPE}__{METRIC}.parquet",
        clin=f"{RESULTS_DIR}/clinical_cleaned__{CANCER_TYPE}__{METRIC}.tsv"
    output:
        km_plot=f"{RESULTS_DIR}/km_PC1_split__{CANCER_TYPE}__{METRIC}.png",
        loadings=f"{RESULTS_DIR}/PC1_gene_loadings__{CANCER_TYPE}__{METRIC}.tsv",
        annotated=f"{RESULTS_DIR}/PC1_gene_loadings_annotated__{CANCER_TYPE}__{METRIC}.tsv",
        contributors=f"{RESULTS_DIR}/top_PC1_contributor_genes__{CANCER_TYPE}__{METRIC}.tsv",
        contrib_plot=f"{RESULTS_DIR}/PC1_top_gene_contributors_annotated__{CANCER_TYPE}__{METRIC}.png"
    script:
        "scripts/pc1_survival.py"

# ----------------------------------------------------------

rule gene_summary:
    input:
        contributors=f"{RESULTS_DIR}/top_PC1_contributor_genes__{CANCER_TYPE}__{METRIC}.tsv"
    output:
        summary=f"{RESULTS_DIR}/PC1_gene_summaries_HF__{CANCER_TYPE}__{METRIC}.tsv"
    script:
        "scripts/gene_summary.py"
        
# ----------------------------------------------------------
rule sementic_clustering:
    input:
        summaries=f"{RESULTS_DIR}/PC1_gene_summaries_HF__{CANCER_TYPE}__{METRIC}.tsv"
    output:
        sem_plot=f"{RESULTS_DIR}/semantic_clusters_stable__{CANCER_TYPE}__{METRIC}.png"
    script:
        "scripts/sementic_clustering.py"

# ----------------------------------------------------------
rule generate_report:
    input:
        km_plot=f"{RESULTS_DIR}/km_PC1_split__{CANCER_TYPE}__{METRIC}.png",
        pc1_plot=f"{RESULTS_DIR}/PC1_top_gene_contributors_annotated__{CANCER_TYPE}__{METRIC}.png",
        qc_plot=f"{RESULTS_DIR}/qc_plots/pca_outlier_detection__{CANCER_TYPE}__{METRIC}.png",
        sem_plot=f"{RESULTS_DIR}/semantic_clusters_stable__{CANCER_TYPE}__{METRIC}.png",
        top_genes=f"{RESULTS_DIR}/top_PC1_contributor_genes__{CANCER_TYPE}__{METRIC}.tsv",
        outliers=f"{RESULTS_DIR}/outliers__{CANCER_TYPE}__{METRIC}.tsv",
        summary=f"{RESULTS_DIR}/PC1_gene_summaries_HF__{CANCER_TYPE}__{METRIC}.tsv"
    output:
        report=f"{RESULTS_DIR}/final_report__{CANCER_TYPE}__{METRIC}.pdf"
    script:
        "scripts/generate_report.py"