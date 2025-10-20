# SurvAIval

A lightweight **AI-driven survival analysis pipeline** integrating TCGA expression, immune metrics, and clinical data.  
Built with **Python + Snakemake**, it demonstrates reproducible data science practices for translational oncology and biomarker discovery.

---

## 🚀 Key Features
- Automated **data integration** of TCGA expression and clinical metadata  
- **Principal Component Analysis (PCA)** and **Kaplan–Meier** stratification  
- **AI-based gene summaries** via Hugging Face LLMs  
- **Semantic clustering** of biological terms  
- End-to-end **PDF report generation**

---

## 🧩 Workflow Overview
```mermaid
graph TD
    A[Data Loading] --> B[Immune Gene Selection]
    B --> C[PCA & Survival Analysis]
    C --> D[AI Gene Summaries]
    D --> E[Semantic Clustering]
    E --> F[Final PDF Report]
