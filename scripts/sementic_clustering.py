import pandas as pd
import matplotlib.pyplot as plt
import umap
import numpy as np
import yaml, os
from sentence_transformers import SentenceTransformer
import hdbscan
from sklearn.metrics.pairwise import cosine_similarity

# === Load config ===
with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

CANCER_TYPE = config["cancer_type"]
METRIC = config["metric"]
RESULTS_DIR = config["results_dir"]

# === Paths ===
summ_path = os.path.join(RESULTS_DIR, f"PC1_gene_summaries_HF__{CANCER_TYPE}__{METRIC}.tsv")
out_plot = os.path.join(RESULTS_DIR, f"semantic_clusters_stable__{CANCER_TYPE}__{METRIC}.pdf")

# === Load summaries ===
df = pd.read_csv(summ_path, sep="\t")
texts = df["Summary"].fillna("").tolist()
genes = df["Gene_Symbol"].tolist()
print(f"🧠 Embedding {len(texts)} gene summaries using BioBERT...")

# === Step 1: BioBERT embeddings (biologically aware model) ===
model = SentenceTransformer("sentence-transformers/all-MiniLM-L12-v2")
embeddings = model.encode(texts, show_progress_bar=True, convert_to_numpy=True)

# === Step 2: Dimensionality reduction with UMAP ===
reducer = umap.UMAP(random_state=42, n_neighbors=8, min_dist=0.2, metric="cosine")
coords = reducer.fit_transform(embeddings)
df["UMAP1"], df["UMAP2"] = coords[:, 0], coords[:, 1]

# === Step 3: Semantic clustering using HDBSCAN ===
print("🔹 Running HDBSCAN for robust clustering...")
clusterer = hdbscan.HDBSCAN(min_cluster_size=3, metric='euclidean')
df["Cluster"] = clusterer.fit_predict(coords)

# If some points unclustered (-1), reassign to nearest cluster centroid
if (df["Cluster"] == -1).any():
    clustered = df[df["Cluster"] != -1]
    centroids = clustered.groupby("Cluster")[["UMAP1", "UMAP2"]].mean()
    for i, row in df[df["Cluster"] == -1].iterrows():
        # Convert centroids and row to numeric arrays for safe math
        centroid_array = centroids[["UMAP1", "UMAP2"]].to_numpy().astype(float)
        point = row[["UMAP1", "UMAP2"]].to_numpy().astype(float)
    
        # Compute Euclidean distances
        dists = np.linalg.norm(centroid_array - point, axis=1)
        nearest_cluster = centroids.index[np.argmin(dists)]
        df.loc[i, "Cluster"] = nearest_cluster
    print("⚠️ Unclustered points reassigned to nearest centroid.")

n_clusters = df["Cluster"].nunique()
print(f"✅ {n_clusters} stable clusters detected.")

# === Step 4: Evaluate intra-cluster coherence ===
for c in sorted(df["Cluster"].unique()):
    idx = df["Cluster"] == c
    mean_sim = cosine_similarity(embeddings[idx]).mean() if idx.sum() > 1 else 0
    print(f"Cluster {c}: {idx.sum()} genes, mean semantic similarity = {mean_sim:.3f}")

# === Step 5: Color palette ===
palette = ["#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69"]
colors = {c: palette[i % len(palette)] for i, c in enumerate(sorted(df["Cluster"].unique()))}

# === Step 6: Plot ===
plt.figure(figsize=(9,7))
for c in sorted(df["Cluster"].unique()):
    subset = df[df["Cluster"] == c]
    plt.scatter(subset["UMAP1"], subset["UMAP2"], s=140,
                color=colors[c], edgecolor="k", alpha=0.9, label=f"Cluster {c+1}")
    # Label cluster centroid
    cx, cy = subset["UMAP1"].mean(), subset["UMAP2"].mean()
    plt.text(cx, cy, f"Cluster {c+1}", fontsize=10, ha="center", va="center",
             bbox=dict(facecolor='white', alpha=0.7, boxstyle="round,pad=0.3"))

# Add gene names
for i, row in df.iterrows():
    plt.text(row["UMAP1"], row["UMAP2"], row["Gene_Symbol"], fontsize=8, ha="center", va="center")

plt.title(f"Semantic Clustering of AI-Generated Summaries ({CANCER_TYPE})", fontsize=14)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.legend(fontsize=9, loc="best")
plt.tight_layout()

# --- Export both formats ---
base_path = out_plot.replace(".pdf", "")
plt.savefig(base_path + ".pdf", dpi=300)
plt.savefig(base_path + ".png", dpi=300)
plt.close()

print(f"🧩 Stable semantic clustering plot saved → {base_path}.pdf and .png")