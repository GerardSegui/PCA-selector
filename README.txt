# PCA Selector 
---

## Table of Contents
1. Introduction
2. Requirements
3. How to Run
4. Overview
5. Label Extraction Algorithm
6. Normalizations
7. Metrics
8. Metabolite Scoring (NAD/FAD/PP)
9. PC1 Loadings Analysis
10. Output

---

## Introduction

PCA Selector is a software tool designed to find the optimal normalization for a given spectral dataset, with the goal of achieving maximum separation between biological conditions in PCA space.

Beyond simply computing a PCA, the software evaluates its quality in a quantitative and automated way. Instead of relying on visual inspection, multiple complementary metrics are combined to assess geometric separation, clustering consistency, and statistical significance. This approach allows for an objective and robust selection of the best data representation. Additionally, PC1 loading profiles are generated to identify which spectral regions contribute most strongly to the observed separation in PCA space.

Optionally, the user can provide a reference file containing the pure spectra of metabolites (NAD, FAD, and optionally PP). When provided, the software computes a fluorophore scoring for each cell and generates additional PCA plots colored by metabolite contribution, including a NAD/FAD ratio map.

The software is fully automated: the user only needs to select the dataset to be analyzed, and the system handles the rest of the process. At the end, a structured folder containing all results is generated automatically.

---

## Requirements

- The dataset must follow a samples x channels structure (rows = samples; columns = channels + label columns). Label columns are detected automatically—no manual specification is required.

- Label conditions must be separated using underscores (_) or spaces. Numerical values are automatically discarded, so if a number carries important information (e.g., week number), it should be attached to its corresponding label:

  - "WEEK 3" or "WEEK_3" -> NOT OK
  - "WEEK3" -> OK

- Labels must not be written as a single continuous string. Groups must be clearly separated; otherwise, the software will interpret them as a single condition and will not distinguish between groups.

- There are no restrictions on the number of channels or samples. The algorithm is designed to handle datasets of varying sizes.

- Supported file formats: .xlsx, .xls, .csv. An "All files" option may appear in the file dialog for navigation purposes only—files in unsupported formats will not be processed.

- If a metabolite reference file is provided, it must follow this structure:
  - First column: label (e.g., NAD, FAD, PP)
  - Remaining columns: channel intensities (must match the number of channels in the main dataset)

- The metabolite reference file must contain raw data. The ZIP file supports up to two raw data files: one containing NAD & FAD, and another including PP.

---

## How to Run

1. Run the main script (`pca_selector.py`).
2. A window appears — click **Select file** to choose the dataset (`.csv` or `.xlsx`).
3. A dialog will ask whether you want to load a **metabolite reference file** (NAD/FAD/PP). This step is optional but enables fluorophore scoring and colored PCA plots.
4. Enter a **name** for the output folder and files.
5. A window displays the **detected label groups**. Select which groups should define the PCA conditions, use **Preview** to verify, and click **OK** to confirm.
6. The analysis runs automatically.
7. A final window confirms the process has finished and shows the output folder path.

---

## Overview

Once launched, the software guides the user through a sequence of dialogs:

**File selection** → **Optional metabolite reference** → **Output name** → **Label group selection** → **Analysis** → **Results**

All intermediate plots and numeric results are saved automatically to a structured output folder. No manual intervention is required during the analysis itself.

---

## Label Extraction Algorithm

The software includes a custom algorithm to extract biological conditions from raw label columns.

**Step 1 — Unification:** if multiple label columns are present, they are merged into a single unified label per sample.

**Step 2 — Tokenization:** each label is split by `_` into tokens, after removing common file extensions (`.tif`, `.png`, `.jpg`).

**Step 3 — Informativeness filtering:** tokens that appear in all samples (constants) or in only one sample (noise) are discarded. Purely numeric tokens are also excluded. Only tokens that show variability across samples are kept.

**Step 4 — Grouping:** tokens that never co-occur in the same label are grouped together, as they likely represent mutually exclusive conditions (e.g., different treatments or time points).

**Step 5 — User selection:** the detected groups are shown in a GUI window. The user selects which groups should define the PCA conditions and can preview the resulting labels before confirming.

This algorithm has been tested on diverse datasets and handles multi-condition, multi-factor label structures robustly.

---

## Normalizations

Four normalization strategies are evaluated in parallel. Each produces a different representation of the data, and the best one is selected based on the metrics described below.

| Name | Description |
|---|---|
| **MyFacts** | Log2 transformation followed by median-of-ratios normalization. Corrects for differences in total intensity across samples while preserving relative channel structure. Falls back to intensity normalization when applied to single reference vectors. |
| **Intensity** | Each spectrum is divided by its total intensity (sum = 1). Removes absolute intensity differences, preserving only spectral shape. |
| **Log1p** | Applies `log(1 + x)` to each value. Compresses dynamic range and reduces the influence of high-intensity channels. |
| **L2** | Each spectrum is normalized to unit Euclidean norm. Suitable when the direction of the spectral vector is more informative than its magnitude. |

---

## Metrics

Five metrics are computed for each normalization. All are normalized to [0, 1] and combined into a single weighted score. Only normalizations that pass the PERMANOVA significance threshold (p < 0.05) are considered valid candidates for the best PCA.

### Separation (weight: 0.30)
A custom score defined as the ratio between inter-cluster distance and intra-cluster distance. Inter-cluster distance is the mean pairwise distance between group centroids; intra-cluster distance is the mean distance of each point to its group centroid. Higher values indicate better-separated, more compact groups. Provides a global geometric view of the cluster structure.

### Silhouette (weight: 0.20)
Evaluates how well each individual point fits within its assigned cluster compared to the nearest alternative cluster. Values range from -1 to 1: values near 1 indicate well-separated clusters, values near 0 indicate overlap, and negative values suggest possible misassignment. Provides local, per-sample information about clustering quality.

### Davies-Bouldin — DB (weight: 0.05)
Measures the average similarity between each cluster and its most similar neighbour, considering both intra-cluster spread and inter-cluster distance. Lower values indicate better clustering. In this software the metric is inverted so that higher values correspond to better results, keeping it consistent with the other metrics.

### Calinski-Harabasz — CH (weight: 0.15)
Defined as the ratio between between-cluster variance and within-cluster variance. Higher values indicate more compact, well-separated clusters. More statistically grounded than the separation score, as it accounts for the full distribution of points.

### PERMANOVA (weight: 0.30)
A permutational multivariate analysis of variance test that evaluates whether the observed separation between conditions is statistically significant or could have arisen by chance. A p-value below 0.05 is required for a normalization to be considered valid. When included in the score, lower p-values contribute proportionally more (via `-log10(p)` transformation).

---

## Metabolite Scoring (NAD/FAD/PP)

When a reference file is provided, the software computes a **fluorophore score** for each cell spectrum using a weighted projection method.

### How it works

For each reference spectrum (NAD, FAD, and optionally PP), a weight vector is computed that emphasizes channels where that fluorophore is dominant and penalizes channels where it overlaps with the others. The weight formula for NAD (with PP present) is:

```
w_NAD = (NAD / total) × ((NAD − (FAD + PP) / 2) / total)
```

The score for each cell is then the dot product of its normalized spectrum with the weight vector. This acts as a fixed spectral filter derived entirely from the reference spectra — it does not change with the dataset being analyzed.

### Colored PCA plots

For each normalization, four additional PCA plots are generated:

| Plot | Description |
|---|---|
| `*_NAD.png` | Each point colored by NAD score (viridis colormap) |
| `*_FAD.png` | Each point colored by FAD score (viridis colormap) |
| `*_PP.png` | Each point colored by PP score (viridis colormap) — only if PP is present |
| `*_NADperFAD.png` | Each point colored by NAD/FAD ratio (RdBu_r: blue = FAD dominant, red = NAD dominant) |

In all colored plots, the marker shape encodes the condition group, allowing simultaneous visualization of metabolite contribution and biological condition per cell.

### Biological validation

The Spearman correlation between the NAD/FAD channel ratio and each principal component (PC1 and PC2) is computed. If the absolute correlation exceeds 0.5 and the p-value is below 0.05, the normalization is considered **biologically valid** — meaning the PCA captures the metabolic state of the cells. Normalizations that are both statistically significant (PERMANOVA) and biologically valid are preferred over those that are only statistically significant.

---

## PC1 Loadings Analysis

In addition to the PCA projections, the software computes and stores the loadings associated with the first principal component (PC1) for each normalization.

The loading values quantify the contribution of each spectral channel to PC1, allowing identification of the regions of the spectrum that drive the separation between biological conditions.

### Generated plots

For every normalization, a loading plot is generated showing:

- Spectral channel (x-axis)
- PC1 loading value (y-axis)

Positive and negative loading regions indicate channels contributing in opposite directions along the principal component axis.

These plots provide an interpretable link between the statistical PCA separation and the underlying spectral features responsible for that separation.

### Biological interpretation

Channels with high absolute loading values are the most influential in defining PC1. Comparing loading profiles across normalizations can help identify robust spectral biomarkers and determine whether the PCA separation is driven by biologically meaningful regions of the spectrum.

---

## Output

For each execution, a folder named `best_pca_<name>` is generated with the following structure:

```
best_pca_<name>/
│
├── best_pca/
│   └── <best_normalization>.png        # PCA of the best normalization
│
├── all_pca/
│   └── <normalization>.png             # PCA for each normalization
│
├── mean_spectra/
│   └── <normalization>_spectra.png     # Mean ± SD spectra per condition
│
├── nadfad_spectra/
│   └── <normalization>_nadfad.png      # NAD/FAD/PP reference spectra (if provided)
│
├── pca_colored/
│   ├── <normalization>_NAD.png         # PCA colored by NAD score
│   ├── <normalization>_FAD.png         # PCA colored by FAD score
│   ├── <normalization>_PP.png          # PCA colored by PP score (if present)
│   └── <normalization>_NADperFAD.png   # PCA colored by NAD/FAD ratio
|
├── pc1_loadings/
│   └── <normalization>_pc1_loadings.png   # PC1 loading profile for each normalization
│
└── <name>.txt                      # Full numeric results and scoring breakdown
```

The txt file contains, for each normalization: separation score, silhouette, Davies-Bouldin, Calinski-Harabasz, PERMANOVA p-value, Spearman correlation, biological validation result, normalized metric values, final weighted score, and an explanation of why the best normalization was selected.