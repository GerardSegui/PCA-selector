#!/usr/bin/env python
# coding: utf-8

# In[2]:


#------ABSTRACT------#
#This is a code that opens a UI when compiled. You will have to choose a file which will be normalized by different normalizations and a PCA
#will be done for each normalization and, after that, by doing some calculations, the best PCA will be chosen and returned to the user.

#------IMPORTS------#
import tkinter as tk
from tkinter import *
from tkinter import filedialog, Toplevel
from tkinter import messagebox
import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms  
import sys
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler, QuantileTransformer, PowerTransformer, normalize
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from scipy.spatial.distance import pdist, cdist
from scipy.stats import spearmanr
from skbio.stats.distance import DistanceMatrix, permanova
from collections import defaultdict
from collections import Counter

NADH_CHANNELS = list(range(0, 8))   # columnes 1-8 (índex 0-7)
FAD_CHANNELS  = list(range(13, 20)) # columnes 14-20 (índex 13-19)
WAVELENGTHS = list(np.arange(415, 696, 9))  # 32 canals: 415, 424, 433, ... 694 nm
SPEARMAN_THRESHOLD = 0.5

#------LOAD DATA------#

def load_csv(filepath):
    return pd.read_csv(filepath)

def load_excel(filepath):
    return pd.read_excel(filepath, sheet_name=None)

def load_data(filepath):
    _, ext = os.path.splitext(filepath)

    if ext == ".csv":
        return load_csv(filepath)
    elif ext in [".xlsx", ".xls"]:
        return load_excel(filepath)
    else:
        raise ValueError("Not a correct file")


#------LABEL EXTRACTION PIPELINE------#

def unify_labels(df_labels):
    if isinstance(df_labels, pd.DataFrame):
        if df_labels.shape[1] > 1:
            labels = df_labels.astype(str).agg("_".join, axis=1)
        else:
            labels = df_labels.iloc[:,0].astype(str)
    else:
        labels = pd.Series(df_labels).astype(str)
    return labels


def clean_label(label):
    for ext in [".tif", ".png", ".jpg"]:
        label = label.replace(ext, "")
    return label


def tokenize_labels(labels):
    return [re.split(r'[_\s]+', clean_label(lbl)) for lbl in labels]


def select_informative_tokens(tokenized, max_unique=10):
    n_samples = len(tokenized)
    token_presence = defaultdict(set)

    for i, tokens in enumerate(tokenized):
        for t in tokens:
            token_presence[t].add(i)

    informative = []
    for token, sample_indices in token_presence.items():
        count = len(sample_indices)
        if token.isdigit():
            continue
        if count == 1 or count == n_samples:
            continue
        informative.append(token)

    return informative


def build_token_groups(tokenized, informative_tokens):
    
    n = len(informative_tokens)
    tok_set = set(informative_tokens)
    idx = {t: i for i, t in enumerate(informative_tokens)}

    cooccur = np.zeros((n, n), dtype=int)

    for tokens in tokenized:
        present = [t for t in tokens if t in tok_set]
        for a in present:
            for b in present:
                if a != b:
                    cooccur[idx[a], idx[b]] += 1

    visited = [False] * n
    groups = []

    for i in range(n):
        if visited[i]:
            continue
        group = [informative_tokens[i]]
        visited[i] = True
        for j in range(i + 1, n):
            if not visited[j] and cooccur[i, j] == 0:
                if all(cooccur[idx[g], j] == 0 for g in group):
                    group.append(informative_tokens[j])
                    visited[j] = True
        groups.append(group)

    return groups


def build_labels_from_groups(tokenized, groups):
    labels = []

    for t in tokenized:
        token_set = set(t)
        parts = []

        for group in groups:
            match = next((tok for tok in group if tok in token_set), None)
            if match is not None:
                parts.append(match)
            else:
                parts.append("ctrl")

        labels.append("_".join(parts) if parts else "UNKNOWN")

    return labels


#------WINDOW FOR CONDITIONS------#

def group_selection_gui(groups, tokenized):

    result = {"groups": groups}

    top = Toplevel()
    top.title("Select Label Groups")
    top.geometry("500x600")

    Label(top, text="Detected token groups:", font=("Arial", 12, "bold")).pack(pady=10)

    vars = []
    for i, g in enumerate(groups):
        var = IntVar(value=1)
        chk = Checkbutton(top, text=f"[{i}] -> {g}", variable=var, anchor="w", justify="left")
        chk.pack(fill="x", padx=10)
        vars.append(var)

    Label(top, text="\nExamples of reconstructed labels:").pack(pady=10)

    preview_text = StringVar()
    Label(top, textvariable=preview_text, wraplength=450, justify="left").pack()

    def update_preview():
        selected = [groups[i] for i in range(len(groups)) if vars[i].get() == 1]
        labels_preview = build_labels_from_groups(tokenized, selected)
        
        seen = {}
        for label in labels_preview:
            if label not in seen:
                seen[label] = label
    
        preview_text.set("\n".join(seen.keys()))

    Button(top, text="Preview", command=update_preview).pack(pady=5)

    def confirm():
        selected = [groups[i] for i in range(len(groups)) if vars[i].get() == 1]
        if not selected:
            selected = groups
        result["groups"] = selected
        top.destroy()

    Button(top, text="OK", command=confirm).pack(pady=10)

    top.grab_set()
    top.wait_window()

    return result["groups"]


#------PIPELINE EXTRACT CONDITIONS------#

def extract_conditions_pipeline(df_labels):

    labels = unify_labels(df_labels)
    labels = labels.apply(clean_label)
    tokenized = tokenize_labels(labels)

    informative_tokens = select_informative_tokens(tokenized)
    print("\n Informative tokens:", informative_tokens)

    groups = build_token_groups(tokenized, informative_tokens)
    selected_groups = group_selection_gui(groups, tokenized)

    final_conditions = np.array(
        build_labels_from_groups(tokenized, selected_groups)
    )

    return final_conditions

#------NORMALIZATIONS------#

def normalize_myfacts(data):
    mat = data.T
    mat = np.log2(mat + 1)
    mns = np.exp(np.mean(np.log(mat + 1), axis=1))
    matm = mat / mns[:, None]
    facts = np.median(matm, axis=0)
    return (mat / facts).T

def normalize_intensity(data):
    return data / data.sum(axis=1, keepdims=True)

def normalize_log1p(data):
    return np.log1p(data)

def normalize_l2(data):
    return normalize(data, norm='l2')

#------APPLY NORMALIZATION TO A SINGLE VECTOR------#

def apply_normalization_to_vector(vec, norm_name, data_ref):
    """
    Applies the normalization fitted on data_ref to a single 1D vector (e.g. NAD or FAD spectrum).
    vec      : 1D array of shape (n_channels,)
    norm_name: name of the normalization
    data_ref : the original cell data matrix (n_cells x n_channels), used to fit scalers
    Returns a normalized 1D array.
    """
    vec = vec.astype(float).reshape(1, -1)  # shape (1, n_channels)

    if norm_name == "MyFacts":
        # MyFacts requires a full matrix to be meaningful — for a single vector,
        # fall back to intensity normalization (sum = 1)
        total = vec.sum()
        return (vec / (total + 1e-12)).flatten()

    elif norm_name == "Intensity":
        total = vec.sum()
        return (vec / (total + 1e-12)).flatten()

    elif norm_name == "Log1p":
        return np.log1p(vec).flatten()

    elif norm_name == "L2":
        return normalize(vec, norm='l2').flatten()

    else:
        return vec.flatten()


#------METRICS------#

def compute_centroid(points):
    return points.mean(axis=0)


def separation_score(X, conditions):
    groups = np.unique(conditions)
    centroids, intra = [], []

    for g in groups:
        idx = np.where(conditions == g)[0]
        subset = X[idx, :2]
        c = compute_centroid(subset)
        centroids.append(c)
        intra.extend(np.linalg.norm(subset - c, axis=1))

    return np.mean(pdist(np.array(centroids))) / np.mean(intra)


def permanova_test(X, conditions):
    dm = DistanceMatrix(cdist(X, X))
    return permanova(dm, conditions)['p-value']


#------BIOLOGICAL VALIDATION------#

def biological_validation(X, data_norm, conditions):

    try:
        nadh = data_norm[:, NADH_CHANNELS].sum(axis=1)
        fad  = data_norm[:, FAD_CHANNELS].sum(axis=1)
        ratio = nadh / (fad + 1e-6)

        rho1, p1 = spearmanr(ratio, X[:, 0])
        rho2, p2 = spearmanr(ratio, X[:, 1])

        if abs(rho1) >= abs(rho2):
            rho, p_val, pc = rho1, p1, "PC1"
        else:
            rho, p_val, pc = rho2, p2, "PC2"

        is_valid = (abs(rho) >= SPEARMAN_THRESHOLD) and (p_val < 0.05)

        return is_valid, round(rho, 4), pc

    except Exception as e:
        print(f"  Biological validation failed: {e}")
        return False, 0.0, "N/A"


#------PCA + ANALYSIS------#

def run_analysis(data, conditions):
    pca = PCA(n_components=2)
    X = pca.fit_transform(data)

    labels_num = pd.factorize(conditions)[0]

    return (
        X, pca,
        separation_score(X, conditions),
        silhouette_score(X, labels_num),
        davies_bouldin_score(X, labels_num),
        calinski_harabasz_score(X, labels_num),
        permanova_test(X, conditions)
    )

#------ELIPSES------#

def draw_ellipse(ax, points, n_std=2.0):

    if len(points) < 2:
        return

    cov = np.cov(points, rowvar=False)
    mean = np.mean(points, axis=0)

    eigvals, eigvecs = np.linalg.eigh(cov)

    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))

    width, height = 2 * n_std * np.sqrt(eigvals)

    ellipse = Ellipse(
        xy=mean,
        width=width,
        height=height,
        angle=angle,
        edgecolor='black',
        facecolor='none',
        linestyle='--',
        linewidth=1.5
    )

    ax.add_patch(ellipse)

#------SAVE PCA------#

def save_pca_plot(folder, filename, X, conditions, title, p_value=None, pca=None):

    os.makedirs(folder, exist_ok=True)
    filepath = os.path.join(folder, filename)

    fig, ax = plt.subplots()

    unique_conditions = np.unique(conditions)
    
    # Colors més diferenciats: combina tab10 i Set1 per evitar similituds
    base_colors = [
        "#1f77b4", "#d62728", "#2ca02c", "#9467bd",
        "#e77f00", "#17becf", "#8c564b", "#e377c2",
        "#000080", "#006400", "#8B0000", "#4B0082"
    ]
    colors_list = (base_colors * ((len(unique_conditions) // len(base_colors)) + 1))[:len(unique_conditions)]

    for i, cond in enumerate(unique_conditions):
        idx = conditions == cond
        points = X[idx, :2]
        color = colors_list[i]

        ax.scatter(points[:, 0], points[:, 1], label=cond, color=color, alpha=0.7)

        centroid = np.mean(points, axis=0)
        ax.scatter(centroid[0], centroid[1], color=color, marker='x', s=120, linewidths=3)

        if len(points) > 2:
            cov = np.cov(points, rowvar=False)
            eigvals, eigvecs = np.linalg.eigh(cov)
            order = eigvals.argsort()[::-1]
            eigvals, eigvecs = eigvals[order], eigvecs[:, order]
            angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
            width, height = 2 * 2.0 * np.sqrt(eigvals)

            ellipse = Ellipse(
                xy=centroid, width=width, height=height, angle=angle,
                facecolor=color, edgecolor=color,
                alpha=0.08,          # ← abans 0.2, ara molt més transparent
                linewidth=1.5
            )
            ax.add_patch(ellipse)

    if pca is not None:
        var1 = pca.explained_variance_ratio_[0] * 100
        var2 = pca.explained_variance_ratio_[1] * 100
        ax.set_xlabel(f"PC1 ({var1:.1f}%)")
        ax.set_ylabel(f"PC2 ({var2:.1f}%)")
    else:
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    ax.set_title(title)
    ax.legend()

    if p_value is not None:
        ax.text(0.99, 0.01, f"p = {p_value:.2e}", transform=ax.transAxes,
                fontsize=9, ha='right', va='bottom')

    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300)
    plt.close()


def save_pca_colored_plot(folder, filename, X, scores, score_label, title, pca=None, conditions=None):
    os.makedirs(folder, exist_ok=True)
    filepath = os.path.join(folder, filename)

    fig, ax = plt.subplots()

    markers = ['o', '^', 's', 'D', 'P', '*', 'X', 'v', '<', '>']
    base_colors = [
        "#1f77b4", "#d62728", "#2ca02c", "#9467bd",
        "#e77f00", "#17becf", "#8c564b", "#e377c2",
        "#000080", "#006400", "#8B0000", "#4B0082"
    ]

    sc = None

    if conditions is not None:
        unique_conditions = np.unique(conditions)
        colors_list = (base_colors * ((len(unique_conditions) // len(base_colors)) + 1))[:len(unique_conditions)]

        for i, cond in enumerate(unique_conditions):
            idx = conditions == cond
            sc = ax.scatter(
                X[idx, 0], X[idx, 1],
                c=scores[idx], cmap='viridis',
                alpha=0.7, s=30,
                marker=markers[i % len(markers)],
                label=cond,
                vmin=scores.min(), vmax=scores.max()  # escala comuna per tots els grups
            )

        ax.legend(fontsize=8, loc='upper right')

    else:
        sc = ax.scatter(X[:, 0], X[:, 1], c=scores, cmap='viridis', alpha=0.7, s=20)

    if sc is not None:
        plt.colorbar(sc, ax=ax, label=score_label)

    if pca is not None:
        var1 = pca.explained_variance_ratio_[0] * 100
        var2 = pca.explained_variance_ratio_[1] * 100
        ax.set_xlabel(f"PC1 ({var1:.1f}%)")
        ax.set_ylabel(f"PC2 ({var2:.1f}%)")
    else:
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    ax.set_title(title)
    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300)
    plt.close()


def save_pca_nadfad_ratio_plot(folder, filename, X, score_nad, score_fad, title, pca=None, conditions=None):
    os.makedirs(folder, exist_ok=True)
    filepath = os.path.join(folder, filename)

    ratio = score_nad / (score_fad + 1e-12)
    p5, p95 = np.percentile(ratio, 5), np.percentile(ratio, 95)
    ratio_clipped = np.clip(ratio, p5, p95)

    fig, ax = plt.subplots()

    markers = ['o', '^', 's', 'D', 'P', '*', 'X', 'v', '<', '>']

    sc = None

    if conditions is not None:
        unique_conditions = np.unique(conditions)

        for i, cond in enumerate(unique_conditions):
            idx = conditions == cond
            sc = ax.scatter(
                X[idx, 0], X[idx, 1],
                c=ratio_clipped[idx], cmap='RdBu_r',
                alpha=0.7, s=30,
                marker=markers[i % len(markers)],
                label=cond,
                vmin=ratio_clipped.min(), vmax=ratio_clipped.max()
            )

        ax.legend(fontsize=8, loc='upper right')

    else:
        sc = ax.scatter(X[:, 0], X[:, 1], c=ratio_clipped, cmap='RdBu_r', alpha=0.7, s=20)

    if sc is not None:
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("NAD / FAD score ratio\n(blau=FAD, vermell=NAD)")

    if pca is not None:
        var1 = pca.explained_variance_ratio_[0] * 100
        var2 = pca.explained_variance_ratio_[1] * 100
        ax.set_xlabel(f"PC1 ({var1:.1f}%)")
        ax.set_ylabel(f"PC2 ({var2:.1f}%)")
    else:
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    ax.set_title(title)
    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300)
    plt.close()

#------SAVE MEAN SPECTRA PLOT (cells)------#

def save_mean_spectra_plot(folder, filename, data_norm, conditions, title):
    """
    For each condition: plots the mean spectrum as a solid line and shades
    the region between (mean - std) and (mean + std). Each condition gets
    a distinct colour. All individual cell spectra are shown as faint lines.
    Bottom X axis: channel numbers. Top X axis: wavelengths in nm (clipped
    to the actual number of channels in the data).
    """
    os.makedirs(folder, exist_ok=True)
    filepath = os.path.join(folder, filename)

    n_channels = data_norm.shape[1]
    wavelengths_used = np.array(WAVELENGTHS[:n_channels])
    channels = np.arange(1, n_channels + 1)

    unique_conditions = np.unique(conditions)
    cmap = plt.colormaps['tab10'].resampled(len(unique_conditions))

    fig, ax = plt.subplots(figsize=(10, 5))

    for i, cond in enumerate(unique_conditions):
        idx = np.where(conditions == cond)[0]
        subset = data_norm[idx]           # shape: (n_cells, n_channels)
        color = cmap(i)

        # Individual cell spectra (faint)
        for cell_spectrum in subset:
            ax.plot(channels, cell_spectrum, color=color, alpha=0.08, linewidth=0.5)

        mean_spec = subset.mean(axis=0)
        std_spec  = subset.std(axis=0)

        ax.fill_between(
            channels,
            mean_spec - std_spec,
            mean_spec + std_spec,
            color=color, alpha=0.25, label=f"{cond} ± 1 SD"
        )
        ax.plot(channels, mean_spec, color=color, linewidth=2, label=f"{cond} mean")

    ax.set_xlabel("Channel")
    ax.set_ylabel("Intensity (normalized)")
    ax.set_title(title)
    ax.set_xlim(channels[0], channels[-1])
    ax.set_xticks(channels)
    ax.tick_params(axis='x', labelsize=7)

    # ---- Top X axis: wavelengths ----
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    
    # Mostrar tots els canals
    tick_positions = channels
    tick_labels = [f"{int(w)}" for w in wavelengths_used]
    
    ax_top.set_xticks(tick_positions)
    ax_top.set_xticklabels(tick_labels, fontsize=7, rotation=90)
    
    ax_top.set_xlabel("Wavelength (nm)")

    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300)
    plt.close()


#------SAVE NAD/FAD REFERENCE SPECTRA PLOT------#

def save_nadfad_spectra_plot(folder, filename, nad_norm, fad_norm, title, pp_norm=None):
    """
    Plots the normalized NAD, FAD and optionally PP reference spectra on the same axes.
    nad_norm, fad_norm: 1D arrays of shape (n_channels,)
    pp_norm: optional 1D array of shape (n_channels,)
    Bottom X axis: channel numbers. Top X axis: wavelengths in nm (clipped
    to the actual number of channels in the data).
    """
    os.makedirs(folder, exist_ok=True)
    filepath = os.path.join(folder, filename)

    n_channels = len(nad_norm)
    wavelengths_used = np.array(WAVELENGTHS[:n_channels])
    channels = np.arange(1, n_channels + 1)

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(channels, nad_norm, color='royalblue',  linewidth=2, label='NAD (normalized)')
    ax.plot(channels, fad_norm, color='darkorange', linewidth=2, label='FAD (normalized)')
    if pp_norm is not None:
        ax.plot(channels, pp_norm, color='firebrick', linewidth=2, label='PP (normalized)')

    ax.set_xlabel("Channel")
    ax.set_ylabel("Intensity (normalized)")
    ax.set_title(title)
    ax.set_xlim(channels[0], channels[-1])
    ax.set_xticks(channels)
    ax.tick_params(axis='x', labelsize=7)

    # ---- Top X axis: wavelengths ----
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    
    # Mostrar tots els canals
    tick_positions = channels
    tick_labels = [f"{int(w)}" for w in wavelengths_used]
    
    ax_top.set_xticks(tick_positions)
    ax_top.set_xticklabels(tick_labels, fontsize=7, rotation=90)
    
    ax_top.set_xlabel("Wavelength (nm)")

    ax.legend(fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300)
    plt.close()


#------SAVE PC1 LOADINGS PLOT------#

def save_pc1_loadings_plot(folder, filename, pca, title):
    """
    Saves a plot of PC1 loadings across channels/wavelengths.
    """

    os.makedirs(folder, exist_ok=True)
    filepath = os.path.join(folder, filename)

    # Loadings de la PC1
    loadings = pca.components_[0]

    n_channels = len(loadings)

    wavelengths_used = np.array(WAVELENGTHS[:n_channels])
    channels = np.arange(1, n_channels + 1)

    fig, ax = plt.subplots(figsize=(10, 5))

    # línia principal
    ax.plot(channels, loadings, linewidth=2)

    # línia y=0
    ax.axhline(0, linestyle='--', linewidth=1)

    # labels
    ax.set_xlabel("Channel")
    ax.set_ylabel("PC1 loading")
    ax.set_title(title)

    ax.set_xlim(channels[0], channels[-1])

    # Mostrar tots els canals
    ax.set_xticks(channels)
    ax.tick_params(axis='x', labelsize=7)

    # ---- Top X axis: wavelengths ----
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())

    ax_top.set_xticks(channels)
    ax_top.set_xticklabels(
        [str(int(w)) for w in wavelengths_used],
        fontsize=7,
        rotation=90
    )

    ax_top.set_xlabel("Wavelength (nm)")

    ax.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig(filepath, dpi=300)
    plt.close()

#------LOAD NAD/FAD REFERENCE FILE------#

def load_nadfad_reference(filepath):
    """
    Loads an Excel/CSV where:
      - First column  : label  (e.g. 'NAD', 'FAD')
      - Remaining cols: channel intensities

    Returns two 1D numpy arrays: nad_raw, fad_raw
    """
    _, ext = os.path.splitext(filepath)
    if ext in [".xlsx", ".xls"]:
        df = pd.read_excel(filepath, header=0)
    else:
        df = pd.read_csv(filepath, header=0)

    # First column is the label
    label_col = df.columns[0]
    df[label_col] = df[label_col].astype(str).str.strip().str.upper()

    nad_rows = df[df[label_col].str.contains("NAD")]
    fad_rows = df[df[label_col].str.contains("FAD")]
    pp_rows = df[df[label_col].str.contains("PP")]

    if nad_rows.empty or fad_rows.empty:
        raise ValueError("Could not find NAD or FAD rows in the reference file. "
                         "Make sure the first column contains 'NAD' and 'FAD' labels.")

    nad_raw = nad_rows.iloc[0, 1:].values.astype(float)
    fad_raw = fad_rows.iloc[0, 1:].values.astype(float)
    pp_raw = pp_rows.iloc[0, 1:].values.astype(float) if not pp_rows.empty else None

    return nad_raw, fad_raw, pp_raw

def compute_weights(nad, fad, pp=None):
    """
    Calcula els vectors de pesos per NAD, FAD i opcionalment PP.
    Els espectres d'entrada han d'estar ja normalitzats (en el mateix espai que les dades).
    """
    if pp is None:
        total = nad + fad
        w_nad = (nad / (total + 1e-12)) * (np.abs(nad - fad) / (total + 1e-12))
        w_fad = (fad / (total + 1e-12)) * (np.abs(fad - nad) / (total + 1e-12))
        return w_nad, w_fad, None
    else:
        total = nad + fad + pp
        w_nad = (nad / (total + 1e-12)) * ((nad - (fad + pp) / 2) / (total + 1e-12))
        w_fad = (fad / (total + 1e-12)) * ((fad - (nad + pp) / 2) / (total + 1e-12))
        w_pp  = (pp  / (total + 1e-12)) * ((pp  - (nad + fad) / 2) / (total + 1e-12))
        return w_nad, w_fad, w_pp


#------ASK SAVED FILE NAME------#

def ask_filename_gui(default="result"):

    result = {"name": default}

    top = Toplevel()
    top.title("Output name")
    top.geometry("300x150")

    Label(top, text="File name:", font=("Arial", 11)).pack(pady=10)

    entry = Entry(top, width=25)
    entry.insert(0, default)
    entry.pack(pady=5)

    def confirm():
        name = entry.get().strip()
        if name == "":
            name = default
        result["name"] = name
        top.destroy()

    Button(top, text="OK", command=confirm).pack(pady=10)

    top.grab_set()
    top.wait_window()

    return result["name"]


#------SPLIT LABELS & DATA------#

def split_labels_and_data(df):
    label_cols = []
    data_cols = []

    for col in df.columns:
        converted = pd.to_numeric(df[col], errors="coerce")

        if converted.isna().any():
            label_cols.append(col)
        else:
            data_cols.append(col)

    labels_raw = df[label_cols]
    data = df[data_cols].values.astype(float)

    return labels_raw, data


#------MAIN PIPELINE------#

def run_full_pipeline(df, name="result", nad_raw=None, fad_raw=None, pp_raw=None):

    # ---------------- LOGGING SETUP ----------------
    folder = f"best_pca_{name}"
    os.makedirs(folder, exist_ok=True)

    log_path = os.path.join(folder, f"{name}.txt")
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    log_file = open(log_path, "w", encoding="utf-8")
    sys.stdout = log_file
    sys.stderr = log_file

    try:
        # ---------------- PIPELINE ----------------

        labels_raw, data = split_labels_and_data(df)
        conditions = extract_conditions_pipeline(labels_raw)

        # ── LOADING ─────────────────────────────────────────────
        loading = Toplevel()
        loading.title("Processing")
        loading.geometry("250x100")
        loading.resizable(False, False)
        Label(loading, text="Processing, please wait...", font=("Arial", 11)).pack(expand=True)
        loading.update()
        # ────────────────────────────────────────────────────────

        normalizations = {
            "MyFacts":     normalize_myfacts,
            "Intensity":   normalize_intensity,
            "Log1p":       normalize_log1p,
            "L2":          normalize_l2,
        }

        results = {}

        # Subfolders
        all_pca_folder     = os.path.join(folder, "all_pca")
        spectra_folder     = os.path.join(folder, "mean_spectra")
        nadfad_folder      = os.path.join(folder, "nadfad_spectra")
        loadings_folder    = os.path.join(folder, "pc1_loadings")

        print("\n=== RESULTS ===")

        for norm_name, func in normalizations.items():
            try:
                data_norm = func(data)
                X, pca, sep, sil, db, ch, p = run_analysis(data_norm, conditions)
                # --- PC1 loadings plot ---
                save_pc1_loadings_plot(
                    loadings_folder,
                    f"{norm_name}_PC1_loadings.png",
                    pca,
                    f"PC1 loadings — {norm_name}"
                )
                is_bio, rho, pc_used = biological_validation(X, data_norm, conditions)

                print(f"{norm_name} → Sep={sep:.3f} | Sil={sil:.3f} | DB={db:.3f} | "
                      f"CH={ch:.1f} | p={p:.4e} | Spearman({pc_used})={rho:.3f} "
                      f"{'✓' if is_bio else '✗'}")

                # --- PCA plot ---
                save_pca_plot(
                    all_pca_folder,
                    f"{norm_name}.png",
                    X, conditions,
                    f"PCA - {norm_name}  |  ρ({pc_used})={rho:.3f}",
                    p_value=p,
                    pca=pca
                )

                # --- Mean spectra plot (cells) ---
                save_mean_spectra_plot(
                    spectra_folder,
                    f"{norm_name}_spectra.png",
                    data_norm, conditions,
                    f"Mean spectra per condition — {norm_name}"
                )

                # --- NAD/FAD reference spectra plot (if reference provided) ---
                if nad_raw is not None and fad_raw is not None:
                    n_ch_data = data.shape[1]
                    n_ch_ref  = len(nad_raw)

                    if n_ch_ref != n_ch_data:
                        print(f"  [WARNING] NAD/FAD reference has {n_ch_ref} channels "
                              f"but data has {n_ch_data}. Skipping NAD/FAD plot for {norm_name}.")
                    else:
                        nad_norm_vec = apply_normalization_to_vector(nad_raw, norm_name, data)
                        fad_norm_vec = apply_normalization_to_vector(fad_raw, norm_name, data)
                        pp_norm_vec_plot = apply_normalization_to_vector(pp_raw, norm_name, data) if pp_raw is not None else None
                        
                        save_nadfad_spectra_plot(
                            nadfad_folder,
                            f"{norm_name}_nadfad.png",
                            nad_norm_vec, fad_norm_vec,
                            f"NAD / FAD / PP reference spectra — {norm_name}",
                            pp_norm=pp_norm_vec_plot
                        )

                    # --- PCA pintada per score NAD/FAD ---
                    if nad_raw is not None and fad_raw is not None:
                        n_ch_data = data.shape[1]
                        n_ch_ref  = len(nad_raw)
    
                        if n_ch_ref == n_ch_data:
                            nad_norm_vec = apply_normalization_to_vector(nad_raw, norm_name, data)
                            fad_norm_vec = apply_normalization_to_vector(fad_raw, norm_name, data)
                            pp_norm_vec = apply_normalization_to_vector(pp_raw, norm_name, data) if pp_raw is not None else None

                            w_nad, w_fad, w_pp = compute_weights(nad_norm_vec, fad_norm_vec, pp_norm_vec)
    
                            score_nad = data_norm @ w_nad
                            score_fad = data_norm @ w_fad
    
                            colored_folder = os.path.join(folder, "pca_colored")
    
                            save_pca_colored_plot(
                                colored_folder,
                                f"{norm_name}_NAD.png",
                                X, score_nad, "NAD score",
                                f"PCA colored by NAD — {norm_name}",
                                pca=pca,
                                conditions=conditions
                            )
                            save_pca_colored_plot(
                                colored_folder,
                                f"{norm_name}_FAD.png",
                                X, score_fad, "FAD score",
                                f"PCA colored by FAD — {norm_name}",
                                pca=pca,
                                conditions=conditions
                            )
                            if w_pp is not None:
                                score_pp = data_norm @ w_pp
                                save_pca_colored_plot(
                                    colored_folder,
                                    f"{norm_name}_PP.png",
                                    X, score_pp, "PP score",
                                    f"PCA colored by PP — {norm_name}",
                                    pca=pca,
                                    conditions=conditions
                                )
                            save_pca_nadfad_ratio_plot(
                                colored_folder,
                                f"{norm_name}_NADperFAD.png",
                                X, score_nad, score_fad,
                                f"PCA colored by NAD/FAD ratio — {norm_name}",
                                pca=pca,
                                conditions=conditions
                            )

                results[norm_name] = (X, pca, sep, sil, db, ch, p, is_bio, rho)

            except Exception as e:
                print(f"{norm_name} failed: {e}")

        # ---------------- SCORING ---------------- #

        def compute_metrics(sep, sil, db, ch, p):
            return {
                "sep":  sep,
                "sil":  sil,
                "db":   -db,
                "ch":   ch,
                "perm": -np.log10(p + 1e-12) if p < 0.05 else 0
            }

        keys  = ["sep", "sil", "db", "ch", "perm"]
        names_list = list(results.keys())

        M = []
        for nm in names_list:
            _, _, sep, sil, db, ch, p, is_bio, rho = results[nm]
            m = compute_metrics(sep, sil, db, ch, p)
            M.append([m[k] for k in keys])

        M = np.array(M)
        scaler_score = MinMaxScaler()
        M_norm = scaler_score.fit_transform(M)

        weights = np.array([0.30, 0.20, 0.05, 0.15, 0.30])
        scores  = np.dot(M_norm, weights)

        p_values  = [results[n][6] for n in names_list]
        bio_valid = [results[n][7] for n in names_list]

        valid_permanova = [p < 0.05  for p in p_values]
        valid_both      = [perm and bio for perm, bio in zip(valid_permanova, bio_valid)]

        if any(valid_both):
            valid_mask = valid_both
            print("\n Ranking with biological validation ✓")
        elif any(valid_permanova):
            valid_mask = valid_permanova
            print("\nNO normalization passed biological validation. Ranking by PERMANOVA only.")
        else:
            valid_mask = [False] * len(names_list)

        if not any(valid_mask):
            best = None
        else:
            valid_scores = np.array(scores, dtype=float)
            valid_scores[~np.array(valid_mask)] = -np.inf
            best_idx = np.argmax(valid_scores)
            best = names_list[best_idx]

        print(f"\n Best normalization: {best}")

        if best is not None:
            rho_best = results[best][8]
            bio_best = results[best][7]
            print(f" Spearman NAD/FAD: ρ={rho_best:.3f} | "
                  f"Biological validation: {'PASSED ✓' if bio_best else 'NOT PASSED ✗'}")

        print("\n=== NORMALIZED METRICS (0–1) ===")
        for i, nm in enumerate(names_list):
            print(f"{nm}: {M_norm[i]}")

        print("\n=== SCORE BREAKDOWN ===")
        for i, nm in enumerate(names_list):
            print(f"{nm} score = {scores[i]:.3f}")

        if best is not None:
            winner_idx = names_list.index(best)
            diff = M_norm[winner_idx] - np.mean(M_norm, axis=0)

            print("\n=== WHY THIS WON ===")
            metric_names = ["Separation", "Silhouette", "DB(inv)", "Calinski", "PERMANOVA"]
            for i in np.argsort(-diff):
                if diff[i] > 0:
                    print(f"+ {metric_names[i]} advantage: {diff[i]:.2f}")

            X_best, pca_best = results[best][0], results[best][1]
            best_folder = os.path.join(folder, "best_pca")
            p_best = results[best][6]

            save_pca_plot(
                best_folder,
                f"{best}.png",
                X_best,
                conditions,
                f"PCA - BEST ({best})",
                p_value=p_best,
                pca=pca_best
            )

        else:
            print("\n No valid PCA")
        
        loading.destroy()
        messagebox.showinfo(
            "Process finished",
            f"Results saved in folder:\n{folder}"
        )

        return best

    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_file.close()


#------ASK FILE UI FUNCTION------#

def openFile():
    filepath = filedialog.askopenfilename(
        filetypes=[
            ("Excel files", "*.xlsx *.xls"),
            ("CSV files", "*.csv"),
            ("All files", "*.*")
        ]
    )

    if not filepath:
        return

    # Ask for optional NAD/FAD reference file
    nad_raw, fad_raw, pp_raw = None, None, None
    use_ref = messagebox.askyesno(
        "Metabolites reference",
        "Do you want to load a metabolites reference file?\n"
    )
    if use_ref:
        ref_path = filedialog.askopenfilename(
            title="Select metabolites reference file",
            filetypes=[
                ("Excel files", "*.xlsx *.xls"),
                ("CSV files", "*.csv"),
                ("All files", "*.*")
            ]
        )
        if ref_path:
            try:
                nad_raw, fad_raw, pp_raw = load_nadfad_reference(ref_path)
                pp_msg = f" + PP" if pp_raw is not None else " (PP not found)"
                messagebox.showinfo("Reference loaded",
                                    f"NAD + FAD{pp_msg} reference loaded ({len(nad_raw)} channels).")
            except Exception as e:
                messagebox.showerror("Error loading reference", str(e))

    data = load_data(filepath)

    if isinstance(data, dict):
        for sheet_name, df in data.items():
            print(f"\n--- Processing: {sheet_name} ---")
            user_name = ask_filename_gui(sheet_name)
            run_full_pipeline(df, user_name, nad_raw=nad_raw, fad_raw=fad_raw, pp_raw=pp_raw)
    else:
        print("\n--- Processing CSV ---")
        user_name = ask_filename_gui("csv_file")
        run_full_pipeline(df, user_name, nad_raw=nad_raw, fad_raw=fad_raw, pp_raw=pp_raw)


#------UI WINDOW------#

window = tk.Tk()
window.title("Best Norm & PCA")
window.geometry("350x200")
window.resizable(False, False)

button = Button(window, text="Select file", command=openFile)
button.pack(expand=True, padx=20, pady=20)

exit_button = Button(window, text="Exit", command=window.destroy)
exit_button.pack(expand=True, padx=20, pady=5)

window.mainloop()


# In[ ]:




