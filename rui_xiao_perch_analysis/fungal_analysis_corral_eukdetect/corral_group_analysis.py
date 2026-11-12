#!/usr/bin/env python3
"""
CORRAL eukaryote CPM comparison: Control vs Survived vs Died (PERCH Mali)

Inputs (in /home/user/perch/corral):
  - 20260313_Mali_PERCH_mapping_file.tsv : metadata (SampleID, study_group, study_group_collapse, death)
  - busco_taxid_link.txt                 : marker -> taxid (2 col); species parsed from marker name
  - cpm.matrix.tsv                       : taxa (rows) x samples (cols), CPM; blank/NA = 0

Groups:
  Control  : study_group_collapse == "Community control"
  Survived : study_group == "Very severe case" & death == "No"
  Died     : study_group == "Very severe case" & death == "Yes"

Colors:
  Control  #BA616C
  Survived #2767B0
  Died     #BEC100

Outputs (in /home/user/perch/corral/analysis):
  group_assignments.tsv, taxon_names.tsv, taxon_group_stats.tsv,
  diff_abundance.tsv, top_signatures.tsv, SUMMARY.md, fig1..fig5 .png
"""
import os, re, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SRC = "/home/user/perch/corral"
OUT = os.path.join(SRC, "analysis")
os.makedirs(OUT, exist_ok=True)

COLORS = {"Control": "#BA616C", "Survived": "#2767B0", "Died": "#BEC100"}
GROUPS = ["Control", "Survived", "Died"]

# ---------------------------------------------------------------- load data
print("loading CPM matrix ...")
cpm_raw = pd.read_csv(os.path.join(SRC, "cpm.matrix.tsv"), sep="\t",
                      index_col=0, dtype=object)
cpm = cpm_raw.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(float)
cpm.index = [str(x) for x in cpm_raw.index]
is_unk_rows = np.asarray(cpm.index.str.startswith("?"))  # marker_alignments "unknown taxon" clusters
print(f"  {cpm.shape[0]} taxa x {cpm.shape[1]} samples, "
      f"{int((cpm.values > 0).sum())} positive cells, max CPM {cpm.values.max():.2f} "
      f"({int(is_unk_rows.sum())} rows are '?' unknown-taxon clusters, "
      f"{int((~is_unk_rows).sum())} identified taxa)")

print("loading metadata ...")
meta = pd.read_csv(os.path.join(SRC, "20260313_Mali_PERCH_mapping_file.tsv"),
                   sep="\t", dtype=str)

def assign(r):
    if r.study_group_collapse == "Community control":
        return "Control"
    if r.study_group == "Very severe case":
        if r.death == "No":  return "Survived"
        if r.death == "Yes": return "Died"
        return "Very severe case (death NA)"
    if r.study_group == "Severe case":
        return "Severe case (not in 3 groups)"
    return "Unknown group"

meta["group"] = meta.apply(assign, axis=1)
meta_i = meta.set_index("SampleID")
g = meta_i[meta_i.group.isin(GROUPS)]["group"].reindex(cpm.columns)
n_grp = g.value_counts()
print("group sizes (in matrix):")
for gr in GROUPS:
    print(f"  {gr:9s} n={n_grp.get(gr,0)}")
print(f"  matrix samples outside the 3 groups (Severe case / unclassified): {int(g.isna().sum())}")

# ------------------------------------------------------------- taxon names
print("building taxid -> species from busco_taxid_link.txt ...")
link = pd.read_csv(os.path.join(SRC, "busco_taxid_link.txt"),
                   sep="\t", header=None, names=["marker", "taxid"],
                   dtype=str, engine="python")
link["taxid"] = link["taxid"].str.strip()

rx_norm = re.compile(r"^([A-Za-z]+)-(.+?)-(\d+at2759-[A-Z]\d+)$")   # domain prefix is capitalised (Fungi-, Protist-...)
rx_coll = re.compile(r"^([A-Za-z]+)-(.*)Collapse")

norm = link["marker"].str.extract(rx_norm.as_pattern() if hasattr(rx_norm, 'as_pattern') else rx_norm.pattern)
coll = link["marker"].str.extract(rx_coll.pattern)
link["domain"] = norm[0].fillna(coll[0]).fillna("?")
link["species"] = norm[1].fillna(coll[1]).fillna(link["marker"])
link["is_coll"] = norm[1].isna().astype(int)

# prefer non-collapsed markers when naming a taxid
link = link.sort_values("is_coll")
names = (link.drop_duplicates("taxid")
           .set_index("taxid")[["species", "domain"]]
           .rename(columns={"species": "taxon_name", "domain": "domain"}))
names["taxon_name"] = names["taxon_name"].str.replace("_", " ", regex=False)
missing = [t for t in cpm.index if not t.startswith("?") and t not in names.index]
for t in missing:
    names.loc[t, ["taxon_name", "domain"]] = [f"taxid_{t}", "?"]
    print(f"  WARN taxid {t} not found in link file")
# label the '?' unknown-taxon clusters from their candidate taxids
for lab in cpm.index[cpm.index.str.startswith("?")]:
    ids = lab.lstrip("?").split(",")
    cn = [names.loc[i, "taxon_name"] for i in ids if i in names.index]
    cn = list(dict.fromkeys(cn))
    if len(cn) > 3:
        disp = f"unknown cluster ({', '.join(cn[:3])} +{len(cn)-3} more)"
    elif cn:
        disp = f"unknown cluster ({', '.join(cn)})"
    else:
        disp = "unknown cluster"
    names.loc[lab, ["taxon_name", "domain"]] = [disp, "unknown cluster"]
names = names.reindex(cpm.index)
names.to_csv(os.path.join(OUT, "taxon_names.tsv"), sep="\t")
cpm.index = [f"{r}  {names.loc[r,'taxon_name']} ({names.loc[r,'domain']})"
             for r in cpm.index]

# ------------------------------------------------------------------- stats
cpm_t = cpm.T  # samples x taxa
gmat = g

def group_block(cols):
    sel = [c for c in cpm.columns if gmat[c] in cols]
    return sel

stats_rows = []
for txid in cpm.index:
    for gr in GROUPS:
        x = cpm.loc[txid, gmat == gr].values
        nz = x > 0
        stats_rows.append({
            "taxid": txid.split()[0],
            "taxon": names.loc[txid.split()[0], "taxon_name"],
            "domain": names.loc[txid.split()[0], "domain"],
            "group": gr,
            "n_samples": int(x.size),
            "n_prevalent": int(nz.sum()),
            "prevalence_pct": 100.0 * nz.mean(),
            "mean_cpm": float(x.mean()),
            "median_cpm": float(np.median(x)),
        })
tgs = pd.DataFrame(stats_rows)
tgs.to_csv(os.path.join(OUT, "taxon_group_stats.tsv"), sep="\t", index=False)

# ------------------------------- differential (pairwise Wilcoxon + BH FDR)
pseudo = 0.1
log1p = np.log1p(cpm.values.T)          # samples x taxa
pairs = [(a, b) for a in GROUPS for b in GROUPS if a != b]   # all 6 ordered pairs

def bh(p):
    """Benjamini-Hochberg FDR, matches R p.adjust(method='BH')."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    raw = p[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(raw[::-1])[::-1]   # running min ALONG sorted-p order
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out

diff_rows = []
for a, b in pairs:
    xa = log1p[gmat == a]                # samples_a x taxa
    xb = log1p[gmat == b]                # samples_b x taxa
    na, nb = xa.shape[0], xb.shape[0]
    mw_p, fe_p = [], []
    for j in range(cpm.shape[0]):        # j indexes TAXA (rows)
        va, vb = xa[:, j], xb[:, j]
        # Mann-Whitney U (= Wilcoxon rank-sum) on log1p CPM; independent groups
        if np.unique(va).size == 1 and np.unique(vb).size == 1:
            mw_p.append(1.0)             # identical constants in both groups
        else:
            try:
                mw_p.append(float(stats.mannwhitneyu(va, vb,
                                                     alternative="two-sided").pvalue))
            except Exception as e:
                print(f"  WARN MW failed for taxon {j}: {e}")
                mw_p.append(1.0)
        # Fisher exact on prevalence (present/absent) - powerful for sparse data
        pa = int((va > 0).sum())
        pb = int((vb > 0).sum())
        table = np.array([[pa, na - pa], [pb, nb - pb]])
        fe_p.append(float(stats.fisher_exact(table).pvalue))
    mw_fdr = bh(np.array(mw_p))
    fe_fdr = bh(np.array(fe_p))
    for j in range(cpm.shape[0]):
        tid = cpm.index[j].split()[0]
        diff_rows.append({
            "pair": f"{a}_vs_{b}", "taxid": tid,
            "taxon": names.loc[tid, "taxon_name"],
            "domain": names.loc[tid, "domain"],
            "mean_cpm_" + a: float(cpm.values.T[gmat == a, j].mean()),
            "mean_cpm_" + b: float(cpm.values.T[gmat == b, j].mean()),
            "n_prevalent_" + a: int((cpm.values.T[gmat == a, j] > 0).sum()),
            "n_prevalent_" + b: int((cpm.values.T[gmat == b, j] > 0).sum()),
            "log2FC_a_vs_b": float(np.log2(cpm.values.T[gmat == a, j].mean() + pseudo)
                                   - np.log2(cpm.values.T[gmat == b, j].mean() + pseudo)),
            "mannwhitney_p": float(mw_p[j]),
            "fdr_bh": float(mw_fdr[j]),
            "fisher_p": float(fe_p[j]),
            "fisher_fdr": float(fe_fdr[j]),
        })
diff = pd.DataFrame(diff_rows)
diff.to_csv(os.path.join(OUT, "diff_abundance.tsv"), sep="\t", index=False)

# ------------------------------------------------- top signatures per group
sig_rows = []
for gr in GROUPS:
    sub = tgs[tgs.group == gr].copy()
    sub["mean_rank"] = sub["mean_cpm"].rank(ascending=False)
    sub["prev_rank"] = sub["prevalence_pct"].rank(ascending=False)
    # specificity: prevalence in this group vs max prevalence in other groups
    spec = np.ones(len(sub))
    for i, (_, row) in enumerate(sub.iterrows()):
        others = tgs[(tgs.taxon == row.taxon) & (tgs.group != gr)]
        spec[i] = row.prevalence_pct / max(others["prevalence_pct"].max(), 1.0)
    sub["specificity"] = spec
    sub["score"] = sub.mean_rank + sub.prev_rank
    top = sub.sort_values("score").head(10)
    for _, r in top.iterrows():
        sig_rows.append({"group": gr, "taxid": r.taxid, "taxon": r.taxon,
                         "domain": r.domain, "mean_cpm": r.mean_cpm,
                         "prevalence_pct": r.prevalence_pct,
                         "specificity_vs_others": r.specificity,
                         "rank_by_score": len(sig_rows) % 10 + 1})
sig = pd.DataFrame(sig_rows)
sig.to_csv(os.path.join(OUT, "top_signatures.tsv"), sep="\t", index=False)

# ------------------------------------------------------------- assignments
meta_out = meta[["SampleID", "study_group", "study_group_collapse", "death", "group"]].copy()
meta_out.to_csv(os.path.join(OUT, "group_assignments.tsv"), sep="\t", index=False)

# ------------------------------------------------------------------- FIG 1
# sample-level heatmap of top variable taxa (ordered by group)
# one-way F per taxon (between-group variation on log1p CPM)
fvals = []
for j in range(cpm.shape[0]):        # j indexes TAXA (rows)
    grs = [np.log1p(cpm.values.T[gmat == gr, j]) for gr in GROUPS]
    try:
        fvals.append(stats.f_oneway(*grs)[0])
    except Exception:
        fvals.append(0.0)
order = np.argsort(fvals)[::-1][:25]
cols_ordered = []
for gr in GROUPS:                             # Control, Survived, Died
    sel = [c for c in cpm.columns if gmat[c] == gr]
    sel.sort(key=lambda c: -cpm.loc[:, [c]].values.sum())
    cols_ordered += sel
pos = np.array([cpm.columns.get_loc(c) for c in cols_ordered])
Hm = np.log1p(cpm.values[order][:, pos])      # (25 taxa, n grouped samples)
Hm = Hm - Hm.mean(axis=1, keepdims=True)      # per-taxon (row) centering
Hm = np.clip(Hm, Hm.min(), np.percentile(Hm, 99.5))

fig, ax = plt.subplots(figsize=(14, 9))
im = ax.imshow(Hm, aspect="auto", cmap="Reds", origin="upper")
ax.set_yticks(range(Hm.shape[0]))
ax.set_yticklabels([cpm.index[i].split("  ")[1].split(" (")[0] for i in order],
                   fontsize=8)
# group strip below the image
x = 0
for gr in GROUPS:
    w = int((gmat == gr).sum())
    ax.add_patch(plt.Rectangle((x + 0.5, Hm.shape[0] + 1.2), max(w - 1, 1), 0.8,
                               facecolor=COLORS[gr], edgecolor="none",
                               clip_on=False))
    ax.text(x + w / 2, Hm.shape[0] + 2.5, gr,
            ha="center", fontsize=12, color=COLORS[gr], fontweight="bold",
            clip_on=False)
    x += w
ax.set_xticks([])
ax.set_title("Top 25 taxa by between-group variation (one-way F on log1p CPM)\n"
             "samples ordered by group, then total eukaryotic CPM (taxon-centred, log1p CPM)",
             fontsize=12)
cb = plt.colorbar(im, fraction=0.02, pad=0.01)
cb.set_label("log1p(CPM), taxon-centred", fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig1_heatmap_top_taxa.png"), dpi=200)
plt.close()
print("fig1 saved")

# ------------------------------------------------------------------- FIG 2
# boxplots of top 12 taxa by max group mean
maxgm = {r: max(cpm.loc[r, gmat == gr].mean() for gr in GROUPS) for r in cpm.index}
top12 = sorted(cpm.index, key=lambda r: -maxgm[r])[:12]
fig, axes = plt.subplots(4, 3, figsize=(13, 15))
for ax, r in zip(axes.ravel(), top12):
    xvals, groups = [], []
    for i, gr in enumerate(GROUPS):
        xvals.append(np.log1p(cpm.loc[r, gmat == gr].values))
        groups.append(gr)
    bp = ax.boxplot(xvals, tick_labels=groups, patch_artist=True, showfliers=False,
                    widths=0.55)
    for i, (art, gr) in enumerate(zip(bp["boxes"], groups), 1):
        art.set_facecolor(COLORS[gr]); art.set_alpha(0.85)
        art.set_edgecolor("0.3")
    for i, gr in enumerate(groups, 1):
        xs = np.random.default_rng(i).normal(i, 0.045, len(xvals[i-1]))
        ax.scatter(xs, xvals[i-1], s=7, color=COLORS[gr], alpha=0.45, edgecolor="none")
    mean_c = {gr: cpm.loc[r, gmat == gr].mean() for gr in GROUPS}
    prev_c = {gr: float((cpm.loc[r, gmat == gr].values > 0).mean() * 100) for gr in GROUPS}
    ax.set_title(f"{names.loc[r.split()[0],'taxon_name']}\n"
                 f"mean CPM: C {mean_c['Control']:.2f} | S {mean_c['Survived']:.2f} | "
                 f"D {mean_c['Died']:.2f}", fontsize=9)
    ax.set_ylabel("log1p(CPM)", fontsize=8)
for ax in axes.ravel():
    ax.grid(axis="y", alpha=0.25)
fig.suptitle("Eukaryote CPM by clinical group - top 12 taxa by maximum group mean CPM",
             fontsize=13, y=0.995)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig2_boxplots_key_taxa.png"), dpi=200)
plt.close()
print("fig2 saved")

# ------------------------------------------------------------------- FIG 3
# top-10 taxa per group (horizontal bars)
fig, axes = plt.subplots(1, 3, figsize=(16, 9), sharey=False)
for ax, gr in zip(axes, GROUPS):
    sel = sig[sig.group == gr]
    sel = sel.iloc[::-1]   # rank 1 on top
    y = np.arange(len(sel))[::-1]
    ax.barh(y, sel.mean_cpm, color=COLORS[gr], alpha=0.9, height=0.7)
    labels = [f"{r.taxon}  ({r.prevalence_pct:.0f}%)" for r in sel.itertuples()]
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("mean CPM in group", fontsize=9)
    ax.set_title(f"{gr} (n={n_grp[gr]}) - top 10 by mean CPM",
                 color=COLORS[gr], fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.25)
    for i, v in enumerate(sel.mean_cpm):
        ax.text(v + 0.02 * sel.mean_cpm.max(), y[i], f"{v:.2f}", va="center", fontsize=7)
fig.suptitle("Top 10 eukaryotic taxa per clinical group (prevalence in parentheses)",
             fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig3_group_top_taxa.png"), dpi=200)
plt.close()
print("fig3 saved")

# ------------------------------------------------------------------- FIG 4
# richness + total eukaryotic CPM per sample
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
rich = (cpm.T > 0).sum(axis=1)
tot = cpm.T.sum(axis=1)
data = [rich.values[gmat.values == gr] for gr in GROUPS]
bp = axes[0].boxplot(data, tick_labels=[f"{gr}\nn={n_grp[gr]}" for gr in GROUPS],
                     patch_artist=True, showfliers=False, widths=0.55)
for art, gr in zip(bp["boxes"], GROUPS):
    art.set_facecolor(COLORS[gr]); art.set_alpha(0.85); art.set_edgecolor("0.3")
axes[0].set_title("Number of eukaryotic taxa detected per sample")
axes[0].set_ylabel("taxa with CPM > 0"); axes[0].grid(axis="y", alpha=0.25)
data2 = [np.log1p(tot.values[gmat.values == gr]) for gr in GROUPS]
bp2 = axes[1].boxplot(data2, tick_labels=[f"{gr}\nn={n_grp[gr]}" for gr in GROUPS],
                      patch_artist=True, showfliers=False, widths=0.55)
for art, gr in zip(bp2["boxes"], GROUPS):
    art.set_facecolor(COLORS[gr]); art.set_alpha(0.85); art.set_edgecolor("0.3")
axes[1].set_title("Total eukaryotic CPM per sample (log1p)")
axes[1].set_ylabel("log1p(sum of CPM)"); axes[1].grid(axis="y", alpha=0.25)
plt.suptitle("Eukaryotic load in PERCH sputum samples by clinical group", fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig4_richness.png"), dpi=200)
plt.close()
print("fig4 saved")

# ------------------------------------------------------------------- FIG 5
# pairwise differential: top |log2FC| taxa, FDR highlighted
fig, axes = plt.subplots(1, 3, figsize=(16, 8))
short = {"Died_vs_Survived": "Died vs Survived",
         "Died_vs_Control": "Died vs Control",
         "Survived_vs_Control": "Survived vs Control"}
FIG5_PAIRS = [("Died", "Survived"), ("Died", "Control"), ("Survived", "Control")]
for ax, (a, b) in zip(axes, FIG5_PAIRS):
    sub = diff[diff.pair == f"{a}_vs_{b}"].copy()
    top = sub.reindex(sub["log2FC_a_vs_b"].abs().sort_values(ascending=False).index).head(15)
    top = top.iloc[::-1]
    y = np.arange(len(top))
    cols = [COLORS[a] if v > 0 else COLORS[b] for v in top["log2FC_a_vs_b"]]
    ax.barh(y, top["log2FC_a_vs_b"], color=cols, height=0.7, alpha=0.9)
    sig_mask = top["fdr_bh"] < 0.05
    ax.plot([0] * len(top), y, color="0.5", lw=0.8, zorder=0)
    labels = [f"{r.taxon}{'*' if r.fdr_bh < 0.05 else ''}" for r in top.itertuples()]
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("log2FC (mean CPM, pseudocount 0.1)")
    ax.set_title(f"{short[f'{a}_vs_{b}']}\n"
                 f"{int(sig_mask.sum())} taxa FDR<0.05 (of {cpm.shape[0]})",
                 fontsize=11)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(handles=[Patch(facecolor=COLORS[a], label=f"higher in {a}"),
                       Patch(facecolor=COLORS[b], label=f"higher in {b}")],
              fontsize=8, loc="lower right")
fig.suptitle("Top 15 taxa by |log2FC| per pairwise comparison "
             "(* = BH FDR < 0.05, Mann-Whitney U on log1p CPM)", fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig5_pairwise_diff.png"), dpi=200)
plt.close()
print("fig5 saved")

# ---------------------------------------------------------------- summary
def top_taxa_table(gr, n=10):
    s = tgs[tgs.group == gr].copy()
    s["score"] = s["mean_cpm"].rank(ascending=False) + s["prevalence_pct"].rank(ascending=False)
    s = s.sort_values("score").head(n)
    lines = [f"{i}. {r.taxon} ({r.domain}) - mean CPM {r.mean_cpm:.2f}, "
             f"prevalence {r.prevalence_pct:.0f}% ({r.n_prevalent}/{r.n_samples})"
             for i, (_, r) in enumerate(s.iterrows(), 1)]
    return "\n".join(lines)

def sig_taxa(gr, n=5):
    rows = []
    for a, b in [(gr, x) for x in GROUPS if x != gr]:
        sub = diff[diff.pair == f"{a}_vs_{b}"]
        sub = sub[(sub.fdr_bh < 0.05) & (sub.log2FC_a_vs_b > 0)]
        sub = sub.sort_values("fdr_bh")
        for r in sub.head(n).itertuples():
            rows.append((b, r.taxon, r.fdr_bh, r.log2FC_a_vs_b, r.fisher_fdr))
    rows.sort(key=lambda t: t[2])
    return "\n".join(f"   - {t[1]}  ({t[0]}: log2FC {t[3]:+.1f}, "
                     f"MW FDR {t[2]:.2e}, Fisher FDR {t[4]:.2e})"
                     for t in rows[:n]) if rows else "   - none at FDR<0.05"

med_known = {gr: float(np.median(cpm[~is_unk_rows].T[gmat.values == gr].sum(axis=1).values))
              for gr in GROUPS}
med_unk = {gr: float(np.median(cpm[is_unk_rows].T[gmat.values == gr].sum(axis=1).values))
           for gr in GROUPS}

summary = f"""# CORRAL eukaryote CPM analysis - PERCH Mali (Control / Survived / Died)

Generated: {pd.Timestamp.now():%Y-%m-%d %H:%M}
Data: {cpm.shape[0]} eukaryotic taxa x {cpm.shape[1]} samples
(int {int((cpm.values>0).sum())} positive CPM cells; empty cells = no detection)

## Groups (from 20260313_Mali_PERCH_mapping_file.tsv)
| Group | Definition | n (in CPM matrix) |
|---|---|---|
| Control | study_group_collapse = "Community control" | {n_grp.get('Control',0)} |
| Survived | study_group = "Very severe case" & death = "No" | {n_grp.get('Survived',0)} |
| Died | study_group = "Very severe case" & death = "Yes" | {n_grp.get('Died',0)} |

(Excluded from comparison: {int((meta['study_group']=='Severe case').sum())} "Severe case"
and {int(meta['study_group'].isna().sum())} unclassified samples that are in the metadata
but outside the requested 3 groups.)

## Eukaryotic load (per sample)
- Known vs unknown signal (median total CPM per sample):
  - Identified taxa ({int((~is_unk_rows).sum())}): Control {med_known['Control']:.2f} | Survived {med_known['Survived']:.2f} | Died {med_known['Died']:.2f}
  - 'Unknown taxon' clusters ({int(is_unk_rows.sum())}): Control {med_unk['Control']:.2f} | Survived {med_unk['Survived']:.2f} | Died {med_unk['Died']:.2f}
- Taxa detected: Control median {np.median(rich.values[gmat.values=='Control']):.0f} |
  Survived median {np.median(rich.values[gmat.values=='Survived']):.0f} |
  Died median {np.median(rich.values[gmat.values=='Died']):.0f}
- Total eukaryotic CPM: Control median {np.median(np.log1p(tot.values[gmat.values=='Control'])):.2f} (log1p) |
  Survived {np.median(np.log1p(tot.values[gmat.values=='Survived'])):.2f} |
  Died {np.median(np.log1p(tot.values[gmat.values=='Died'])):.2f}

## Top 10 taxa per group (rank by mean CPM + prevalence)
### Control
{top_taxa_table('Control')}
### Survived
{top_taxa_table('Survived')}
### Died
{top_taxa_table('Died')}

## Statistically significant enrichments
(Tests: Mann-Whitney U on log1p CPM; Fisher exact on prevalence; both BH-corrected. Shown at FDR<0.05 on Mann-Whitney.)
### Enriched in Died
{sig_taxa('Died')}
### Enriched in Survived
{sig_taxa('Survived')}
### Enriched in Control
{sig_taxa('Control')}

## Files
- group_assignments.tsv - per-sample group labels
- taxon_names.tsv - taxid -> species (from busco_taxid_link.txt marker names)
- taxon_group_stats.tsv - prevalence / mean / median CPM per taxon x group
- diff_abundance.tsv - pairwise Mann-Whitney U + Fisher exact on prevalence + log2FC + BH FDR (all 84 taxa x 6 ordered pairs)
- top_signatures.tsv - top-10 signature taxa per group
- fig1_heatmap_top_taxa.png - top-25 variable taxa heatmap, samples ordered by group
- fig2_boxplots_key_taxa.png - boxplots of the 12 most abundant taxa per group
- fig3_group_top_taxa.png - top-10 taxa per group
- fig4_richness.png - taxa richness and total eukaryotic CPM per group
- fig5_pairwise_diff.png - top differential taxa per pairwise comparison
"""
with open(os.path.join(OUT, "SUMMARY.md"), "w") as f:
    f.write(summary)
print("SUMMARY.md written")
print("ALL DONE")
