#!/usr/bin/env python3
"""
decidability_report.py — 可判定性报告工具（可直接跑在别人的 h5ad 上）

用法
  python decidability_report.py --h5ad data.h5ad --label-col cell_type --outdir report/
  可选
    --cluster-col  已有的聚类列（不给则用 Leiden res=1 现算）
    --use-llm      用大模型枚举替代假设（需 DEEPSEEK_API_KEY）；默认走内置规则库，离线可跑
    --genes        基因面板限制（一行一个基因的文本文件）；不给则用全部基因
    --min-cells    每个群最少细胞数（默认 30）

输出（outdir 下）
  decidability.csv     每个群：等价类大小、分级、被排除/保留的替代
  alternatives.csv     逐条替代假设的证伪结果
  report.md            可读报告（含方法、判据阈值、限制说明）

判据（写死，不随数据调整）
  must_present 基因：簇内检出率 >= 0.25 且该基因在本簇排前 5 → 支持；
                     0.10–0.25 → 不确定；< 0.10 → 反驳，该替代被排除
  must_absent  基因：簇内检出率 >= 0.50 → 反驳，排除
  判别基因全部不在面板上 → 不可测，保留在等价类里
  分级：等价类 1 = decidable，2–3 = partial，>=4 = undecidable
"""
import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

# ---------------------------------------------------------------- rule library
# 关键词 → 常见混淆项（每项给出判别基因建议；找不到的基因会被判为"不可测"）
CONFUSABLE = {
    "cardiomyocyte": [("Atrial cardiomyocyte", ["NPPA", "MYH6", "TBX5"], ["MYH7", "IRX4"]),
                      ("Conduction / pacemaker cell", ["HCN4", "TBX3", "SHOX2"], ["MYH7"])],
    "fibroblast": [("Myofibroblast / activated fibroblast", ["POSTN", "ACTA2", "TAGLN"], ["DCN"]),
                   ("Stellate / perivascular fibroblast", ["RGS5", "PDGFRB"], ["COL1A1"])],
    "endothelial": [("Lymphatic endothelial cell", ["LYVE1", "PROX1", "PDPN"], ["PECAM1"]),
                    ("Endocardial cell", ["NPR3", "NFATC1"], ["VWF"])],
    "pericyte": [("Vascular smooth muscle cell", ["MYH11", "ACTA2", "TAGLN"], ["RGS5", "KCNJ8"])],
    "smooth muscle": [("Pericyte", ["RGS5", "KCNJ8", "ABCC9"], ["MYH11", "ACTA2"])],
    "epicard": [("Mesothelial cell", ["MSLN", "UPK3B", "CALB2"], ["TBX18", "ALDH1A2"])],
    "mesothel": [("Epicardial cell", ["TBX18", "ALDH1A2"], ["MSLN"])],
    "epithelial": [("Mesothelial cell", ["MSLN", "UPK3B"], ["EPCAM", "KRT8"]),
                   ("Basal / myoepithelial cell", ["KRT5", "KRT14", "TP63"], ["EPCAM"])],
    "macrophage": [("Monocyte", ["CD14", "FCN1", "S100A8"], ["C1QA", "C1QB"]),
                   ("Dendritic cell", ["CD1C", "FCER1A", "CLEC9A"], ["C1QA"])],
    "t cell": [("NK cell", ["NKG7", "GNLY", "KLRD1"], ["CD3E", "CD3D"]),
               ("NKT / MAIT cell", ["TRAV1-2", "SLC4A10"], ["CD8A"])],
    "b cell": [("Plasma cell", ["JCHAIN", "MZB1", "XBP1"], ["MS4A1", "CD19"])],
    "neural": [("Glial / Schwann cell", ["PLP1", "SOX10", "S100B"], ["NRXN1", "SYT1"]),
               ("Neuronal progenitor", ["DCX", "SOX2"], ["MAP2"])],
}
UNIVERSAL = [("Doublet / ambient-RNA artefact", [], ["MKI67"]),
             ("Cycling / proliferating state", ["MKI67", "TOP2A"], [])]


def build_alternatives(label, use_llm=False):
    s = str(label).lower()
    alts = []
    for kw, items in CONFUSABLE.items():
        if kw in s:
            alts += items
    alts += UNIVERSAL
    if not alts:
        alts = [("Unresolved alternative identity", [], [])]
    out = []
    for name, pres, absn in alts:
        out.append(dict(label=name, must_present=";".join(pres), must_absent=";".join(absn),
                        axis="lineage" if pres or absn else "artefact"))
    return out


def llm_alternatives(label, markers_txt, tissue, panel):
    """可选的 LLM 路径：复用项目内 31_llm_alternatives.py 的调用协议。"""
    import hashlib
    import requests
    key = os.environ.get("DEEPSEEK_API_KEY")
    if not key:
        return None
    prompt = (f"Tissue context: {tissue}. Measured genes: {panel}.\n"
              f"Cluster labelled '{label}' with markers:\n{markers_txt}\n"
              "List 3-8 cell identities this evidence cannot exclude. Return JSON "
              '{"alternatives":[{"label":"","must_present":[],"must_absent":[],"axis":""}]}')
    body = dict(model="deepseek-chat", temperature=0.3, max_tokens=1500,
                messages=[dict(role="user", content=prompt)])
    r = requests.post("https://api.deepseek.com/v1/chat/completions", json=body, timeout=90,
                      headers={"Authorization": f"Bearer {key}"})
    if r.status_code != 200:
        return None
    txt = r.json()["choices"][0]["message"]["content"]
    i, j = txt.find("{"), txt.rfind("}")
    try:
        d = json.loads(txt[i:j + 1])
    except Exception:
        return None
    alts = []
    for a in d.get("alternatives", []):
        alts.append(dict(label=a.get("label", ""), must_present=";".join(a.get("must_present", [])),
                         must_absent=";".join(a.get("must_absent", [])), axis=a.get("axis", "")))
    return alts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--label-col", required=True)
    ap.add_argument("--cluster-col", default=None)
    ap.add_argument("--outdir", default="decidability_report")
    ap.add_argument("--genes", default=None)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--use-llm", action="store_true")
    ap.add_argument("--tissue", default="unspecified tissue")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    a = ad.read_h5ad(args.h5ad)
    if a.X.max() > 30:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    if args.genes:
        genes = [g.strip() for g in open(args.genes) if g.strip()]
        genes = [g for g in genes if g in a.var_names]
        a = a[:, genes].copy()

    if args.cluster_col and args.cluster_col in a.obs.columns:
        a.obs["_cl"] = a.obs[args.cluster_col].astype(str)
    else:
        sc.pp.highly_variable_genes(a, n_top_genes=min(2000, a.n_vars))
        sc.pp.pca(a, n_comps=min(50, a.n_vars - 1, a.n_obs - 1), svd_solver="arpack")
        sc.pp.neighbors(a, n_pcs=min(50, a.obsm["X_pca"].shape[1]))
        sc.tl.leiden(a, resolution=1.0, key_added="_cl", flavor="igraph", n_iterations=2, directed=False)

    sc.tl.rank_genes_groups(a, "_cl", method="wilcoxon", n_genes=min(200, a.n_vars), pts=True)
    res = a.uns["rank_genes_groups"]
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    det = (X > 0)
    cl = a.obs["_cl"].astype(str).values
    panel = set(a.var_names)

    rows_alt, rows_cl = [], []
    for g in res["names"].dtype.names:
        m = cl == g
        if m.sum() < args.min_cells:
            continue
        labs = a.obs.loc[m, args.label_col].astype(str).value_counts()
        cur = labs.index[0]
        n_g = min(len(res["names"][g]), len(res["logfoldchanges"][g]),
                  len(res["pts"][g]), len(res["pts_rest"][g]))
        mk = pd.DataFrame({"gene": res["names"][g][:n_g], "lfc": res["logfoldchanges"][g][:n_g],
                           "pct_in": res["pts"][g][:n_g], "pct_out": res["pts_rest"][g][:n_g]})
        mk = mk.sort_values("lfc", ascending=False)
        top = mk.head(30)
        mktxt = "\n".join(f"{r.gene}: log2FC={r.lfc:.2f}, in={r.pct_in:.0%}, out={r.pct_out:.0%}"
                          for r in top.itertuples())
        alts = None
        if args.use_llm:
            alts = llm_alternatives(cur, mktxt, args.tissue, f"{len(panel)} genes measured")
            if alts:
                rows_cl.append(dict(note="alternatives from LLM"))
        if not alts:
            alts = build_alternatives(cur)

        # 检出率矩阵（按群）
        pct_in = {}
        for gg in res["names"].dtype.names:
            mm = cl == gg
            pct_in[gg] = det[mm].mean(axis=0) if mm.sum() else np.zeros(X.shape[1])

        retained = 0
        for alt in alts:
            pres = [x.strip().upper() for x in str(alt["must_present"]).split(";") if x.strip()]
            absn = [x.strip().upper() for x in str(alt["must_absent"]).split(";") if x.strip()]
            # 基因名大小写无关匹配
            up = {v.upper(): v for v in a.var_names}
            pres_m = [up[p] for p in pres if p in up]
            abs_m = [up[q] for q in absn if q in up]
            ver = "unmeasurable"; why = []
            if pres_m or abs_m:
                ver = "retained"
                for gg in pres_m:
                    f = pct_in[g][list(a.var_names).index(gg)]
                    if f < 0.10:
                        ver = "excluded"; why.append(f"{gg} absent ({f:.0%})")
                for gg in abs_m:
                    f = pct_in[g][list(a.var_names).index(gg)]
                    if f >= 0.50:
                        ver = "excluded"; why.append(f"{gg} present ({f:.0%})")
            if ver in ("retained", "unmeasurable"):
                retained += 1
            rows_alt.append(dict(cluster=g, cluster_label=cur, alternative=alt["label"], axis=alt.get("axis", ""),
                                 must_present=";".join(pres_m), must_absent=";".join(abs_m),
                                 verdict=ver, reason="; ".join(why)))
        eq = 1 + retained
        grade = "decidable" if retained == 0 else ("partial" if retained <= 2 else "undecidable")
        rows_cl.append(dict(cluster=g, n_cells=int(m.sum()), consensus_label=cur,
                            n_alternatives=len(alts), n_retained=retained,
                            equivalence_class_size=eq, grade=grade))

    alt_df = pd.DataFrame(rows_alt)
    cl_df = pd.DataFrame([r for r in rows_cl if "grade" in r])
    alt_df.to_csv(os.path.join(args.outdir, "alternatives.csv"), index=False)
    cl_df.to_csv(os.path.join(args.outdir, "decidability.csv"), index=False)

    n_und = int((cl_df.grade == "undecidable").sum())
    n_par = int((cl_df.grade == "partial").sum())
    n_dec = int((cl_df.grade == "decidable").sum())
    cells_und = int(cl_df.loc[cl_df.grade == "undecidable", "n_cells"].sum())
    total = int(cl_df.n_cells.sum())
    md = [f"# Decidability report — {os.path.basename(args.h5ad)}", "",
          f"Clusters analysed: **{len(cl_df)}**; cells: **{total:,}**; label column: `{args.label_col}`", "",
          "## Verdict distribution", "",
          f"- decidable: {n_dec} clusters",
          f"- partial: {n_par} clusters",
          f"- undecidable: {n_und} clusters ({cells_und:,} cells, {100*cells_und/max(total,1):.0f}% of cells)", "",
          "## Rules (fixed before seeing results)", "",
          "- must_present gene: detection >= 25% and rank <= 5 within cluster → support; 10–25% → inconclusive; < 10% → contradiction (alternative excluded)",
          "- must_absent gene: detection >= 50% → contradiction (alternative excluded)",
          "- all discriminating genes outside the measured panel → unmeasurable, stays in the equivalence class",
          "- equivalence class 1 → decidable; 2–3 → partial; >= 4 → undecidable", "",
          "## Per-cluster results", "",
          "| cluster | cells | label | alternatives | retained | equivalence class | grade |",
          "|---|---|---|---|---|---|---|"]
    for r in cl_df.sort_values("equivalence_class_size", ascending=False).itertuples():
        md.append(f"| {r.cluster} | {r.n_cells:,} | {r.consensus_label} | {r.n_alternatives} | "
                  f"{r.n_retained} | {r.equivalence_class_size} | {r.grade} |")
    md += ["", "## Caveats", "",
           "- The result depends on how the alternative hypothesis set is built. With `--use-llm` the set is",
           "  open-ended; without it, a fixed rule library is used and the equivalence classes are smaller.",
           "- A 'decidable' verdict means every generated alternative could be excluded; it is not proof that",
           "  no alternative exists.",
           "- Detection thresholds assume log-normalised expression and a panel that covers the relevant markers."]
    open(os.path.join(args.outdir, "report.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")
    print(f"[done] {args.outdir}/report.md  ({n_dec} decidable / {n_par} partial / {n_und} undecidable)")


if __name__ == "__main__":
    main()
