#!/usr/bin/env python
"""Auditoría de ChEMBL: ¿qué dianas tienen IC50 pareados WT/mutante?"""
from pathlib import Path
import pandas as pd
from chembl_webresource_client.new_client import new_client

# Los IDs se verifican al descargar (se imprime el nombre de cada diana)
CANDIDATES = {
    "EGFR": "CHEMBL203",
    "ABL1": "CHEMBL1862",
    "ALK": "CHEMBL4247",
    "KIT": "CHEMBL1936",
    "BRAF": "CHEMBL5145",
}
RAW = Path("data/raw")
TAB = Path("results/tables")
RAW.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

FIELDS = [
    "molecule_chembl_id", "canonical_smiles", "standard_type",
    "standard_relation", "standard_value", "standard_units",
    "pchembl_value", "assay_chembl_id", "assay_type", "assay_description",
    "document_chembl_id", "assay_variant_accession",
    "assay_variant_mutation", "data_validity_comment",
]
WT = "WT_o_sin_anotar"


def download(gene, chembl_id):
    path = RAW / f"chembl_{gene}_IC50.csv"
    if path.exists():
        return pd.read_csv(path)
    t = new_client.target.get(chembl_id)
    print(f"{gene}: {chembl_id} -> {t['pref_name']} ({t['organism']})")
    res = new_client.activity.filter(
        target_chembl_id=chembl_id, standard_type="IC50").only(FIELDS)
    df = pd.DataFrame(list(res))
    df.to_csv(path, index=False)
    return df


def normalize_mutation(x):
    if pd.isna(x) or not str(x).strip():
        return WT
    parts = str(x).replace(";", ",").replace("/", ",").split(",")
    return "/".join(sorted(p.strip() for p in parts if p.strip()))


def clean(df):
    df = df[df["standard_relation"] == "="]
    df = df[df["standard_units"] == "nM"]
    df = df[df["data_validity_comment"].isna()]
    df = df[df["pchembl_value"].notna()].copy()
    df["pchembl_value"] = df["pchembl_value"].astype(float)
    df["variant"] = df["assay_variant_mutation"].map(normalize_mutation)
    return df[df["assay_type"] == "B"]   # solo ensayos bioquímicos


def paired(df):
    g = (df.groupby(["molecule_chembl_id", "document_chembl_id", "variant"])
           ["pchembl_value"].median().reset_index())
    wt = (g[g.variant == WT].drop(columns="variant")
            .rename(columns={"pchembl_value": "pchembl_wt"}))
    mut = (g[g.variant != WT]
           .rename(columns={"pchembl_value": "pchembl_mut"}))
    p = mut.merge(wt, on=["molecule_chembl_id", "document_chembl_id"])
    p["delta_pchembl"] = p["pchembl_mut"] - p["pchembl_wt"]
    return p


rows, all_pairs = [], []
for gene, cid in CANDIDATES.items():
    df = clean(download(gene, cid))
    p = paired(df)
    p["gene"] = gene
    all_pairs.append(p)
    rows.append({
        "gene": gene,
        "n_measurements_B": len(df),
        "n_paired_points": len(p),
        "n_compounds": p["molecule_chembl_id"].nunique(),
        "n_documents": p["document_chembl_id"].nunique(),
        "n_mutations": p["variant"].nunique(),
        "n_|delta|>=1": int((p["delta_pchembl"].abs() >= 1).sum()),
        "sd_delta": round(p["delta_pchembl"].std(), 2),
    })

summary = pd.DataFrame(rows)
summary.to_csv(TAB / "target_audit_summary.csv", index=False)
pd.concat(all_pairs).to_csv(TAB / "target_audit_pairs.csv", index=False)

print("\n=== RESUMEN ===")
print(summary.to_string(index=False))
for p in all_pairs:
    print(f"\n--- {p['gene'].iloc[0] if len(p) else '?'}: top mutaciones ---")
    print(p["variant"].value_counts().head(8).to_string())

