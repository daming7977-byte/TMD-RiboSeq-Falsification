from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

FORMAL = Path(
    "tmd_anchor_preregister/bridge/"
    "FORMAL_TMD_ANCHOR_CANDIDATES_v1p1.tsv"
)

PAIR = Path(
    "crossmap_window_preregister/output/"
    "GSE300977_PAIRED_REAL_PSEUDO_CROSSMAP_MASK_AUDIT_v1.tsv"
)

CATALOG = Path(
    "tmd_spatial_secondary_scan_v1/output/"
    "SECONDARY_11CENTER_GEOMETRY_CATALOG_v1.tsv"
)

COUNTDIR = Path(
    "tmd_spatial_secondary_scan_v1/geometry_counts"
)

OUTDIR = Path(
    "tmd_spatial_topology_compartment_interaction_v1/output"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

CENTERS = [-30,-15,0,15,30,45,60,75,90,105,120]
CLASSES = ["FIRST","INTERNAL"]
LOCATIONS = ["ER","mitochondria","cytosol"]

LIBRARIES = [
    "B01","B02","B03","B04","B05","B06",
    "B07","B08","B09","B10","B11","B12",
]

CONTRASTS = [
    ("ER",1,"B01","B02"),
    ("ER",2,"B11","B12"),
    ("mitochondria",1,"B09","B10"),
    ("mitochondria",2,"B07","B08"),
    ("cytosol",1,"B05","B06"),
    ("cytosol",2,"B03","B04"),
]

INTERACTIONS = [
    ("MITO_MINUS_ER","mitochondria","ER"),
    ("CYTOSOL_MINUS_ER","cytosol","ER"),
    ("MITO_MINUS_CYTOSOL","mitochondria","cytosol"),
]

BOOT_N = 10000
BOOT_BATCH = 200
BOOT_SEED_BASE = 30097775

# ============================================================
# 1. Transcript topology
# ============================================================

f = pd.read_csv(FORMAL, sep="\t", dtype=str)

ft = (
    f[["target_transcript","pass_class"]]
    .dropna()
    .drop_duplicates()
)

nclass = (
    ft.groupby("target_transcript")["pass_class"]
    .nunique()
)

if int((nclass > 1).sum()) != 0:
    raise SystemExit("TOPOLOGY_CONFLICT_GATE_FAIL")

tx_class = (
    ft.drop_duplicates("target_transcript")
    .set_index("target_transcript")["pass_class"]
    .to_dict()
)

print("TOPOLOGY_TRANSCRIPTS =", len(tx_class))
print("TOPOLOGY_MAPPING_GATE = PASS")

# ============================================================
# 2. Frozen paired eligible units
# ============================================================

d = pd.read_csv(PAIR, sep="\t", dtype=str)

d["center"] = pd.to_numeric(d["center"]).astype(int)
d["real_crossmap_unmasked_nt"] = pd.to_numeric(
    d["real_crossmap_unmasked_nt"]
)
d["pseudo_min_unmasked_nt"] = pd.to_numeric(
    d["pseudo_min_unmasked_nt"]
)

d = d[
    (d["real_crossmap_unmasked_nt"] >= 1)
    &
    (d["pseudo_min_unmasked_nt"] >= 1)
].copy()

def splitcsv(x):
    return [z for z in str(x).split(",") if z != ""]

def topology(row):
    txs = splitcsv(row["target_transcripts"])
    classes = sorted({tx_class[x] for x in txs})

    if classes == ["MULTI_PASS"]:
        return "MULTI_PASS"
    if classes == ["SINGLE_PASS"]:
        return "SINGLE_PASS"
    return "MIXED_TOPOLOGY"

d["topology_class"] = d.apply(topology, axis=1)

m = d[
    d["topology_class"].eq("MULTI_PASS")
    &
    d["position_class_consensus"].isin(CLASSES)
].copy()

print("MULTIPASS_FIRST_INTERNAL_ROWS =", len(m))
print("MULTIPASS_UNIVERSE_GATE = PASS")

# ============================================================
# 3. Geometry catalog
# ============================================================

cat = pd.read_csv(CATALOG, sep="\t", dtype=str)

key_to_gid = {}

for _,r in cat.iterrows():
    key = (
        str(r["strand"]),
        str(r["geometry"]),
        str(r["mask_bits"]),
    )

    if key in key_to_gid:
        raise SystemExit("CATALOG_DUPLICATE_GATE_FAIL")

    key_to_gid[key] = str(r["geometry_id"])

print("GEOMETRY_CATALOG_GATE = PASS")

# ============================================================
# 4. Reuse completed counts
# ============================================================

density = {}

for lib in LIBRARIES:

    p = (
        COUNTDIR /
        f"{lib}.SECONDARY_GEOMETRY_DENSITIES_v1.tsv"
    )

    if not p.exists():
        raise SystemExit(f"MISSING_COUNT_CHECKPOINT:{lib}")

    x = pd.read_csv(p, sep="\t")

    density[lib] = dict(
        zip(
            x["geometry_id"].astype(str),
            x["density"].astype(float),
        )
    )

print("SECONDARY_COUNT_CHECKPOINTS = REUSED")

# ============================================================
# Helpers
# ============================================================

def splitsemi(x):
    return [z for z in str(x).split(";") if z != ""]

COLS = []

for lib in LIBRARIES:
    COLS.extend([
        f"{lib}_real",
        f"{lib}_pseudo",
    ])

COLIDX = {x:i for i,x in enumerate(COLS)}

def build_matrix(z):

    rows = []

    for _,r in z.iterrows():

        rkey = (
            str(r["strand"]),
            str(r["real_window_geometry"]),
            str(r["real_crossmap_mask_bits"]),
        )

        rgid = key_to_gid[rkey]

        pgeoms = splitsemi(r["pseudo_window_geometries"])
        pbits = splitsemi(r["pseudo_crossmap_mask_bits"])

        if len(pgeoms) != len(pbits):
            raise SystemExit("PSEUDO_CARDINALITY_GATE_FAIL")

        pgids = [
            key_to_gid[(str(r["strand"]),g,b)]
            for g,b in zip(pgeoms,pbits)
        ]

        row = []

        for lib in LIBRARIES:

            rd = float(density[lib][rgid])

            pdens = float(np.mean(
                [density[lib][gid] for gid in pgids]
            ))

            row.extend([rd,pdens])

        rows.append(row)

    X = np.asarray(rows, dtype=float)

    if (
        X.ndim != 2
        or X.shape[1] != 24
        or not np.isfinite(X).all()
        or (X < 0).any()
    ):
        raise SystemExit("MATRIX_GATE_FAIL")

    return X

def Dvec(M,pdlib,inlib):

    Rpd = M[:,COLIDX[f"{pdlib}_real"]]
    Rin = M[:,COLIDX[f"{inlib}_real"]]
    Ppd = M[:,COLIDX[f"{pdlib}_pseudo"]]
    Pin = M[:,COLIDX[f"{inlib}_pseudo"]]

    if (
        (Rpd <= 0).any()
        or (Rin <= 0).any()
        or (Ppd <= 0).any()
        or (Pin <= 0).any()
    ):
        raise SystemExit("ZERO_AGGREGATE_GATE_FAIL")

    return (
        np.log2(Rpd/Rin)
        -
        np.log2(Ppd/Pin)
    )

# ============================================================
# 5. Class matrices
# ============================================================

matrices = {}
account = []

for c in CENTERS:

    for cls in CLASSES:

        z = m[
            m["center"].eq(c)
            &
            m["position_class_consensus"].eq(cls)
        ].copy()

        if len(z) == 0:
            raise SystemExit(f"EMPTY_STRATUM:{c}:{cls}")

        X = build_matrix(z)

        matrices[(c,cls)] = X

        account.append({
            "center":c,
            "TMD_class":cls,
            "N_units":len(z),
        })

        print(
            f"CENTER={c:+d}",
            f"CLASS={cls}",
            f"N={len(z)}"
        )

pd.DataFrame(account).to_csv(
    OUTDIR /
    "TOPOLOGY_COMPARTMENT_INTERACTION_ACCOUNTING_v1.tsv",
    sep="\t",
    index=False
)

print("CLASS_MATRIX_GATE = PASS")

# ============================================================
# 6. Observed D and FI
# ============================================================

Dobs = {}

for c in CENTERS:

    for cls in CLASSES:

        mean = matrices[(c,cls)].mean(
            axis=0,
            keepdims=True
        )

        for loc,rep,pdlib,inlib in CONTRASTS:

            Dobs[(c,cls,loc,rep)] = float(
                Dvec(mean,pdlib,inlib)[0]
            )

observed_rows = []

for c in CENTERS:

    FI = {}

    for loc in LOCATIONS:

        f1 = Dobs[(c,"FIRST",loc,1)]
        f2 = Dobs[(c,"FIRST",loc,2)]

        i1 = Dobs[(c,"INTERNAL",loc,1)]
        i2 = Dobs[(c,"INTERNAL",loc,2)]

        FI[(loc,1)] = f1-i1
        FI[(loc,2)] = f2-i2

    for name,A,B in INTERACTIONS:

        r1 = FI[(A,1)] - FI[(B,1)]
        r2 = FI[(A,2)] - FI[(B,2)]

        observed_rows.append({
            "center":c,
            "interaction":name,
            "rep1":r1,
            "rep2":r2,
            "mean":(r1+r2)/2,
            "replicate_concordance":
                (
                    "SAME_SIGN"
                    if np.sign(r1)==np.sign(r2)
                    else "OPPOSITE_SIGN"
                ),
        })

observed = pd.DataFrame(observed_rows)

print("OBSERVED_INTERACTION_GATE = PASS")

# ============================================================
# 7. Bootstrap direct interaction
# ============================================================

boot_store = {
    (c,name):[]
    for c in CENTERS
    for name,_,_ in INTERACTIONS
}

for center_index,c in enumerate(CENTERS):

    XF = matrices[(c,"FIRST")]
    XI = matrices[(c,"INTERNAL")]

    NF = XF.shape[0]
    NI = XI.shape[0]

    rngF = np.random.default_rng(
        BOOT_SEED_BASE + 1000*center_index + 11
    )

    rngI = np.random.default_rng(
        BOOT_SEED_BASE + 1000*center_index + 29
    )

    pF = np.full(NF,1.0/NF)
    pI = np.full(NI,1.0/NI)

    done = 0

    while done < BOOT_N:

        b = min(BOOT_BATCH,BOOT_N-done)

        WF = rngF.multinomial(
            NF,
            pF,
            size=b
        )

        WI = rngI.multinomial(
            NI,
            pI,
            size=b
        )

        MF = (WF @ XF) / NF
        MI = (WI @ XI) / NI

        DF = {}
        DI = {}

        for loc,rep,pdlib,inlib in CONTRASTS:

            DF[(loc,rep)] = Dvec(
                MF,pdlib,inlib
            )

            DI[(loc,rep)] = Dvec(
                MI,pdlib,inlib
            )

        FI = {}

        for loc in LOCATIONS:

            FI[(loc,1)] = (
                DF[(loc,1)]
                -
                DI[(loc,1)]
            )

            FI[(loc,2)] = (
                DF[(loc,2)]
                -
                DI[(loc,2)]
            )

        for name,A,B in INTERACTIONS:

            r1 = FI[(A,1)] - FI[(B,1)]
            r2 = FI[(A,2)] - FI[(B,2)]

            boot_store[(c,name)].append(
                (r1+r2)/2
            )

        done += b

    print(
        "BOOTSTRAP_CENTER_COMPLETED =",
        f"{c:+d}"
    )

print("INTERACTION_BOOTSTRAP_GATE = PASS")

# ============================================================
# 8. Final table
# ============================================================

rows = []

for _,r in observed.iterrows():

    c = int(r["center"])
    name = r["interaction"]

    vals = np.concatenate(
        boot_store[(c,name)]
    )

    lo = float(
        np.percentile(vals,2.5)
    )

    hi = float(
        np.percentile(vals,97.5)
    )

    rows.append({
        **r.to_dict(),
        "CI_2.5":lo,
        "CI_97.5":hi,
        "CI_excludes_zero":
            (
                "YES"
                if lo > 0 or hi < 0
                else "NO"
            ),
    })

out = pd.DataFrame(rows)

out.to_csv(
    OUTDIR /
    "MULTIPASS_TOPOLOGY_COMPARTMENT_INTERACTIONS_v1.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 9. Reveal
# ============================================================

print()
print(
    "===================================================="
)
print(
    "DIRECT MULTI-PASS TOPOLOGY × COMPARTMENT INTERACTION"
)
print(
    "===================================================="
)

for name,_,_ in INTERACTIONS:

    print()
    print("INTERACTION =",name)

    z = out[
        out["interaction"].eq(name)
    ]

    for _,r in z.iterrows():

        print(
            f"CENTER={int(r['center']):+d}",
            f"REP1={r['rep1']:.12f}",
            f"REP2={r['rep2']:.12f}",
            f"MEAN={r['mean']:.12f}",
            f"CI95=[{r['CI_2.5']:.12f},"
            f"{r['CI_97.5']:.12f}]",
            f"CONCORDANCE={r['replicate_concordance']}",
            f"CI_EXCLUDES_0={r['CI_excludes_zero']}",
        )

print()
print(
    "TOPOLOGY_COMPARTMENT_INTERACTION_GATE = PASS"
)

# ============================================================
# 10. SHA
# ============================================================

for path in [
    OUTDIR /
    "TOPOLOGY_COMPARTMENT_INTERACTION_ACCOUNTING_v1.tsv",

    OUTDIR /
    "MULTIPASS_TOPOLOGY_COMPARTMENT_INTERACTIONS_v1.tsv",
]:

    sha = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()

    print(
        "SHA256",
        sha,
        path
    )
