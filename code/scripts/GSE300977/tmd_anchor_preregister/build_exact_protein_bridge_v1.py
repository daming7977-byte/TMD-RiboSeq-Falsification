from pathlib import Path
from collections import defaultdict
import hashlib
import re
import subprocess
import pandas as pd

ANCHOR = Path(
    "REPRO_WORKSPACE/GSE133111_stage0C/reference/tmd/"
    "primary_tmd_anchors.tsv"
)

GTF = Path(
    "mapping_reference_blinded/"
    "Homo_sapiens.GRCh38.104.gtf"
)

GENOME = Path(
    "mapping_reference_blinded/"
    "Homo_sapiens.GRCh38.release104.primary_assembly.fa"
)

OUTDIR = Path("tmd_anchor_preregister/bridge")
AUDIT  = Path("tmd_anchor_preregister/audit")

OUTDIR.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

EXPECTED_ANCHOR_SHA = (
    "7257a07c51158ea4235c1c45a01ada899fcd78bcaddba68518bbb9a2146a92df"
)

EXPECTED_GTF_SHA = (
    "f3de2bccc02877709272a0ce66af59a2922d407a4076541bc5a65c9d3518581e"
)

EXPECTED_GENOME_SHA = (
    "78777b0886e8dfa5e14e4957fbbaa53736fcbaa5668d59e09b6b7945fca93d8c"
)

def sha256_file(path, block=1024*1024):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            b = fh.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

print("=== SHA GATES ===")

anchor_sha = sha256_file(ANCHOR)
gtf_sha = sha256_file(GTF)
genome_sha = sha256_file(GENOME)

print("ANCHOR_SHA =", anchor_sha)
print("GTF_SHA    =", gtf_sha)
print("GENOME_SHA =", genome_sha)

if anchor_sha != EXPECTED_ANCHOR_SHA:
    raise SystemExit("ANCHOR_SHA_GATE = FAIL")

if gtf_sha != EXPECTED_GTF_SHA:
    raise SystemExit("GTF_SHA_GATE = FAIL")

if genome_sha != EXPECTED_GENOME_SHA:
    raise SystemExit("GENOME_SHA_GATE = FAIL")

print("REFERENCE_SHA_GATES = PASS")

# ----------------------------------------------------------
# FASTA indexed random access
# ----------------------------------------------------------

FAI = Path(str(GENOME) + ".fai")

if not FAI.exists():
    subprocess.run(
        ["samtools", "faidx", str(GENOME)],
        check=True
    )

fai = {}

with FAI.open() as fh:
    for line in fh:
        chrom, length, offset, line_bases, line_width = line.rstrip().split("\t")[:5]
        fai[chrom] = tuple(
            map(int, (length, offset, line_bases, line_width))
        )

genome_fh = GENOME.open("rb")

def fetch(chrom, start0, end0):
    length, offset, line_bases, line_width = fai[chrom]

    if not (0 <= start0 <= end0 <= length):
        raise ValueError(
            f"OUT_OF_RANGE {chrom}:{start0}-{end0}"
        )

    def bytepos(pos):
        return (
            offset
            + (pos // line_bases) * line_width
            + (pos % line_bases)
        )

    a = bytepos(start0)
    b = bytepos(end0)

    genome_fh.seek(a)
    raw = genome_fh.read(b-a + 8)

    seq = raw.replace(b"\n", b"").replace(b"\r", b"")
    return seq[:end0-start0].decode().upper()

COMP = str.maketrans(
    "ACGTNacgtn",
    "TGCANtgcan"
)

def revcomp(s):
    return s.translate(COMP)[::-1]

# ----------------------------------------------------------
# Standard translation
# ----------------------------------------------------------

CODON = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def translate(cds):
    aa = []
    upto = len(cds) - len(cds) % 3

    for i in range(0, upto, 3):
        aa.append(CODON.get(cds[i:i+3], 'X'))

    p = ''.join(aa)

    # Ensembl CDS normally excludes stop_codon, but make handling explicit
    if p.endswith("*"):
        p = p[:-1]

    return p

# ----------------------------------------------------------
# Parse Ensembl104 CDS by transcript
# ----------------------------------------------------------

attr_re = re.compile(r'(\S+) "([^"]*)"')

cds_by_tx = defaultdict(list)
tx_meta = {}

with GTF.open() as fh:
    for line in fh:
        if line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")

        if len(f) != 9 or f[2] != "CDS":
            continue

        chrom = f[0]
        start1 = int(f[3])
        end1 = int(f[4])
        strand = f[6]
        attrs = dict(attr_re.findall(f[8]))

        tx = attrs.get("transcript_id")

        if tx is None:
            continue

        cds_by_tx[tx].append(
            (chrom, start1-1, end1, strand)
        )

        tx_meta[tx] = {
            "gene_id": attrs.get("gene_id", ""),
            "gene_name": attrs.get("gene_name", ""),
            "transcript_biotype": attrs.get(
                "transcript_biotype", ""
            ),
        }

print()
print("TARGET_TRANSCRIPTS_WITH_CDS =", len(cds_by_tx))

# ----------------------------------------------------------
# Build translated protein SHA table
# ----------------------------------------------------------

target_rows = []
sha_to_targets = defaultdict(list)

bad_chrom = 0
non_multiple_three = 0
internal_stop = 0
empty_protein = 0

for tx, segs in cds_by_tx.items():

    strand = segs[0][3]

    if strand == "+":
        segs = sorted(segs, key=lambda x: (x[1], x[2]))
    else:
        segs = sorted(
            segs,
            key=lambda x: (x[1], x[2]),
            reverse=True
        )

    pieces = []

    failed = False

    for chrom, start0, end0, st in segs:
        if chrom not in fai:
            failed = True
            bad_chrom += 1
            break

        s = fetch(chrom, start0, end0)

        if strand == "-":
            s = revcomp(s)

        pieces.append(s)

    if failed:
        continue

    cds = ''.join(pieces)

    if len(cds) % 3 != 0:
        non_multiple_three += 1

    protein = translate(cds)

    if not protein:
        empty_protein += 1
        continue

    if "*" in protein:
        internal_stop += 1

    protein_sha = hashlib.sha256(
        protein.encode()
    ).hexdigest()

    meta = tx_meta[tx]

    target_rows.append({
        "target_transcript": tx,
        "gene_id": meta["gene_id"],
        "gene_name": meta["gene_name"],
        "transcript_biotype": meta["transcript_biotype"],
        "cds_length": len(cds),
        "protein_length": len(protein),
        "protein_sha256": protein_sha,
        "internal_stop": "YES" if "*" in protein else "NO",
    })

    sha_to_targets[protein_sha].append(tx)

target = pd.DataFrame(target_rows)

target.to_csv(
    OUTDIR / "ENSEMBL104_TRANSLATED_PROTEIN_SHA_v1.tsv",
    sep="\t",
    index=False
)

print("TARGET_PROTEINS_BUILT =", len(target))
print("BAD_CHROM =", bad_chrom)
print("CDS_LENGTH_NOT_MULTIPLE_3 =", non_multiple_three)
print("PROTEIN_WITH_INTERNAL_STOP =", internal_stop)
print("EMPTY_PROTEIN =", empty_protein)

# ----------------------------------------------------------
# Read frozen source TMD universe
# ----------------------------------------------------------

a = pd.read_csv(
    ANCHOR,
    sep="\t",
    dtype=str
)

required = {
    "transcript",
    "protein_length",
    "protein_sha256",
    "TMD_index",
    "n_primary_TMDs",
    "TMD_position_class",
    "aa_start0",
    "aa_end0_exclusive",
}

missing = required - set(a.columns)

if missing:
    raise SystemExit(
        f"ANCHOR_REQUIRED_COLUMNS_FAIL: {sorted(missing)}"
    )

source = (
    a[
        [
            "transcript",
            "protein_length",
            "protein_sha256",
            "n_primary_TMDs",
        ]
    ]
    .drop_duplicates()
)

if source["transcript"].duplicated().any():
    raise SystemExit(
        "SOURCE_TRANSCRIPT_PROTEIN_IDENTITY_GATE = FAIL"
    )

print()
print("SOURCE_TMD_TRANSCRIPTS =", len(source))
print("SOURCE_PRIMARY_TMD_ROWS =", len(a))

# ----------------------------------------------------------
# Exact bridge classification
# ----------------------------------------------------------

bridge = []

for _, r in source.iterrows():

    src_tx = r["transcript"]
    sha = r["protein_sha256"]

    candidates = sorted(
        sha_to_targets.get(sha, [])
    )

    same = src_tx in candidates

    if same:
        cls = "SAME_ENST_EXACT"
        target_tx = src_tx

    elif len(candidates) == 1:
        cls = "UNIQUE_EXACT_OTHER_ENST"
        target_tx = candidates[0]

    elif len(candidates) > 1:
        cls = "MULTIPLE_EXACT"
        target_tx = ""

    else:
        cls = "NO_EXACT"
        target_tx = ""

    bridge.append({
        "source_transcript": src_tx,
        "source_protein_length": int(r["protein_length"]),
        "source_protein_sha256": sha,
        "source_n_primary_TMDs": int(r["n_primary_TMDs"]),
        "match_class": cls,
        "target_transcript": target_tx,
        "n_exact_target_candidates": len(candidates),
        "exact_target_candidates": ",".join(candidates),
    })

bridge = pd.DataFrame(bridge)

bridge.to_csv(
    OUTDIR / "APPRIS19192_TO_ENSEMBL104_EXACT_PROTEIN_BRIDGE_v1.tsv",
    sep="\t",
    index=False
)

summary = (
    bridge["match_class"]
    .value_counts()
    .rename_axis("match_class")
    .reset_index(name="transcripts")
)

summary.to_csv(
    OUTDIR / "EXACT_PROTEIN_BRIDGE_SUMMARY_v1.tsv",
    sep="\t",
    index=False
)

print()
print("=== BRIDGE SUMMARY ===")
print(summary.to_string(index=False))

admissible = bridge["match_class"].isin([
    "SAME_ENST_EXACT",
    "UNIQUE_EXACT_OTHER_ENST",
])

print()
print(
    "FORMALLY_ADMISSIBLE_TRANSCRIPTS =",
    int(admissible.sum())
)

print(
    "FORMALLY_EXCLUDED_TRANSCRIPTS =",
    int((~admissible).sum())
)

print("EXACT_PROTEIN_BRIDGE_EXECUTION = PASS")

genome_fh.close()
