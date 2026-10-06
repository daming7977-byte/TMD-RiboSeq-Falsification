"""
Stage0C IO utilities.
"""

from pathlib import Path


def read_features(path):
    """
    Read human_features.tsv.

    Returns
    -------
    dict:
        transcript ->
        {
            utr5,
            utr3,
            cds,
            gc_utr5,
            gc_utr3,
            gc_cds
        }
    """

    path = Path(path)

    features = {}

    with path.open() as fh:

        next(fh)

        for line in fh:

            f = line.rstrip("\n").split("\t")

            if len(f) < 7:
                continue

            features[f[0]] = {

                "utr5": int(f[1]),

                "utr3": int(f[2]),

                "cds": int(f[3]),

                "gc_utr5": float(f[4]),

                "gc_utr3": float(f[5]),

                "gc_cds": float(f[6]),

            }

    return features


def cds_start(features, transcript):

    if transcript == "DENGUE":
        return 96

    return features[transcript]["utr5"]


def write_tsv(rows, outfile):

    from csv import DictWriter

    rows = list(rows)

    if len(rows) == 0:
        raise ValueError("No rows.")

    outfile = Path(outfile)

    outfile.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with outfile.open("w", newline="") as fh:

        writer = DictWriter(
            fh,
            fieldnames=rows[0].keys(),
            delimiter="\t",
        )

        writer.writeheader()

        writer.writerows(rows)
