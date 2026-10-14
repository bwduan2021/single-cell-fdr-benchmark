"""Apply a statistic and FDR rule to two supplied count matrices."""

import argparse
from pathlib import Path
import numpy as np
from methods import ALPHA, STAT_FUNC, run_one


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="NPZ file with cell-by-gene arrays XA and XB")
    parser.add_argument("--statistic", choices=tuple(STAT_FUNC), required=True)
    parser.add_argument("--rule", choices=("BH", "eBH", "DS", "TD"), required=True)
    parser.add_argument("--method-seed", type=int,
                        help="Recorded method seed; required for DS and TD")
    parser.add_argument("--output", type=Path,
                        help="Optional output NPZ file containing the rejection mask")
    args = parser.parse_args()
    if args.rule in ("DS", "TD") and args.method_seed is None:
        parser.error("DS and TD require --method-seed")
    if args.method_seed is not None and args.method_seed < 0:
        parser.error("--method-seed must be nonnegative")

    with np.load(args.input, allow_pickle=False) as data:
        if not {"XA", "XB"}.issubset(data.files):
            parser.error("the input file must contain XA and XB")
        XA = np.asarray(data["XA"], dtype=float)
        XB = np.asarray(data["XB"], dtype=float)
    minimum = 4 if args.rule == "DS" else 2
    if XA.ndim != 2 or XB.ndim != 2:
        parser.error("XA and XB must be two-dimensional arrays")
    if XA.shape[1] == 0 or XA.shape[1] != XB.shape[1]:
        parser.error("XA and XB must have the same nonempty gene columns")
    if min(XA.shape[0], XB.shape[0]) < minimum:
        parser.error(f"each group requires at least {minimum} observations")
    if any(not np.isfinite(x).all() or (x < 0).any() for x in (XA, XB)):
        parser.error("counts must be finite and nonnegative")

    rng = np.random.default_rng(args.method_seed)
    rejected = run_one(XA, XB, args.statistic, args.rule, rng)
    print(f"Statistic: {args.statistic}; rule: {args.rule}; alpha: {ALPHA}")
    print(f"Discoveries: {int(rejected.sum())}")
    print("Selected gene indices (zero-based):", np.flatnonzero(rejected).tolist())
    if args.output is not None:
        with args.output.open("wb") as handle:
            np.savez_compressed(handle, rejected=rejected)


if __name__ == "__main__":
    main()
