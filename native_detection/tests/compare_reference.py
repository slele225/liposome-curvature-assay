#!/usr/bin/env python3
"""Compare cme_detect --dump-dir output against the MATLAB golden dump
(reference_tools/dump_reference.m), stage by stage, and the end-to-end tables
(detection_cpp.tsv vs detection_matlab.tsv).

Usage:
  python compare_reference.py --cpp <cpp dump dir> --matlab <matlab dump dir> [--atol 1e-8 --rtol 1e-6]
  python compare_reference.py --e2e-cpp <detection_cpp.tsv> --e2e-matlab <detection_matlab.tsv>
"""
import argparse
import glob
import math
import os
import sys

import numpy as np

DISCRETE_COLS = {"hval_Ar", "hval_AD", "isPSF", "x_init", "y_init", "mask_Ar", "maskN", "keepNonNaN", "keepFinal", "useLoc", "idx", "frame", "index", "iters", "component", "k", "converged", "image", "count", "channel"}


def read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
    if not lines:
        return [], {}
    header = lines[0].split("\t")
    cols = {h: [] for h in header}
    for l in lines[1:]:
        if not l.strip():
            continue
        parts = l.split("\t")
        for h, v in zip(header, parts):
            cols[h].append(v)
    out = {}
    for h in header:
        try:
            out[h] = np.array([float(v) for v in cols[h]], dtype=float)
        except ValueError:
            out[h] = np.array(cols[h], dtype=object)
    return header, out


def read_bin(path):
    with open(path + ".hdr") as f:
        cls, ny, nx, _ = f.read().split()
    ny, nx = int(ny), int(nx)
    dt = np.float64 if cls == "float64" else np.uint8
    a = np.fromfile(path, dtype=dt)
    return a.reshape((nx, ny)).T  # column-major


class Report:
    def __init__(self, atol, rtol):
        self.atol, self.rtol = atol, rtol
        self.fail = 0
        self.pass_ = 0
        self.worst = {}
        self.msgs = []

    def cmp_array(self, name, a, b, discrete=False):
        a = np.asarray(a); b = np.asarray(b)
        if a.shape != b.shape:
            self.fail += 1
            self.msgs.append(f"FAIL {name}: shape {a.shape} vs {b.shape}")
            return False
        if a.dtype == object or b.dtype == object:
            ok = np.all(a == b)
        else:
            nan_a, nan_b = np.isnan(a), np.isnan(b)
            if not np.array_equal(nan_a, nan_b):
                self.fail += 1
                self.msgs.append(f"FAIL {name}: NaN pattern differs ({nan_a.sum()} vs {nan_b.sum()})")
                return False
            m = ~nan_a
            if discrete:
                ok = np.array_equal(a[m], b[m])
                if not ok:
                    self.msgs.append(f"FAIL {name}: {int(np.sum(a[m] != b[m]))} discrete mismatches of {int(m.sum())}")
            else:
                diff = np.abs(a[m] - b[m])
                tol = self.atol + self.rtol * np.abs(b[m])
                bad = diff > tol
                rel = diff / np.maximum(np.abs(b[m]), 1e-300)
                if m.sum():
                    key = name.split("/")[-1]
                    self.worst[key] = max(self.worst.get(key, 0.0), float(np.max(diff)), )
                ok = not np.any(bad)
                if not ok:
                    i = int(np.argmax(diff - tol))
                    self.msgs.append(f"FAIL {name}: {int(bad.sum())}/{int(m.sum())} beyond tol; worst abs {diff[i]:.3g} rel {rel[i]:.3g} (cpp {a[m][i]:.17g} vs matlab {b[m][i]:.17g})")
        if ok:
            self.pass_ += 1
        else:
            self.fail += 1
        return ok

    def cmp_tsv(self, name, pa, pb, cols=None):
        if not os.path.exists(pa) or not os.path.exists(pb):
            self.fail += 1
            self.msgs.append(f"FAIL {name}: missing file ({os.path.exists(pa)}, {os.path.exists(pb)})")
            return
        ha, a = read_tsv(pa)
        hb, b = read_tsv(pb)
        common = [c for c in hb if c in a] if cols is None else cols
        if not common:
            return
        na = len(next(iter(a.values()))) if a else 0
        nb = len(next(iter(b.values()))) if b else 0
        if na != nb:
            self.fail += 1
            self.msgs.append(f"FAIL {name}: row count {na} vs {nb}")
            return
        for c in common:
            base = c.rsplit("_", 1)[0] if c[-1].isdigit() and "_" in c else c
            discrete = (c in DISCRETE_COLS) or (base in DISCRETE_COLS)
            self.cmp_array(f"{name}/{c}", a[c], b[c], discrete=discrete)


def compare_dumps(cpp, mat, rep):
    # sigma
    rep.cmp_tsv("sigma.tsv", os.path.join(cpp, "sigma.tsv"), os.path.join(mat, "sigma.tsv"))
    for f in sorted(glob.glob(os.path.join(mat, "sigma_ch*_svect.tsv"))):
        rep.cmp_tsv(os.path.basename(f), os.path.join(cpp, os.path.basename(f)), f)
    for f in sorted(glob.glob(os.path.join(mat, "sigma_ch*_svect_per_image.tsv"))):
        rep.cmp_tsv(os.path.basename(f), os.path.join(cpp, os.path.basename(f)), f)
    for f in sorted(glob.glob(os.path.join(mat, "sigma_ch*_gmm.tsv"))):
        rep.cmp_tsv(os.path.basename(f), os.path.join(cpp, os.path.basename(f)), f)
    for f in sorted(glob.glob(os.path.join(mat, "sigma_ch*_refit_img*.tsv"))):
        rep.cmp_tsv(os.path.basename(f), os.path.join(cpp, os.path.basename(f)), f)
    # movies
    for md in sorted(glob.glob(os.path.join(mat, "movie*"))):
        name = os.path.basename(md)
        cd = os.path.join(cpp, name)
        if not os.path.isdir(cd):
            rep.fail += 1
            rep.msgs.append(f"FAIL {name}: missing in cpp dump")
            continue
        for f in sorted(glob.glob(os.path.join(md, "frame*_*.bin"))):
            b = read_bin(f)
            cf = os.path.join(cd, os.path.basename(f))
            if not os.path.exists(cf):
                rep.fail += 1; rep.msgs.append(f"FAIL {name}/{os.path.basename(f)}: missing")
                continue
            a = read_bin(cf)
            discrete = "mask" in os.path.basename(f)
            rep.cmp_array(f"{name}/{os.path.basename(f)}", a.astype(float), b.astype(float), discrete=discrete)
        for f in sorted(glob.glob(os.path.join(md, "frame*_*.tsv"))):
            rep.cmp_tsv(f"{name}/{os.path.basename(f)}", os.path.join(cd, os.path.basename(f)), f)
        for f in sorted(glob.glob(os.path.join(md, "frame*_logThreshold.txt"))):
            cf = os.path.join(cd, os.path.basename(f))
            if os.path.exists(cf):
                a = float(open(cf).read()); b = float(open(f).read())
                rep.cmp_array(f"{name}/{os.path.basename(f)}", np.array([a]), np.array([b]))


def compare_e2e(pc, pm, rep):
    hc, a = read_tsv(pc)
    hm, b = read_tsv(pm)
    na = len(a["frame"]) if "frame" in a else 0
    nb = len(b["frame"]) if "frame" in b else 0
    print(f"rows: cpp {na}, matlab {nb}")
    if na != nb:
        rep.fail += 1
        rep.msgs.append(f"FAIL e2e: row count {na} vs {nb}")
        return
    for c in hm:
        if c in ("movie", "source_image"):
            continue
        if c not in a:
            continue
        base = c.rsplit("_", 1)[0]
        discrete = c in DISCRETE_COLS or base in DISCRETE_COLS
        rep.cmp_array(f"e2e/{c}", a[c], b[c], discrete=discrete)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpp")
    ap.add_argument("--matlab")
    ap.add_argument("--e2e-cpp")
    ap.add_argument("--e2e-matlab")
    ap.add_argument("--atol", type=float, default=1e-8)
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--max-msgs", type=int, default=60)
    args = ap.parse_args()
    rep = Report(args.atol, args.rtol)
    if args.cpp and args.matlab:
        compare_dumps(args.cpp, args.matlab, rep)
    if args.e2e_cpp and args.e2e_matlab:
        compare_e2e(args.e2e_cpp, args.e2e_matlab, rep)
    for m in rep.msgs[: args.max_msgs]:
        print(m)
    if len(rep.msgs) > args.max_msgs:
        print(f"... {len(rep.msgs) - args.max_msgs} more")
    print(f"{rep.pass_} comparisons passed, {rep.fail} failed (atol {args.atol}, rtol {args.rtol})")
    worst = sorted(rep.worst.items(), key=lambda kv: -kv[1])[:15]
    print("largest absolute deviations per column:")
    for k, v in worst:
        print(f"  {k:24s} {v:.3g}")
    sys.exit(0 if rep.fail == 0 else 1)


if __name__ == "__main__":
    main()
