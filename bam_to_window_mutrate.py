#!/usr/bin/env python3
"""
Run the BAM -> window-level mutation-rate pipeline on a single BAM or every
BAM in a folder.

For each input BAM <name>.bam the pipeline produces:
  <out>/<name>.counts.txt.gz                                 (bam-readcount | pigz output;
                                                              removed once mutrate is built)
  <out>/<name>.mutrate.txt.bgz                               (sorted bgzip of get_mutrate.v4.1.py
                                                              output; raw .mutrate.txt.gz is
                                                              produced as an intermediate and
                                                              removed after bgzip+index)
  <out>/<name>.mutrate.txt.bgz.tbi                           (tabix index: col 1 = gene,
                                                              col 2/3 = start/end; mirrors
                                                              bgzip_tabix_mutrate.sh)
  <out>/<name>.win{WIN_LEN}_cov{COV_THRESHOLD}.window.csv.gz (mutrate2window.py output;
                                                              filename encodes window
                                                              params so re-runs with
                                                              different settings don't
                                                              overwrite each other)

Each step is skipped if its output already exists.

Usage:
  python bam_to_window_mutrate.py \\
      --input <BAM_or_FOLDER> \\
      --ref-fasta <ref.fa> \\
      --output-dir <outdir> \\
      [--regions-bed <regions.bed>] \\
      [--win-len 10] [--step 10] \\
      [--coverage-threshold 200] [--mutrate-threshold 0.2] \\
      [--gene-readcounts featureCounts.txt] \\
      [--threads 4] [--bam-glob '*.bam'] [--ref-symbol]
"""

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from glob import glob
from pathlib import Path


# Default location of the two helper scripts (override with --shapemap-util-dir).
DEFAULT_SHAPEMAP_UTIL_DIR = "/home/users/astar/gis/linxy/code/shapemap_util"


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True,
                   help="Input BAM file OR a folder containing BAM files")
    p.add_argument("--ref-fasta", required=True,
                   help="Reference FASTA used by bam-readcount (must be indexed: .fai)")
    p.add_argument("--regions-bed", default=None,
                   help="Optional BED of regions to count (passed to bam-readcount -l). "
                        "If omitted, bam-readcount runs without -l (slower; counts every "
                        "covered position).")
    p.add_argument("--output-dir", required=True,
                   help="Directory to write counts/mutrate/window outputs")

    # Pipeline tunables (defaults match run_cellsnp_batched.sh / mutrate2window.py)
    p.add_argument("--win-len", type=int, default=10,
                   help="Window length for mutrate2window.py -w (default: 10)")
    p.add_argument("--step", type=int, default=10,
                   help="Step size for mutrate2window.py -s (default: 10)")
    p.add_argument("--coverage-threshold", type=int, default=200,
                   help="Per-window coverage threshold (default: 200)")
    p.add_argument("--mutrate-threshold", type=float, default=0.2,
                   help="Per-position mutrate threshold (default: 0.2)")
    p.add_argument("--gene-readcounts", default=None,
                   help="Optional featureCounts file for normalized coverage in get_mutrate")
    p.add_argument("--ref-symbol", action="store_true",
                   help="Pass --ref_symbol to get_mutrate (matched nt labeled as '=')")
    p.add_argument("--get-mutrate-w", type=int, default=0,
                   help="Per-gene win threshold passed to get_mutrate -w (default: 0)")

    # Orchestration
    p.add_argument("--threads", type=int, default=4,
                   help="Number of BAMs to process in parallel (default: 4)")
    p.add_argument("--bam-glob", default="*.bam",
                   help="Glob used when --input is a folder (default: '*.bam')")
    p.add_argument("--shapemap-util-dir", default=DEFAULT_SHAPEMAP_UTIL_DIR,
                   help=f"Directory holding get_mutrate.v4.1.py and mutrate2window.py "
                        f"(default: {DEFAULT_SHAPEMAP_UTIL_DIR})")
    p.add_argument("--python", default=sys.executable,
                   help="Python interpreter to use for the helper scripts "
                        "(default: this interpreter)")
    return p.parse_args()


def find_bams(input_path, bam_glob):
    """Resolve --input into a list of absolute BAM paths."""
    p = Path(input_path)
    if not p.exists():
        sys.exit(f"ERROR: --input not found: {input_path}")
    if p.is_file():
        if p.suffix != ".bam":
            sys.exit(f"ERROR: --input is a file but not a .bam: {input_path}")
        return [p.resolve()]
    if p.is_dir():
        bams = sorted(Path(b).resolve() for b in glob(str(p / bam_glob)))
        if not bams:
            sys.exit(f"ERROR: no files matching '{bam_glob}' in {input_path}")
        return bams
    sys.exit(f"ERROR: --input is neither a file nor a directory: {input_path}")


def check_dependencies(util_dir):
    """Validate external tools and helper scripts are present and reachable."""
    if shutil.which("bam-readcount") is None:
        sys.exit("ERROR: bam-readcount not found on PATH. "
                 "Activate the conda env that provides it before running.")
    if shutil.which("samtools") is None:
        sys.exit("ERROR: samtools not found on PATH (needed to auto-index BAMs).")
    if shutil.which("pigz") is None:
        sys.exit("ERROR: pigz not found on PATH (needed to compress bam-readcount output).")
    if shutil.which("bgzip") is None:
        sys.exit("ERROR: bgzip not found on PATH (needed to bgzip the mutrate output).")
    if shutil.which("tabix") is None:
        sys.exit("ERROR: tabix not found on PATH (needed to index the mutrate output).")

    util = Path(util_dir)
    get_mutrate = util / "get_mutrate.v4.1.py"
    mutrate2win = util / "mutrate2window.py"
    for s in (get_mutrate, mutrate2win):
        if not s.is_file():
            sys.exit(f"ERROR: helper script not found: {s}")
    return get_mutrate, mutrate2win


def ensure_bam_indexed(bam_path):
    """bam-readcount needs an indexed BAM. Auto-index if missing."""
    bai = bam_path.with_suffix(bam_path.suffix + ".bai")
    alt_bai = bam_path.with_suffix(".bai")
    if bai.exists() or alt_bai.exists():
        return
    print(f"[{bam_path.name}] no .bai found, running samtools index...", flush=True)
    subprocess.run(["samtools", "index", str(bam_path)], check=True)


def sort_bgzip_index(mutrate_raw, mutrate_bgz, log_prefix):
    """Sort the raw .mutrate.txt.gz by gene/start, bgzip, and tabix-index it.

    Mirrors bgzip_tabix_mutrate.sh: zcat | LC_ALL=C sort -k1,1 -k2,2n | bgzip,
    then `tabix -s 1 -b 2 -e 3` on the resulting .bgz. Writes through a .tmp
    sibling so a partial bgz never appears at the final path.
    """
    tmp_bgz = Path(str(mutrate_bgz) + ".tmp")
    if tmp_bgz.exists():
        tmp_bgz.unlink()

    print(f"{log_prefix} Step 2b sort | bgzip -> {mutrate_bgz.name}", flush=True)

    env = os.environ.copy()
    env["LC_ALL"] = "C"

    zcat = subprocess.Popen(
        ["zcat", str(mutrate_raw)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    sort_proc = subprocess.Popen(
        ["sort", "-k1,1", "-k2,2n", "-S", "2G"],
        stdin=zcat.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, env=env,
    )
    zcat.stdout.close()
    with open(tmp_bgz, "wb") as out_fh:
        bgzip = subprocess.Popen(
            ["bgzip", "-c"], stdin=sort_proc.stdout, stdout=out_fh, stderr=subprocess.DEVNULL,
        )
        sort_proc.stdout.close()
        bgzip_rc = bgzip.wait()
    sort_rc = sort_proc.wait()
    zcat_rc = zcat.wait()

    for name, rc in (("zcat", zcat_rc), ("sort", sort_rc), ("bgzip", bgzip_rc)):
        if rc != 0:
            if tmp_bgz.exists():
                try:
                    tmp_bgz.unlink()
                except OSError:
                    pass
            raise RuntimeError(f"{name} failed (exit {rc})")

    os.replace(tmp_bgz, mutrate_bgz)
    run_step(
        "tabix index",
        ["tabix", "-f", "-s", "1", "-b", "2", "-e", "3", str(mutrate_bgz)],
        log_prefix,
    )


def run_step(label, cmd, log_prefix):
    """Run one pipeline step; subprocess stdout/stderr are discarded."""
    print(f"{log_prefix} {label}: {' '.join(str(c) for c in cmd)}", flush=True)
    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed (exit {result.returncode})")


def process_bam(bam_path, args, get_mutrate, mutrate2win):
    """Run the 3-step pipeline on a single BAM. Returns (bam, ok, error_str)."""
    base = bam_path.stem  # filename without .bam
    out = Path(args.output_dir)
    counts_path = out / f"{base}.counts.txt.gz"
    # get_mutrate.v4.1.py writes a plain .gz; we then sort + bgzip it so single
    # genes can be fetched via tabix. The .gz is removed once the .bgz/.tbi
    # pair is in place.
    mutrate_raw_path = out / f"{base}.mutrate.txt.gz"
    mutrate_path = out / f"{base}.mutrate.txt.bgz"
    mutrate_tbi_path = out / f"{base}.mutrate.txt.bgz.tbi"
    # Encode the windowing params in the Step-3 output filename so a re-run
    # with different --win-len / --coverage-threshold doesn't overwrite earlier
    # results, and Step-1/Step-2 outputs are still reused across param sweeps.
    window_path = out / f"{base}.win{args.win_len}_cov{args.coverage_threshold}.window.csv.gz"
    log_prefix = f"[{base}]"

    try:
        ensure_bam_indexed(bam_path)

        mutrate_ready = (
            mutrate_path.exists() and mutrate_path.stat().st_size > 0
            and mutrate_tbi_path.exists() and mutrate_tbi_path.stat().st_size > 0
        )
        mutrate_raw_ready = (
            mutrate_raw_path.exists() and mutrate_raw_path.stat().st_size > 0
        )

        # Step 1: bam-readcount
        # Skip if counts already exist OR if any downstream mutrate artifact is
        # already present (counts is deleted after Step 2a, and the raw .gz is
        # deleted after Step 2b, so a missing counts file alongside an existing
        # mutrate .gz/.bgz is the expected steady state).
        if counts_path.exists() and counts_path.stat().st_size > 0:
            print(f"{log_prefix} Step 1 skipped (exists): {counts_path.name}", flush=True)
        elif mutrate_ready:
            print(f"{log_prefix} Step 1 skipped (mutrate exists): {mutrate_path.name}", flush=True)
        elif mutrate_raw_ready:
            print(f"{log_prefix} Step 1 skipped (raw mutrate exists): {mutrate_raw_path.name}", flush=True)
        else:
            # bam-readcount writes per-position counts to stdout, piped through pigz.
            brc_cmd = ["bam-readcount", "-w", "0", "-f", args.ref_fasta]
            if args.regions_bed:
                brc_cmd += ["-l", args.regions_bed]
            brc_cmd += [str(bam_path)]
            pigz_cmd = ["pigz", "-p", str(args.threads), "-c"]
            with open(counts_path, "wb") as out_fh:
                print(f"{log_prefix} Step 1 bam-readcount | pigz -> {counts_path.name}", flush=True)
                brc = subprocess.Popen(brc_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
                pigz = subprocess.Popen(pigz_cmd, stdin=brc.stdout, stdout=out_fh, stderr=subprocess.DEVNULL)
                brc.stdout.close()  # let brc receive SIGPIPE if pigz exits
                pigz_rc = pigz.wait()
                brc_rc = brc.wait()
            if brc_rc != 0:
                raise RuntimeError(f"bam-readcount failed (exit {brc_rc})")
            if pigz_rc != 0:
                raise RuntimeError(f"pigz failed (exit {pigz_rc})")

        # Step 2a: get_mutrate.v4.1.py -> raw .mutrate.txt.gz
        if mutrate_ready:
            print(f"{log_prefix} Step 2 skipped (exists): {mutrate_path.name}", flush=True)
        else:
            if mutrate_raw_ready:
                print(f"{log_prefix} Step 2a skipped (raw exists): {mutrate_raw_path.name}", flush=True)
            else:
                cmd = [
                    args.python, str(get_mutrate),
                    "-i", str(counts_path),
                    "-o", str(mutrate_raw_path),
                    "-w", str(args.get_mutrate_w),
                ]
                if args.gene_readcounts:
                    cmd += ["-e", args.gene_readcounts]
                if args.ref_symbol:
                    cmd += ["--ref_symbol"]
                run_step("get_mutrate.v4.1", cmd, log_prefix)

                # Drop the (often huge) per-position counts file once mutrate is built.
                if counts_path.exists():
                    try:
                        counts_path.unlink()
                        print(f"{log_prefix} removed intermediate {counts_path.name}", flush=True)
                    except OSError as e:
                        print(f"{log_prefix} WARN: could not remove {counts_path.name}: {e}", flush=True)

            # Step 2b: sort + bgzip + tabix index -> .mutrate.txt.bgz(+.tbi)
            sort_bgzip_index(mutrate_raw_path, mutrate_path, log_prefix)

            # Drop the unsorted .gz now that the indexed .bgz is in place.
            if mutrate_raw_path.exists():
                try:
                    mutrate_raw_path.unlink()
                    print(f"{log_prefix} removed intermediate {mutrate_raw_path.name}", flush=True)
                except OSError as e:
                    print(f"{log_prefix} WARN: could not remove {mutrate_raw_path.name}: {e}", flush=True)

        # Step 3: mutrate2window.py
        if window_path.exists() and window_path.stat().st_size > 0:
            print(f"{log_prefix} Step 3 skipped (exists): {window_path.name}", flush=True)
        else:
            cmd = [
                args.python, str(mutrate2win),
                "-i", str(mutrate_path),
                "-o", str(window_path),
                "-w", str(args.win_len),
                "-s", str(args.step),
                "--coverage-threshold", str(args.coverage_threshold),
                "--mutrate-threshold", str(args.mutrate_threshold),
            ]
            run_step("mutrate2window", cmd, log_prefix)

        print(f"{log_prefix} DONE -> {window_path}", flush=True)
        return bam_path, True, None

    except Exception as e:
        print(f"{log_prefix} FAILED: {e}", flush=True)
        return bam_path, False, str(e)


def main():
    args = parse_args()

    # Resolve everything to absolute paths so subprocesses don't depend on cwd
    args.ref_fasta = str(Path(args.ref_fasta).resolve())
    args.output_dir = str(Path(args.output_dir).resolve())
    if args.regions_bed:
        args.regions_bed = str(Path(args.regions_bed).resolve())
    if args.gene_readcounts:
        args.gene_readcounts = str(Path(args.gene_readcounts).resolve())

    if not Path(args.ref_fasta).is_file():
        sys.exit(f"ERROR: --ref-fasta not found: {args.ref_fasta}")
    if args.regions_bed and not Path(args.regions_bed).is_file():
        sys.exit(f"ERROR: --regions-bed not found: {args.regions_bed}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    get_mutrate, mutrate2win = check_dependencies(args.shapemap_util_dir)
    bams = find_bams(args.input, args.bam_glob)

    print("=" * 60)
    print("BAM -> window mutation rate pipeline")
    print(f"  inputs:           {len(bams)} BAM file(s)")
    print(f"  output dir:       {args.output_dir}")
    print(f"  ref fasta:        {args.ref_fasta}")
    print(f"  regions bed:      {args.regions_bed if args.regions_bed else '(none; whole BAM)'}")
    print(f"  win/step:         {args.win_len}/{args.step}")
    print(f"  cov threshold:    {args.coverage_threshold}")
    print(f"  mutrate threshold:{args.mutrate_threshold}")
    print(f"  parallel BAMs:    {args.threads}")
    print("=" * 60)

    successes, failures = [], []

    if args.threads <= 1 or len(bams) == 1:
        for bam in bams:
            _, ok, err = process_bam(bam, args, get_mutrate, mutrate2win)
            (successes if ok else failures).append((bam, err))
    else:
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = {ex.submit(process_bam, bam, args, get_mutrate, mutrate2win): bam
                    for bam in bams}
            for fut in as_completed(futs):
                bam, ok, err = fut.result()
                (successes if ok else failures).append((bam, err))

    print("=" * 60)
    print(f"Done: {len(successes)} succeeded, {len(failures)} failed")
    if failures:
        for bam, err in failures:
            print(f"  FAIL  {bam}: {err}")
        sys.exit(1)


if __name__ == "__main__":
    main()
