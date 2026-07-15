#!/usr/bin/env python
"""
Per-read mutation vs. RNA-structure classifier (sparse output).

For every read in a BAM, find its mutations (mismatch + short indels, same
definition as get_mutrate.v4.1: indel runs longer than --max-indel-len are
ignored) and record, per reference position, where each read is mutated.
Output is a SPARSE reads x positions matrix per contig, so it scales to very
large read counts.

A reference position is a POSITIVE SHAPE site if it is:
    - unpaired ('.'), OR
    - a junction: a paired base within --junction-window nt of a '.'
Positivity is a property of the position (column), stored in the per-position
annotation file, so it is not baked into the mutation matrix.

Mutations are read from the MD tag (mismatch/deletion) and the CIGAR
(insertion). No reference FASTA is needed as long as reads carry MD tags
(samtools calmd / most aligners add these). Soft-clips (S) from local
alignment are correctly NOT counted as insertions.

Structure coordinates must match the BAM reference. For mouse 18S the
insmatchref CSV (1870 nt) matches the reference; the raw dot-bracket (1695 nt)
does NOT -- use type 'csv' for 18S and 'dotbracket' for 12S/16S.

Contig names on the left of --structure must match your BAM header exactly
(check with: samtools idxstats your.bam). In the 7th_Round spatial data these
are gene symbols: mt-Rnr1 (12S), mt-Rnr2 (16S), 18S_mouse_7cpu (18S),
7SK_ENSMUSG00002076161.1 (7SK).

Outputs (under --output-dir, prefixed with --output-prefix):
    {prefix}.{contig}.mut.npz        scipy CSR int8, rows=reads, cols=positions
                                     values: 1=mismatch, 2=insertion, 3=deletion
    {prefix}.{contig}.reads.txt.gz   read IDs, one per line, in matrix row order
    {prefix}.{contig}.positions.tsv  pos(1-based), ref_base, symbol,
                                     is_unpaired, is_junction, is_positive
    {prefix}.read_summary.tsv.gz     one row per read: mutation & positive counts

Usage:
    python bam_mutation_structure.py CID_found_R1_filter40_dedup.bam \\
        --dotbracket-file "/home/users/astar/gis/linxy/scratch/data/references/Dot bracket for ribo for mouse.txt" \\
        --structure "18S_mouse_7cpu=csv:/home/users/astar/gis/linxy/scratch/data/references/cliffzhang_mouse18s_insmatchref.csv:dbrt" \\
        --structure "mt-Rnr1=dotbracket:12S_mt" \\
        --structure "mt-Rnr2=dotbracket:16S_mt" \\
        --structure "7SK_ENSMUSG00002076161.1=dotbracket:7SK" \\
        -o Chow_biorep1_2A3 --output-dir ./mut_structure

Load the result:
    import scipy.sparse as sp, numpy as np, pandas as pd, gzip
    M   = sp.load_npz('Chow_biorep1_2A3.mt-Rnr1.mut.npz')      # reads x positions
    pos = pd.read_csv('Chow_biorep1_2A3.mt-Rnr1.positions.tsv', sep='\\t')
    ids = [l.strip() for l in gzip.open('Chow_biorep1_2A3.mt-Rnr1.reads.txt.gz','rt')]
    # positive mutations per read = (M != 0) restricted to positive columns
    # (is_positive is stored as 0/1, so cast to bool before masking):
    pos_cols = np.flatnonzero(pos['is_positive'].to_numpy().astype(bool))
    n_pos_per_read = (M[:, pos_cols] != 0).sum(axis=1).A1
"""

import argparse
import csv
import gzip
import os
import sys

import numpy as np
import pysam
import scipy.sparse as sp

# mutation type codes stored in the sparse matrix
MISMATCH, INSERTION, DELETION = 1, 2, 3
# pysam CIGAR op codes
CIG_MATCH, CIG_INS, CIG_DEL, CIG_REFSKIP = 0, 1, 2, 3
CIG_SOFT, CIG_HARD, CIG_PAD, CIG_EQUAL, CIG_DIFF = 4, 5, 6, 7, 8
CIG_CONSUMES_REF = {CIG_MATCH, CIG_DEL, CIG_REFSKIP, CIG_EQUAL, CIG_DIFF}
PROGRESS_INTERVAL = 500000


def read_dot_bracket(file_path, gene_name):
    """Extract the dot-bracket string for `gene_name` from a >header file."""
    with open(file_path) as f:
        lines = f.readlines()
    target = f'>{gene_name}'
    for i, line in enumerate(lines):
        if line.strip().startswith(target):
            db = ''
            for j in range(i + 1, len(lines)):
                nxt = lines[j].strip()
                if nxt.startswith('>'):
                    break
                db += nxt
            return db or None
    return None


def load_structures(structure_specs, dotbracket_file):
    """Parse --structure specs into {contig: dot_bracket_string}."""
    import pandas as pd

    structures = {}
    for spec in structure_specs:
        if '=' not in spec:
            raise ValueError(f"Bad --structure (missing '='): {spec}")
        contig, rhs = spec.split('=', 1)
        kind = rhs.split(':', 1)[0]
        if kind == 'dotbracket':
            gene = rhs.split(':', 1)[1]
            if dotbracket_file is None:
                raise ValueError("--dotbracket-file required for a 'dotbracket' structure")
            struct = read_dot_bracket(dotbracket_file, gene)
            if struct is None:
                raise ValueError(f"Gene '{gene}' not found in {dotbracket_file}")
        elif kind == 'csv':
            _, path, column = rhs.split(':', 2)
            struct = ''.join(pd.read_csv(path)[column].astype(str).tolist())
        else:
            raise ValueError(f"Unknown structure type '{kind}' in: {spec}")
        structures[contig] = struct
    return structures


def build_positive_site_mask(structure, junction_window=1,
                             unpaired_symbol='.', paired_symbols="()<>[]{}"):
    """Classify each position as unpaired / junction / positive (0-based arrays)."""
    n = len(structure)
    symbols = list(structure)
    is_unpaired = np.array([c == unpaired_symbol for c in symbols])
    is_paired = np.array([c in paired_symbols for c in symbols])

    is_junction = np.zeros(n, dtype=bool)
    unpaired_idx = np.flatnonzero(is_unpaired)
    for i in np.flatnonzero(is_paired):
        lo, hi = i - junction_window, i + junction_window
        if np.any((unpaired_idx >= lo) & (unpaired_idx <= hi)):
            is_junction[i] = True

    is_positive = is_unpaired | is_junction
    return symbols, is_unpaired, is_junction, is_positive


def extract_read_mutations(read, n_struct, ref_seq, max_indel_len, min_base_qual,
                           count_mismatch, count_insertion, count_deletion):
    """
    Return {ref_pos_0based: type_code} for one read (one code/position, taking
    the larger code on collision: deletion > insertion > mismatch).

    Mismatches/deletions come from the MD tag via get_aligned_pairs(with_seq)
    when present (falls back to `ref_seq` comparison otherwise). Insertions come
    strictly from the CIGAR, so soft-clips (S) are never miscounted as inserts.
    """
    muts = {}
    qqual = read.query_qualities
    qseq = read.query_sequence

    ## --- mismatches + deletions ---
    if count_mismatch or count_deletion:
        if read.has_tag("MD"):
            pairs = read.get_aligned_pairs(with_seq=True)
            i, npairs = 0, len(pairs)
            while i < npairs:
                qpos, rpos, rb = pairs[i]
                if qpos is not None and rpos is not None:
                    # with_seq lowercases the ref base at a mismatch
                    if count_mismatch and rb is not None and rb.islower() \
                            and rpos < n_struct:
                        if qqual is None or qqual[qpos] >= min_base_qual:
                            muts[rpos] = max(muts.get(rpos, 0), MISMATCH)
                    i += 1
                elif qpos is None and rpos is not None:
                    run = []
                    while i < npairs and pairs[i][0] is None and pairs[i][1] is not None:
                        run.append(pairs[i][1])
                        i += 1
                    if count_deletion and len(run) <= max_indel_len:
                        for r in run:
                            if r < n_struct:
                                muts[r] = max(muts.get(r, 0), DELETION)
                else:
                    i += 1  # insertion/soft-clip in with_seq; handled below/skipped
        elif ref_seq is not None:
            qp, rp = 0, read.reference_start
            for op, length in read.cigartuples or []:
                if op in (CIG_MATCH, CIG_EQUAL, CIG_DIFF):
                    if count_mismatch:
                        for k in range(length):
                            rpos = rp + k
                            if rpos < n_struct and rpos < len(ref_seq):
                                rbb = ref_seq[rpos]
                                qbb = qseq[qp + k].upper()
                                if qbb != rbb and rbb != 'N' and qbb != 'N':
                                    if qqual is None or qqual[qp + k] >= min_base_qual:
                                        muts[rpos] = max(muts.get(rpos, 0), MISMATCH)
                    qp += length
                    rp += length
                elif op == CIG_DEL:
                    if count_deletion and length <= max_indel_len:
                        for k in range(length):
                            if rp + k < n_struct:
                                muts[rp + k] = max(muts.get(rp + k, 0), DELETION)
                    rp += length
                elif op == CIG_REFSKIP:
                    rp += length
                elif op in (CIG_INS, CIG_SOFT):
                    qp += length
        else:
            raise ValueError(
                f"read {read.query_name} lacks an MD tag and no --fasta given; "
                "run `samtools calmd` or pass --fasta")

    ## --- insertions strictly from the CIGAR (never a soft-clip) ---
    if count_insertion and read.cigartuples:
        rp = read.reference_start
        for op, length in read.cigartuples:
            if op in CIG_CONSUMES_REF:
                rp += length
            elif op == CIG_INS:
                if length <= max_indel_len:
                    anchor = rp - 1 if rp - 1 >= 0 else rp  # ref base before insert
                    if 0 <= anchor < n_struct:
                        muts[anchor] = max(muts.get(anchor, 0), INSERTION)
            # soft/hard/pad consume no reference -> nothing to advance
    return muts


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam', help='input BAM (any sort order)')
    parser.add_argument('-f', '--fasta', default=None,
                        help='OPTIONAL reference FASTA (.fai); only used for the '
                             'ref_base annotation column and as a fallback when a '
                             'read has no MD tag')
    parser.add_argument('--structure', action='append', required=True,
                        metavar='CONTIG=SPEC',
                        help="repeatable; 'CONTIG=dotbracket:GENE' or 'CONTIG=csv:PATH:COL'")
    parser.add_argument('--dotbracket-file', default=None,
                        help="source file for 'dotbracket' structure specs")
    parser.add_argument('-o', '--output-prefix', default='mut_structure')
    parser.add_argument('--output-dir', default='.')
    parser.add_argument('--junction-window', type=int, default=1,
                        help='paired base within this many nt of a "." is a junction (default 1)')
    parser.add_argument('--paired-symbols', default="()<>[]{}",
                        help='symbols counted as paired (default "()<>[]{}")')
    parser.add_argument('--max-indel-len', type=int, default=3,
                        help='indel runs longer than this are ignored (default 3)')
    parser.add_argument('--min-base-qual', type=int, default=0,
                        help='minimum base quality for a mismatch (default 0)')
    parser.add_argument('--no-mismatch', action='store_true')
    parser.add_argument('--no-insertion', action='store_true')
    parser.add_argument('--no-deletion', action='store_true')
    parser.add_argument('--keep-secondary', action='store_true',
                        help='also process secondary/supplementary alignments')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    prefix = os.path.join(args.output_dir, args.output_prefix)

    structures = load_structures(args.structure, args.dotbracket_file)

    fasta = pysam.FastaFile(args.fasta) if args.fasta else None
    bam = pysam.AlignmentFile(args.bam, 'rb')
    bam_lengths = dict(zip(bam.references, bam.lengths))

    # precompute positive-site masks; keep only contigs present in the BAM
    masks, ref_seq = {}, {}
    for contig, struct in structures.items():
        if contig not in bam_lengths:
            print(f"Warning: contig '{contig}' not in BAM header, skipping",
                  file=sys.stderr)
            continue
        if len(struct) != bam_lengths[contig]:
            print(f"Warning: {contig} structure length {len(struct)} != BAM "
                  f"reference length {bam_lengths[contig]}; coordinates may be "
                  "misaligned", file=sys.stderr)
        masks[contig] = build_positive_site_mask(
            struct, junction_window=args.junction_window,
            paired_symbols=args.paired_symbols)
        if fasta is not None and contig in fasta.references:
            ref_seq[contig] = fasta.fetch(contig).upper()
        else:
            ref_seq[contig] = None
    if not masks:
        sys.exit("No requested contigs found in the BAM header.")

    # per-contig sparse triplets and row-ordered read IDs
    triplets = {c: ([], [], []) for c in masks}
    read_ids = {c: [] for c in masks}
    row_counter = {c: 0 for c in masks}

    summary_fh = gzip.open(f"{prefix}.read_summary.tsv.gz", 'wt', newline='')
    summary_writer = csv.writer(summary_fh, delimiter='\t')
    summary_writer.writerow(['read_id', 'contig', 'ref_start', 'ref_end',
                             'n_mutation', 'n_positive', 'n_unpaired',
                             'n_junction', 'n_negative', 'has_positive_mut'])

    count_mismatch = not args.no_mismatch
    count_insertion = not args.no_insertion
    count_deletion = not args.no_deletion

    processed = 0
    for read in bam:
        if read.is_unmapped:
            continue
        if not args.keep_secondary and (read.is_secondary or read.is_supplementary):
            continue
        contig = read.reference_name
        if contig not in masks:
            continue
        if read.query_sequence is None:
            continue

        symbols, is_unpaired, is_junction, is_positive = masks[contig]
        n_struct = len(symbols)
        muts = extract_read_mutations(
            read, n_struct, ref_seq[contig], args.max_indel_len,
            args.min_base_qual, count_mismatch, count_insertion, count_deletion)

        row = row_counter[contig]
        rows, cols, vals = triplets[contig]
        n_pos = n_unp = n_jun = 0
        for rpos, code in muts.items():
            rows.append(row)
            cols.append(rpos)
            vals.append(code)
            if is_positive[rpos]:
                n_pos += 1
            if is_unpaired[rpos]:
                n_unp += 1
            if is_junction[rpos]:
                n_jun += 1

        read_ids[contig].append(read.query_name)
        row_counter[contig] += 1

        n_mut = len(muts)
        summary_writer.writerow([
            read.query_name, contig, read.reference_start + 1, read.reference_end,
            n_mut, n_pos, n_unp, n_jun, n_mut - n_pos, int(n_pos > 0)])

        processed += 1
        if processed % PROGRESS_INTERVAL == 0:
            print(f"processed {processed:,} reads...", file=sys.stderr)

    summary_fh.close()
    bam.close()

    # write per-contig sparse matrix, read IDs, and position annotation
    for contig in masks:
        symbols, is_unpaired, is_junction, is_positive = masks[contig]
        n_pos_cols = len(symbols)
        n_rows = row_counter[contig]
        rows, cols, vals = triplets[contig]
        mat = sp.coo_matrix(
            (np.array(vals, dtype=np.int8),
             (np.array(rows, dtype=np.int64), np.array(cols, dtype=np.int64))),
            shape=(n_rows, n_pos_cols)).tocsr()
        sp.save_npz(f"{prefix}.{contig}.mut.npz", mat)

        with gzip.open(f"{prefix}.{contig}.reads.txt.gz", 'wt') as fh:
            fh.write('\n'.join(read_ids[contig]))
            if read_ids[contig]:
                fh.write('\n')

        rseq = ref_seq[contig]
        with open(f"{prefix}.{contig}.positions.tsv", 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['pos', 'ref_base', 'symbol',
                        'is_unpaired', 'is_junction', 'is_positive'])
            for i in range(n_pos_cols):
                base = rseq[i] if (rseq is not None and i < len(rseq)) else 'N'
                w.writerow([i + 1, base, symbols[i], int(is_unpaired[i]),
                            int(is_junction[i]), int(is_positive[i])])

        print(f"{contig}: {n_rows:,} reads, {mat.nnz:,} mutations, "
              f"{int(is_positive.sum()):,}/{n_pos_cols:,} positive positions")

    if fasta is not None:
        fasta.close()
    print(f"Done. {processed:,} reads processed. Outputs under: {prefix}.*")


if __name__ == "__main__":
    main()
