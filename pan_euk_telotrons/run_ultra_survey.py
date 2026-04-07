#!/usr/bin/env python3
"""
Pan-eukaryotic telotron survey v2 — ULTRA-based microsatellite detection.

For each NCBI genome:
  1. Download genome + GFF3
  2. Extract introns from GFF3
  3. Run ULTRA (-p 10) on intron sequences
  4. Post-filter for telomeric repeat consensus
  5. Checkpoint and cleanup

All tandem repeats are captured; telomeric filtering is post-hoc.
"""

import json, os, re, sys, shutil, subprocess, tempfile, csv
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASEDIR     = Path(__file__).resolve().parent
DATASETS    = BASEDIR.parent / "pan_euk_telotrons_datasets"
ULTRA       = BASEDIR / "bin" / "ultra"
CHECKPOINT  = BASEDIR / "ultra_checkpoint.json"
WORKERS     = 3       # genomes in parallel (reduce to avoid download contention)
ULTRA_THREADS = 16    # threads per ULTRA instance

DIR_GENOMES   = BASEDIR / "genomes"
DIR_RESULTS   = BASEDIR / "ultra_results"     # per-species telotron TSVs
DIR_ALL_TR    = BASEDIR / "ultra_all_tr"      # optional: all tandem repeats

# ---------------------------------------------------------------------------
# Known telomeric repeats (all rotations + reverse complement)
# ---------------------------------------------------------------------------
def _rotations(base):
    s = set()
    for i in range(len(base)):
        s.add(base[i:] + base[:i])
    rc = base[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    for i in range(len(rc)):
        s.add(rc[i:] + rc[:i])
    return s

def _fwd_rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

def _rev_rotations(base):
    rc = base[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    return {rc[i:] + rc[:i] for i in range(len(rc))}

TELO_REPEATS = {
    'TTAGGG':     {'all': _rotations('TTAGGG'),     'fwd': _fwd_rotations('TTAGGG'),     'rev': _rev_rotations('TTAGGG')},
    'TTTAGGG':    {'all': _rotations('TTTAGGG'),    'fwd': _fwd_rotations('TTTAGGG'),    'rev': _rev_rotations('TTTAGGG')},
    'TTTTAGGG':   {'all': _rotations('TTTTAGGG'),   'fwd': _fwd_rotations('TTTTAGGG'),   'rev': _rev_rotations('TTTTAGGG')},
    'TTTGGG':     {'all': _rotations('TTTGGG'),     'fwd': _fwd_rotations('TTTGGG'),     'rev': _rev_rotations('TTTGGG')},
    'TTGGGG':     {'all': _rotations('TTGGGG'),     'fwd': _fwd_rotations('TTGGGG'),     'rev': _rev_rotations('TTGGGG')},
    'TTTGGGG':    {'all': _rotations('TTTGGGG'),    'fwd': _fwd_rotations('TTTGGGG'),    'rev': _rev_rotations('TTTGGGG')},
    'TTTTGGGG':   {'all': _rotations('TTTTGGGG'),   'fwd': _fwd_rotations('TTTTGGGG'),   'rev': _rev_rotations('TTTTGGGG')},
    'TTAGG':      {'all': _rotations('TTAGG'),      'fwd': _fwd_rotations('TTAGG'),      'rev': _rev_rotations('TTAGG')},
    'TTAGGC':     {'all': _rotations('TTAGGC'),     'fwd': _fwd_rotations('TTAGGC'),     'rev': _rev_rotations('TTAGGC')},
    'TAGGG':      {'all': _rotations('TAGGG'),      'fwd': _fwd_rotations('TAGGG'),      'rev': _rev_rotations('TAGGG')},
    'TTCAGG':     {'all': _rotations('TTCAGG'),     'fwd': _fwd_rotations('TTCAGG'),     'rev': _rev_rotations('TTCAGG')},
    'TTGGG':      {'all': _rotations('TTGGG'),      'fwd': _fwd_rotations('TTGGG'),      'rev': _rev_rotations('TTGGG')},
    'TTTTTTAGGG': {'all': _rotations('TTTTTTAGGG'), 'fwd': _fwd_rotations('TTTTTTAGGG'), 'rev': _rev_rotations('TTTTTTAGGG')},
}

# Flatten all telo kmers for quick membership check
ALL_TELO_KMERS = set()
for rsets in TELO_REPEATS.values():
    ALL_TELO_KMERS |= rsets['all']

# ---------------------------------------------------------------------------
# Checkpoint
# ---------------------------------------------------------------------------
def load_checkpoint():
    if CHECKPOINT.exists():
        return json.loads(CHECKPOINT.read_text())
    return {}

def save_checkpoint(ckpt):
    tmp = CHECKPOINT.with_suffix('.tmp')
    tmp.write_text(json.dumps(ckpt, indent=1))
    tmp.replace(CHECKPOINT)

# ---------------------------------------------------------------------------
# FASTA / GFF parsing
# ---------------------------------------------------------------------------
def load_fasta(path):
    seqs = {}
    curr = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if curr:
                    seqs[curr] = ''.join(parts)
                curr = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip().upper())
    if curr:
        seqs[curr] = ''.join(parts)
    return seqs

def extract_introns(gff_path, fna_path, max_intron=100000):
    """Extract intron sequences from GFF + FASTA. Returns list of dicts."""
    genome = load_fasta(fna_path)

    mrna_features = defaultdict(list)
    mrna_strand = {}
    mrna_contig = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] not in ('exon', 'CDS'):
                continue
            attrs = dict(re.findall(r'(\w+)=([^;]+)', parts[8]))
            parent = attrs.get('Parent', '')
            contig = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            mrna_features[parent].append((start, end))
            mrna_strand[parent] = strand
            mrna_contig[parent] = contig

    introns = []
    seen = set()
    for mrna, exons in mrna_features.items():
        exons.sort()
        contig = mrna_contig[mrna]
        strand = mrna_strand[mrna]
        if contig not in genome:
            continue
        seq = genome[contig]
        for i in range(len(exons) - 1):
            istart = exons[i][1]
            iend = exons[i + 1][0]
            if iend <= istart:
                continue
            length = iend - istart
            if length < 20 or length > max_intron:
                continue
            key = (contig, istart, iend)
            if key in seen:
                continue
            seen.add(key)
            intron_seq = seq[istart:iend]
            if not intron_seq:
                continue
            introns.append({
                'contig': contig,
                'start': istart,
                'end': iend,
                'length': length,
                'strand': strand,
                'mrna': mrna,
                'donor': intron_seq[:2],
                'acceptor': intron_seq[-2:],
                'seq': intron_seq,
            })
    return introns

# ---------------------------------------------------------------------------
# Download genome
# ---------------------------------------------------------------------------
def download_genome(acc):
    out_dir = DIR_GENOMES / acc
    zip_path = DIR_GENOMES / f"{acc}.zip"

    if out_dir.exists():
        gff = _find(out_dir, "*.gff")
        fna = _find(out_dir, "*.fna")
        if gff and fna:
            return (gff, fna)

    try:
        subprocess.run(
            [str(DATASETS), "download", "genome", "accession", acc,
             "--include", "genome,gff3", "--filename", str(zip_path)],
            check=True, capture_output=True, timeout=3600,
        )
    except Exception:
        return None

    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(["unzip", "-qo", str(zip_path), "-d", str(out_dir)],
                        check=True, capture_output=True, timeout=600)
    except Exception:
        return None
    finally:
        zip_path.unlink(missing_ok=True)

    gff = _find(out_dir, "*.gff")
    fna = _find(out_dir, "*.fna")
    return (gff, fna) if gff and fna else None

def _find(root, pattern):
    matches = list(root.rglob(pattern))
    return str(matches[0]) if matches else None

# ---------------------------------------------------------------------------
# Run ULTRA on intron sequences
# ---------------------------------------------------------------------------
def run_ultra_on_introns(introns, tmp_dir):
    """Write introns to FASTA, run ULTRA, parse output.
    Returns list of (intron_idx, ultra_hit_dict) tuples."""

    # Write FASTA
    fa_path = os.path.join(tmp_dir, 'introns.fa')
    idx_map = {}  # fasta_id -> intron_index
    with open(fa_path, 'w') as f:
        for i, intr in enumerate(introns):
            fid = f"intron_{i}__{intr['contig']}:{intr['start']}-{intr['end']}"
            idx_map[fid] = i
            f.write(f">{fid}\n{intr['seq']}\n")

    # Run ULTRA
    try:
        result = subprocess.run(
            [str(ULTRA), '-p', '10', '-t', str(ULTRA_THREADS),
             '--tsv', '-c', '--show_seq', fa_path],
            capture_output=True, text=True, timeout=1800,
        )
    except subprocess.TimeoutExpired:
        return []

    if result.returncode != 0:
        return []

    # Parse TSV output
    hits = []
    header = None
    for line in result.stdout.split('\n'):
        if line.startswith('{') or line.startswith('---') or line.startswith('*') or not line.strip():
            continue
        if line.startswith('SeqID'):
            header = line.strip().split('\t')
            continue
        if header is None:
            continue

        fields = line.strip().split('\t')
        if len(fields) < len(header):
            continue

        row = dict(zip(header, fields))
        seq_id = row.get('SeqID', '')

        # Match back to intron
        if seq_id in idx_map:
            intron_idx = idx_map[seq_id]
        else:
            continue

        try:
            period = int(row.get('Period', 0))
            score = float(row.get('Score', 0))
            start = int(row.get('Start', 0))
            end = int(row.get('End', 0))
            consensus = row.get('Consensus', '')
            copies = float(row.get('#copies', 0))
            subs = int(row.get('#substitutions', 0))
            ins = int(row.get('#insertions', 0))
            dels = int(row.get('#deletions', 0))
        except (ValueError, TypeError):
            continue

        tr_len = end - start
        if tr_len < 10:
            continue

        hits.append({
            'intron_idx': intron_idx,
            'tr_start': start,
            'tr_end': end,
            'tr_len': tr_len,
            'period': period,
            'consensus': consensus.upper(),
            'copies': copies,
            'score': score,
            'substitutions': subs,
            'insertions': ins,
            'deletions': dels,
        })

    return hits

# ---------------------------------------------------------------------------
# Classify ULTRA hits as telomeric or not
# ---------------------------------------------------------------------------
def classify_telomeric(consensus, period):
    """Check if ULTRA consensus matches any known telomeric repeat."""
    cons = consensus.upper()

    # Check all rotations of consensus against known telo kmers
    for rtype, rsets in TELO_REPEATS.items():
        if period != len(rtype):
            continue
        # Check if consensus (or any rotation) is in the set
        for i in range(len(cons)):
            rotated = cons[i:] + cons[:i]
            if rotated in rsets['all']:
                # Determine orientation
                if rotated in rsets['fwd']:
                    return rtype, 'coding-strand'
                elif rotated in rsets['rev']:
                    return rtype, 'template-strand'
                else:
                    return rtype, 'unknown'

    # Fuzzy match: allow 1 mismatch in consensus vs any telo repeat
    for rtype, rsets in TELO_REPEATS.items():
        if period != len(rtype):
            continue
        for telo_kmer in rsets['all']:
            if len(cons) != len(telo_kmer):
                continue
            mismatches = sum(a != b for a, b in zip(cons, telo_kmer))
            if mismatches <= 1:
                orient = 'coding-strand' if telo_kmer in rsets['fwd'] else 'template-strand'
                return rtype + '~1mm', orient

    return None, None

# ---------------------------------------------------------------------------
# Process one genome
# ---------------------------------------------------------------------------
def process_one_genome(acc, name):
    result = {
        'acc': acc, 'name': name,
        'status': 'failed', 'error': None,
        'timestamp': datetime.now().isoformat(),
    }

    tmp_dir = tempfile.mkdtemp(prefix=f'ultra_{acc}_')

    try:
        # Download
        paths = download_genome(acc)
        if paths is None:
            result['error'] = 'download_failed'
            return result
        gff_path, fna_path = paths

        # Extract introns
        introns = extract_introns(gff_path, fna_path)
        result['n_introns'] = len(introns)

        if len(introns) == 0:
            result['status'] = 'ok'
            result['n_telotrons'] = 0
            result['n_tr_introns'] = 0
            return result

        # Run ULTRA
        ultra_hits = run_ultra_on_introns(introns, tmp_dir)
        result['n_ultra_hits'] = len(ultra_hits)

        # For each intron, find the dominant TR and check if telomeric
        intron_tr = defaultdict(list)
        for hit in ultra_hits:
            intron_tr[hit['intron_idx']].append(hit)

        result['n_tr_introns'] = len(intron_tr)

        # Classify telomeric hits
        telotrons = []
        for intron_idx, tr_hits in intron_tr.items():
            intr = introns[intron_idx]
            intron_len = intr['length']

            # Find the best (longest) TR hit for this intron
            best = max(tr_hits, key=lambda h: h['tr_len'])

            # Check if telomeric
            telo_type, orient = classify_telomeric(best['consensus'], best['period'])
            if telo_type is None:
                continue

            # Calculate coverage: total TR bases / intron length
            tr_coverage = best['tr_len'] / intron_len

            # Also check for converging: multiple TR hits with different orientations
            orientations_seen = set()
            for h in tr_hits:
                tt, ot = classify_telomeric(h['consensus'], h['period'])
                if tt and ot:
                    orientations_seen.add(ot)

            if len(orientations_seen) >= 2:
                final_orient = 'converging'
            else:
                final_orient = orient

            telotrons.append({
                'contig': intr['contig'],
                'start': intr['start'],
                'end': intr['end'],
                'length': intron_len,
                'strand': intr['strand'],
                'mrna': intr['mrna'],
                'donor': intr['donor'],
                'acceptor': intr['acceptor'],
                'repeat_type': telo_type,
                'period': best['period'],
                'consensus': best['consensus'],
                'copies': best['copies'],
                'tr_coverage': round(tr_coverage, 4),
                'orientation': final_orient,
                'score': best['score'],
                'substitutions': best['substitutions'],
                'insertions': best['insertions'],
                'deletions': best['deletions'],
                'intron_seq': intr['seq'],
            })

        result['n_telotrons'] = len(telotrons)

        # Write telotron TSV
        if telotrons:
            out_path = DIR_RESULTS / f"{acc}_{name}_telotrons.tsv"
            cols = ['contig', 'start', 'end', 'length', 'strand', 'mrna',
                    'donor', 'acceptor', 'repeat_type', 'period', 'consensus',
                    'copies', 'tr_coverage', 'orientation', 'score',
                    'substitutions', 'insertions', 'deletions', 'intron_seq']
            with open(out_path, 'w') as f:
                f.write('genome_acc\tgenome_name\t' + '\t'.join(cols) + '\n')
                for t in telotrons:
                    vals = [str(t.get(c, '')) for c in cols]
                    f.write(f'{acc}\t{name}\t' + '\t'.join(vals) + '\n')

        # Summary stats
        if telotrons:
            repeat_counts = Counter(t['repeat_type'] for t in telotrons)
            orient_counts = Counter(t['orientation'] for t in telotrons)
            result['repeat_types'] = dict(repeat_counts)
            result['orientations'] = dict(orient_counts)
            result['mean_coverage'] = round(
                sum(t['tr_coverage'] for t in telotrons) / len(telotrons), 3)

        result['status'] = 'ok'

    except Exception as e:
        import traceback
        result['error'] = str(e)
        result['traceback'] = traceback.format_exc()

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        # Cleanup genome
        if '--keep-genomes' not in sys.argv:
            genome_dir = DIR_GENOMES / acc
            if genome_dir.exists():
                shutil.rmtree(genome_dir, ignore_errors=True)
            (DIR_GENOMES / f"{acc}.zip").unlink(missing_ok=True)

    return result

# ---------------------------------------------------------------------------
# Build genome list from previous checkpoint
# ---------------------------------------------------------------------------
def get_genome_list():
    """Get list of (acc, name) from the previous survey checkpoint."""
    old_ckpt = BASEDIR / "checkpoint.json"
    if old_ckpt.exists():
        ckpt = json.loads(old_ckpt.read_text())
        genomes = []
        for acc, r in ckpt.items():
            if r.get('status') == 'ok':
                genomes.append((acc, r.get('name', '?')))
        return genomes

    # Fallback: scan NCBI for annotated eukaryotic genomes
    print("No previous checkpoint found. Use checkpoint.json from v1 survey.")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    for d in [DIR_GENOMES, DIR_RESULTS]:
        d.mkdir(parents=True, exist_ok=True)

    ckpt = load_checkpoint()
    all_genomes = get_genome_list()

    # Filter to unprocessed
    if '--force' in sys.argv:
        todo = all_genomes
    elif '--only' in sys.argv:
        idx = sys.argv.index('--only') + 1
        target = sys.argv[idx] if idx < len(sys.argv) else ''
        todo = [(a, n) for a, n in all_genomes
                if target.lower() in a.lower() or target.lower() in n.lower()]
    else:
        todo = [(a, n) for a, n in all_genomes if a not in ckpt]

    if not todo:
        print(f"All {len(all_genomes)} genomes processed. Use --force to rerun.")
        return

    print(f"{'='*70}")
    print(f" ULTRA-based Pan-Eukaryotic Telotron Survey v2")
    print(f" Genomes: {len(all_genomes)} total, {len(ckpt)} done, {len(todo)} to process")
    print(f" ULTRA: period ≤10, {ULTRA_THREADS} threads per genome")
    print(f" Telomeric repeats: {len(TELO_REPEATS)} types searched")
    print(f" Workers: {WORKERS}")
    print(f"{'='*70}")
    print()

    t0 = datetime.now()
    n_done = len(ckpt)
    n_submitted = 0

    import threading
    ckpt_lock = threading.Lock()

    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = {}
        for acc, name in todo:
            fut = pool.submit(process_one_genome, acc, name)
            futures[fut] = (acc, name)
            n_submitted += 1

        for fut in as_completed(futures):
            acc, name = futures[fut]
            try:
                result = fut.result()
            except Exception as e:
                result = {'acc': acc, 'name': name, 'status': 'crashed',
                          'error': str(e), 'n_introns': 0, 'n_telotrons': 0,
                          'n_tr_introns': 0}

            with ckpt_lock:
                ckpt[acc] = result
                save_checkpoint(ckpt)
                n_done += 1

            elapsed = (datetime.now() - t0).total_seconds() / 3600
            rate = (n_done - len(load_checkpoint()) + n_submitted) / elapsed if elapsed > 0.001 else 0

            status = result.get('status', '?')
            n_int = result.get('n_introns', 0)
            n_tel = result.get('n_telotrons', 0)
            n_tr = result.get('n_tr_introns', 0)

            if n_tel > 0:
                tag = f"HIT({n_tel})"
                reps = result.get('repeat_types', {})
                rep_str = ','.join(f"{k}:{v}" for k, v in sorted(reps.items(), key=lambda x: -x[1])[:5])
            else:
                tag = 'OK' if status == 'ok' else 'FAIL'
                rep_str = ''

            print(f"  [{n_done:>4d}/{len(all_genomes)}] "
                  f"[{tag:<20s}] "
                  f"{name:.<45s} "
                  f"introns={n_int:>9,}  tr={n_tr:>6,}  telo={n_tel:>4,}  "
                  f"{rep_str}",
                  flush=True)

    # Final summary
    ckpt = load_checkpoint()
    n_ok = sum(1 for v in ckpt.values() if v.get('status') == 'ok')
    n_telo = sum(1 for v in ckpt.values() if v.get('n_telotrons', 0) > 0)
    total_introns = sum(v.get('n_introns', 0) for v in ckpt.values())
    total_telo = sum(v.get('n_telotrons', 0) for v in ckpt.values())

    print(f"\n{'='*70}")
    print(f" COMPLETE")
    print(f" Genomes processed: {n_ok}")
    print(f" Total introns: {total_introns:,}")
    print(f" Species with telotrons: {n_telo}")
    print(f" Total telotron candidates: {total_telo}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
