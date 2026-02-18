import pandas as pd
import urllib.request
import urllib.parse
import json
import time
import os

API_TIMEOUT = 15  # seconds — prevents indefinite hangs on stalled Ensembl connections

# ── Internal HTTP helpers ──────────────────────────────────────────────────────
# We use urllib (stdlib) instead of requests to avoid a Python 3.14 + urllib3
# SSL incompatibility that causes requests to hang on macOS/certain platforms.

def _get(url, timeout=API_TIMEOUT):
    """GET url, return parsed JSON or None on failure."""
    req = urllib.request.Request(
        url,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        # Re-raise with status so callers can inspect it
        e.status_code = e.code
        raise
    except Exception:
        raise


def _post(url, payload, timeout=API_TIMEOUT):
    """POST JSON payload, return parsed JSON or raise."""
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        url,
        data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        e.status_code = e.code
        raise
    except Exception:
        raise


# ── GC content ────────────────────────────────────────────────────────────────

def calculate_gc(seq):
    if not seq: return 0.0
    g = seq.count('G') + seq.count('g')
    c = seq.count('C') + seq.count('c')
    return (g + c) / len(seq) * 100

def fetch_transcript_gc_batch(transcript_ids, build='hg38'):
    """
    Fetches transcript sequences in batches and calculates GC content.
    Returns a dict: {transcript_id: gc_percent}
    """
    server = "https://rest.ensembl.org" if build == 'hg38' else "https://grch37.rest.ensembl.org"
    ext = "/sequence/id"

    results = {}

    # Process in chunks of 50 (Ensembl POST limit)
    chunk_size = 50
    for i in range(0, len(transcript_ids), chunk_size):
        chunk = transcript_ids[i:i+chunk_size]

        # Clean IDs: remove versions for lookup
        chunk_clean = [tid.split('.')[0] for tid in chunk]

        payload = {
            "ids": chunk_clean,
            "type": "cdna",
        }

        try:
            data = _post(server + ext, payload)

            # Guard: if the response is not a list (e.g. an error dict), skip gracefully.
            if not isinstance(data, list):
                print(f"Unexpected GC response format: {str(data)[:200]}")
                for tid in chunk:
                    results[tid] = None
                continue

            # Build lookup: unversioned_id -> sequence
            seq_map = {item['id']: item.get('seq') for item in data if 'id' in item}

            for original_id, clean_id in zip(chunk, chunk_clean):
                seq = seq_map.get(clean_id)
                results[original_id] = calculate_gc(seq) if seq else None

        except urllib.error.HTTPError as e:
            print(f"Batch GC fetch failed (HTTP {e.code}) for chunk starting at index {i}")
            for tid in chunk:
                results[tid] = None
            continue
        except Exception as e:
            print(f"Error in batch GC fetch: {e}")
            for tid in chunk:
                results[tid] = None

        # Be nice to API
        time.sleep(0.1)

    return results


# ── MANE data loading ──────────────────────────────────────────────────────────

MANE_SUMMARY_FILE = "MANE.summary.txt"

def load_mane_data():
    """Bypassed if file doesn't exist, returns None."""
    if not os.path.exists(MANE_SUMMARY_FILE):
        return None
    df = pd.read_csv(MANE_SUMMARY_FILE, sep="\t")
    return df

def get_mane_transcripts(gene_list, mane_df):
    """
    Filter MANE data for the given genes with prioritization:
    1. MANE Select
    2. MANE Plus Clinical (only if Select is missing)

    Returns a DataFrame with relevant columns.
    """
    valid_genes = set(g.upper().strip() for g in gene_list)

    mask = mane_df['symbol'].str.upper().isin(valid_genes)
    subset = mane_df[mask].copy()

    subset['MANE_status'] = pd.Categorical(subset['MANE_status'], categories=['MANE Select', 'MANE Plus Clinical'], ordered=True)
    subset = subset.sort_values('MANE_status')

    result = subset.drop_duplicates(subset=['symbol'], keep='first')
    result = result[['symbol', 'Ensembl_nuc', 'RefSeq_nuc', 'MANE_status', 'Ensembl_Gene', 'GRCh38_chr']]
    return result

def get_ensembl_canonical(gene_list, build='hg38'):
    """
    Fetch canonical transcripts for genes not in MANE via Ensembl API.
    Returns DataFrame matching MANE structure.
    """
    if not gene_list:
        return pd.DataFrame()

    server = "https://rest.ensembl.org" if build == 'hg38' else "https://grch37.rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens"

    # Special manual mappings for loci or ambiguous symbols
    MANUAL_ALIASES = {
        'IGH': 'IGHG1',
        'IGK': 'IGKC',
        'IGL': 'IGLC1',
    }

    records = []

    def resolve_alias(symbol, server):
        try:
            xref_url = f"{server}/xrefs/symbol/homo_sapiens/{symbol}?object_type=gene"
            r = _get(xref_url)
            if not r:
                return None
            gene_id = r[0]['id']
            lookup_url = f"{server}/lookup/id/{gene_id}?expand=1"
            return _get(lookup_url)
        except Exception as e:
            print(f"Error resolving alias for {symbol}: {e}")
        return None

    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    for batch in chunks(list(gene_list), 200):
        try:
            data = _post(server + ext, {"symbols": batch, "expand": 1})

            for gene in batch:
                info = data.get(gene)

                if not info:
                    if gene in MANUAL_ALIASES:
                        target_alias = MANUAL_ALIASES[gene]
                        try:
                            info = _get(f"{server}/lookup/symbol/homo_sapiens/{target_alias}?expand=1")
                        except Exception:
                            info = None

                    if not info:
                        info = resolve_alias(gene, server)
                        if not info:
                            continue

                canonical_transcript = None
                transcripts = info.get('Transcript', [])

                for t in transcripts:
                    if t.get('is_canonical') == 1:
                        canonical_transcript = t
                        break

                if not canonical_transcript and transcripts:
                    canonical_transcript = transcripts[0]

                if canonical_transcript:
                    records.append({
                        'original_input': gene,
                        'symbol': info.get('display_name', gene),
                        'Ensembl_nuc': canonical_transcript['id'],
                        'RefSeq_nuc': '-',
                        'MANE_status': 'Ensembl Canonical',
                        'Ensembl_Gene': info['id'],
                        'GRCh38_chr': info['seq_region_name']
                    })

        except Exception as e:
            print(f"Error in bulk lookup: {e}")

    if records:
        return pd.DataFrame(records)
    else:
        return pd.DataFrame(columns=['symbol', 'Ensembl_nuc', 'RefSeq_nuc', 'MANE_status', 'Ensembl_Gene', 'GRCh38_chr'])


# ── hg19 gene symbol aliases ───────────────────────────────────────────────────

HG19_ALIASES = {
    'H2BC6': 'HIST1H2BC',
    'H3C2': 'HIST1H3B',
    'H2AC11': 'HIST1H2AC',
    'H2AC6': 'HIST1H2AB',
}

def find_canonical_on_build(gene_symbol, build='hg38'):
    """
    Find the canonical transcript ID for a given gene on a specific build.
    """
    if build == 'hg38':
        server = "https://rest.ensembl.org"
    else:
        server = "https://grch37.rest.ensembl.org"
        if gene_symbol in HG19_ALIASES:
            gene_symbol = HG19_ALIASES[gene_symbol]

    url = f"{server}/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1&content-type=application/json"

    try:
        data = _get(url)
        if not data:
            return None

        canonical = None
        longest = None
        max_len = -1

        transcripts = data.get('Transcript', [])
        for t in transcripts:
            if t.get('is_canonical') == 1:
                canonical = t
                break
            length = t.get('end', 0) - t.get('start', 0)
            if length > max_len:
                max_len = length
                longest = t

        target = canonical if canonical else longest

        if target:
            return target['id']

    except Exception as e:
        print(f"Error lookup canonical on {build} for {gene_symbol}: {e}")

    return None

def fetch_ensembl_exons(transcript_id_version, build='hg38', gene_symbol=None):
    """
    Fetch exon structure for a given transcript ID from Ensembl REST API.
    transcript_id_version: e.g. ENST00000357654.9
    build: 'hg38' or 'hg19'
    gene_symbol: Optional, used for fallback if ID lookup fails.

    Returns a dict with transcript structure, or None on failure.
    """
    transcript_id = transcript_id_version.split('.')[0]

    if build == 'hg38':
        server = "https://rest.ensembl.org"
    else:
        server = "https://grch37.rest.ensembl.org"

    url = f"{server}/lookup/id/{transcript_id}?expand=1&content-type=application/json"

    fetched_data = None

    # Retry logic for INITIAL ID
    for attempt in range(5):
        try:
            fetched_data = _get(url)
            break
        except urllib.error.HTTPError as e:
            if e.code == 429:
                # Rate limit: longer backoff
                time.sleep(2 + attempt * 2)
            elif e.code >= 500:
                # Server error: short backoff
                time.sleep(1 + attempt)
            elif e.code in (400, 404):
                # Not found — no point retrying
                break
            else:
                time.sleep(1)
        except Exception as e:
            print(f"Request failed (attempt {attempt+1}): {e}")
            time.sleep(1 + attempt)

    # FALLBACK: If failed and we have gene symbol, try to find the GENE on this build
    if not fetched_data and gene_symbol:
        alt_id = find_canonical_on_build(gene_symbol, build)
        if alt_id and alt_id != transcript_id:
            return fetch_ensembl_exons(alt_id, build=build, gene_symbol=None)

    if not fetched_data:
        return None

    data = fetched_data

    # If no Translation (CDS), try fallback to canonical
    if 'Translation' not in data and gene_symbol:
        alt_id = find_canonical_on_build(gene_symbol, build)
        if alt_id and alt_id != transcript_id:
             return fetch_ensembl_exons(alt_id, build=build, gene_symbol=None)

    # Extract Exons
    exons = []
    if 'Exon' in data:
        for ex in data['Exon']:
            exons.append({
                'id': ex['id'],
                'start': ex['start'],
                'end': ex['end'],
                'strand': ex['strand'],
                'seq_region_name': ex['seq_region_name'],
                'object_type': 'Exon'
            })

    cds_start = None
    cds_end = None
    if 'Translation' in data:
        cds_start = data['Translation']['start']
        cds_end = data['Translation']['end']

    return {
        'transcript_id': transcript_id,
        'exons': exons,
        'cds_start': cds_start,
        'cds_end': cds_end,
        'chr': data['seq_region_name'],
        'strand': data['strand']
    }

def generate_bed_records(transcript_info, include_5utr=False, include_3utr=False, include_intron=False, flank_5_prime=0, flank_3_prime=0, exon_filter=None, intron_filter=None):
    """
    Generate BED lines from transcript info.
    transcript_info: output from fetch_ensembl_exons
    exon_filter:   set of ints (1-based transcript-order exon numbers) to include.
                   None = include all exons.
    intron_filter: set of ints (1-based transcript-order intron numbers) to include.
                   None = include all introns (subject to include_intron flag).
                   Intron N lies between exon N and exon N+1.
    Flanks are always emitted regardless of filters (gene-level regions).
    """
    records = []
    chrom = transcript_info['chr']
    strand = transcript_info['strand']

    if not chrom.startswith('chr'):
        chrom = f"chr{chrom}"

    exons = sorted(transcript_info['exons'], key=lambda x: x['start'])
    cds_start = transcript_info['cds_start']
    cds_end = transcript_info['cds_end']

    min_start = exons[0]['start']
    max_end = exons[-1]['end']

    # 5' Flank (Promoter — upstream of TSS)
    if flank_5_prime > 0:
        if strand == 1:
            f5_start_1base = max(1, min_start - flank_5_prime)
            f5_end_1base = min_start - 1
            if f5_end_1base >= f5_start_1base:
                records.append((chrom, f5_start_1base - 1, f5_end_1base, "promoter_5prime"))
        else:
            f5_start_1base = max_end + 1
            f5_end_1base = max_end + flank_5_prime
            records.append((chrom, f5_start_1base - 1, f5_end_1base, "promoter_5prime"))

    # 3' Flank (Downstream — past TES)
    if flank_3_prime > 0:
        if strand == 1:
            f3_start_1base = max_end + 1
            f3_end_1base = max_end + flank_3_prime
            records.append((chrom, f3_start_1base - 1, f3_end_1base, "downstream_3prime"))
        else:
            f3_start_1base = max(1, min_start - flank_3_prime)
            f3_end_1base = min_start - 1
            if f3_end_1base >= f3_start_1base:
                records.append((chrom, f3_start_1base - 1, f3_end_1base, "downstream_3prime"))

    total_exons = len(exons)

    final_regions = []

    for i, ex in enumerate(exons):
        start = ex['start']
        end = ex['end']  # 1-based Ensembl

        if strand == 1:
            exon_num = i + 1
        else:
            exon_num = total_exons - i

        region_base_name = f"exon{exon_num}"

        # Exon filter: skip if not in the requested set
        if exon_filter and exon_num not in exon_filter:
            continue

        # Non-coding transcript: emit whole exon
        if cds_start is None:
            final_regions.append((start - 1, end, region_base_name))
            continue

        # Part 1: Pre-CDS
        if start < cds_start:
            u_end = min(end, cds_start - 1)
            if u_end >= start:
                type_ = '5UTR' if strand == 1 else '3UTR'
                if (type_ == '5UTR' and include_5utr) or (type_ == '3UTR' and include_3utr):
                    final_regions.append((start - 1, u_end, f"{region_base_name}_{type_}"))

        # Part 2: CDS
        c_start = max(start, cds_start)
        c_end = min(end, cds_end)
        if c_end >= c_start:
            final_regions.append((c_start - 1, c_end, f"{region_base_name}_CDS"))

        # Part 3: Post-CDS
        if end > cds_end:
            u_start = max(start, cds_end + 1)
            if end >= u_start:
                type_ = '3UTR' if strand == 1 else '5UTR'
                if (type_ == '5UTR' and include_5utr) or (type_ == '3UTR' and include_3utr):
                    final_regions.append((u_start - 1, end, f"{region_base_name}_{type_}"))

    # Introns
    if include_intron:
        for i in range(len(exons) - 1):
            bed_intron_start = exons[i]['end']          # 0-based: first intronic position
            bed_intron_end = exons[i+1]['start'] - 1    # 0-based exclusive: one past last intronic position

            if strand == 1:
                intron_num = i + 1
                left_exon_num  = i + 1
                right_exon_num = i + 2
            else:
                intron_num = total_exons - i - 1
                left_exon_num  = total_exons - i
                right_exon_num = total_exons - i - 1

            if intron_filter and intron_num not in intron_filter:
                continue
            if exon_filter and not (left_exon_num in exon_filter and right_exon_num in exon_filter):
                continue

            if bed_intron_end > bed_intron_start:
                records.append((chrom, bed_intron_start, bed_intron_end, f'intron{intron_num}'))

    # Add exon regions
    for s, e, type_ in final_regions:
        records.append((chrom, s, e, type_))

    return records
