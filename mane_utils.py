import pandas as pd
import requests
import time
import os

API_TIMEOUT = 15  # seconds — prevents indefinite hangs on stalled Ensembl connections

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
    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    
    results = {}
    
    # Process in chunks of 50 (Ensembl POST limit)
    chunk_size = 50
    for i in range(0, len(transcript_ids), chunk_size):
        chunk = transcript_ids[i:i+chunk_size]
        
        # Clean IDs: remove versions for lookup
        chunk_clean = [tid.split('.')[0] for tid in chunk]
        
        # Map cleaned ID back to original for result storage
        # CAUTION: If multiple original IDs map to same clean ID (rare in this context), we must handle it.
        # But here we just assume 1-to-1 or standard usage.
        
        # Ensembl /sequence/id POST: ids + type must be top-level keys.
        # The endpoint returns a JSON list: [{id, seq, ...}, ...]
        payload = {
            "ids": chunk_clean,
            "type": "cdna",
        }

        try:
            r = requests.post(server + ext, headers=headers, json=payload, timeout=API_TIMEOUT)
            if not r.ok:
                print(f"Batch GC fetch failed (HTTP {r.status_code}) for chunk starting at index {i}: {r.text[:200]}")
                for tid in chunk:
                    results[tid] = None
                continue

            data = r.json()

            # Ensembl returns a list of objects: [{id: <unversioned_id>, seq: <sequence>}, ...]
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

        except Exception as e:
            print(f"Error in batch GC fetch: {e}")
            for tid in chunk:
                results[tid] = None
                
        # Be nice to API
        time.sleep(0.1)
        
    return results
MANE_SUMMARY_FILE = "MANE.summary.txt"

def load_mane_data():
    """Bypassed if file doesn't exist, returns None."""
    if not os.path.exists(MANE_SUMMARY_FILE):
        return None
    # Parse the MANE summary file
    df = pd.read_csv(MANE_SUMMARY_FILE, sep="\t")
    return df

def get_mane_transcripts(gene_list, mane_df):
    """
    Filter MANE data for the given genes with prioritization:
    1. MANE Select
    2. MANE Plus Clinical (only if Select is missing)
    
    Returns a DataFrame with relevant columns.
    """
    # Create valid set for lookup
    valid_genes = set(g.upper().strip() for g in gene_list)
    
    # Filter DF for these genes
    mask = mane_df['symbol'].str.upper().isin(valid_genes)
    subset = mane_df[mask].copy()
    
    # Sort to prioritize MANE Select
    subset['MANE_status'] = pd.Categorical(subset['MANE_status'], categories=['MANE Select', 'MANE Plus Clinical'], ordered=True)
    subset = subset.sort_values('MANE_status')
    
    # Drop duplicates keeping the first
    result = subset.drop_duplicates(subset=['symbol'], keep='first')
    
    # Normalize result
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
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # Special manual mappings for loci or ambiguous symbols
    MANUAL_ALIASES = {
        'IGH': 'IGHG1', # Representative Constant Region (Gamma 1)
        'IGK': 'IGKC',  # Kappa Constant Region
        'IGL': 'IGLC1', # Lambda Constant Region 1 (Representative)
    }
    
    records = []
    
    # Helper to resolve aliases for genes not found in strict lookup
    def resolve_alias(symbol, server):
        try:
            # 1. Look for cross-references (often finds old symbols)
            xref_ext = f"/xrefs/symbol/homo_sapiens/{symbol}?object_type=gene"
            # print(f"DEBUG: resolving alias for {symbol} at {server+xref_ext}")
            r = requests.get(server + xref_ext, headers=headers, timeout=API_TIMEOUT)
            if not r.ok:
                print(f"Alias lookup failed for {symbol}: {r.status_code}")
                return None
            if not r.json():
                print(f"Alias lookup empty for {symbol}")
                return None

            # Use the first gene match
            gene_id = r.json()[0]['id']

            # 2. Lookup full details by ID (gives us the new Symbol and Transcripts)
            lookup_ext = f"/lookup/id/{gene_id}?expand=1"
            r2 = requests.get(server + lookup_ext, headers=headers, timeout=API_TIMEOUT)
            if r2.ok:
                return r2.json()
        except Exception as e:
            print(f"Error resolving alias for {symbol}: {e}")
        return None

    # Helper to chunk list
    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    for batch in chunks(list(gene_list), 200):
        try:
            # POST request
            r = requests.post(server+ext, headers=headers, json={"symbols": batch, "expand": 1}, timeout=API_TIMEOUT)
            
            if not r.ok:
                print(f"Bulk lookup failed: {r.status_code}")
                continue
                
            data = r.json()
            
            for gene in batch:
                info = data.get(gene)
                
                # If strict lookup failed, try alias resolution
                if not info:
                    if gene in MANUAL_ALIASES:
                        # Direct lookup of manual override
                        target_alias = MANUAL_ALIASES[gene]
                        r_lookup = requests.get(server+f"/lookup/symbol/homo_sapiens/{target_alias}?expand=1", headers=headers, timeout=API_TIMEOUT)
                        if r_lookup.ok:
                            info = r_lookup.json()
                            
                    if not info:
                        info = resolve_alias(gene, server)
                        if not info:
                            continue
                    
                # Find canonical
                canonical_transcript = None
                transcripts = info.get('Transcript', [])
                
                # First pass: Check is_canonical flag
                for t in transcripts:
                    if t.get('is_canonical') == 1:
                        canonical_transcript = t
                        break
                
                # Fallback: if no is_canonical (rare), take longest? or first?
                if not canonical_transcript and transcripts:
                    canonical_transcript = transcripts[0]
                    
                if canonical_transcript:
                    records.append({
                        'original_input': gene, # Track original for filtering/reporting
                        'symbol': info.get('display_name', gene),
                        'Ensembl_nuc': canonical_transcript['id'],
                        'RefSeq_nuc': '-', # Hard to get easily without more lookups
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

# Common HGNC symbol changes (New -> Old on hg19)
HG19_ALIASES = {
    'H2BC6': 'HIST1H2BC',
    'H3C2': 'HIST1H3B', 
    'H2AC11': 'HIST1H2AC',
    'H2AC6': 'HIST1H2AB',
    # Add others as discovered or implement a more robust alias lookup if needed
}

def find_canonical_on_build(gene_symbol, build='hg38'):
    """
    Find the canonical transcript ID for a given gene on a specific build.
    """
    if build == 'hg38':
        server = "https://rest.ensembl.org"
    else:
        server = "https://grch37.rest.ensembl.org"
        # Check aliases for hg19
        if gene_symbol in HG19_ALIASES:
            # print(f"Using alias {HG19_ALIASES[gene_symbol]} for {gene_symbol} on hg19")
            gene_symbol = HG19_ALIASES[gene_symbol]
        
    ext = f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1;content-type=application/json"

    try:
        r = requests.get(server+ext, headers={"Content-Type": "application/json"}, timeout=API_TIMEOUT)
        if not r.ok:
            return None
        data = r.json()
        
        # Find canonical
        canonical = None
        longest = None
        max_len = -1
        
        transcripts = data.get('Transcript', [])
        for t in transcripts:
            if t.get('is_canonical') == 1:
                canonical = t
                break
            # Track longest just in case
            length = t.get('end') - t.get('start')
            if length > max_len:
                max_len = length
                longest = t
                
        # If no explicit canonical tag (common in older builds?), use longest
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
    
    Returns a list of dicts: [{'start': 100, 'end': 200, 'rank': 1, 'id': 'ENSE...'}]
    """
    
    # Strip version for safety when query cross-builds
    # Often for hg19 we need the generic ID or a different ID entirely
    transcript_id = transcript_id_version.split('.')[0]
    
    if build == 'hg38':
        server = "https://rest.ensembl.org"
    else:
        server = "https://grch37.rest.ensembl.org"
        
    ext = f"/lookup/id/{transcript_id}?expand=1;content-type=application/json"
    
    fetched_data = None
    
    # Retry logic for INITIAL ID
    for attempt in range(2):
        try:
            r = requests.get(server+ext, headers={"Content-Type": "application/json"}, timeout=API_TIMEOUT)
            if r.ok:
                fetched_data = r.json()
                break
            elif r.status_code == 429:
                time.sleep(1)
        except Exception:
            pass
            
    # FALLBACK: If failed and we have gene symbol, try to find the GENE on this build
    if not fetched_data and gene_symbol:
        # print(f"Fallback: fetching canonical for {gene_symbol} on {build}")
        alt_id = find_canonical_on_build(gene_symbol, build)
        if alt_id and alt_id != transcript_id:
            return fetch_ensembl_exons(alt_id, build=build, gene_symbol=None)

    if not fetched_data:
        return None
        
    data = fetched_data
    
    # NEW CHECK: Does it have Translation (CDS)?
    # If not, and we haven't already fallen back (based on ID), try fallback to canonical
    # This fixes cases like MAP3K14 where the MANE ID exists on hg19 but is non-coding/incomplete.
    if 'Translation' not in data and gene_symbol:
        # Avoid infinite recursion if we already are looking up canonical
        # We only try fallback if the current transcript seems "broken" (no CDS)
        # But we must ensure the canonical one IS different.
        alt_id = find_canonical_on_build(gene_symbol, build)
        if alt_id and alt_id != transcript_id:
             # print(f"Fallback (Missing CDS): fetching canonical {alt_id} for {gene_symbol}")
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

def generate_bed_records(transcript_info, include_5utr=False, include_3utr=False, include_intron=False, flank_5_prime=0, flank_3_prime=0):
    """
    Generate BED lines from transcript info.
    transcript_info: output from fetch_ensembl_exons
    """
    records = []
    chrom = transcript_info['chr']
    strand = transcript_info['strand']
    
    if not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
        
    exons = sorted(transcript_info['exons'], key=lambda x: x['start'])
    cds_start = transcript_info['cds_start']
    cds_end = transcript_info['cds_end']
    
    # Calculate Flanks (Strand Aware)
    # The 'exons' list gives us the full transcribed span [min_start, max_end]
    min_start = exons[0]['start']
    max_end = exons[-1]['end']
    
    # Flanks use 1-based Ensembl coordinates internally and convert to BED (0-based start, exclusive end)
    # via (1base_start - 1, 1base_end).

    # 5' Flank (Promoter — upstream of TSS)
    if flank_5_prime > 0:
        if strand == 1:
            # (+) strand: TSS is min_start. Upstream = lower coordinates.
            # 1-based region: [min_start - flank, min_start - 1]
            f5_start_1base = max(1, min_start - flank_5_prime)
            f5_end_1base = min_start - 1
            if f5_end_1base >= f5_start_1base:  # guard against TSS at chr start
                records.append((chrom, f5_start_1base - 1, f5_end_1base, "promoter_5prime"))
        else:
            # (-) strand: TSS is max_end. Upstream = higher coordinates.
            # 1-based region: [max_end + 1, max_end + flank]
            f5_start_1base = max_end + 1
            f5_end_1base = max_end + flank_5_prime
            records.append((chrom, f5_start_1base - 1, f5_end_1base, "promoter_5prime"))

    # 3' Flank (Downstream — past TES)
    if flank_3_prime > 0:
        if strand == 1:
            # (+) strand: TES is max_end. Downstream = higher coordinates.
            # 1-based region: [max_end + 1, max_end + flank]
            f3_start_1base = max_end + 1
            f3_end_1base = max_end + flank_3_prime
            records.append((chrom, f3_start_1base - 1, f3_end_1base, "downstream_3prime"))
        else:
            # (-) strand: TES is min_start. Downstream = lower coordinates.
            # 1-based region: [min_start - flank, min_start - 1]
            f3_start_1base = max(1, min_start - flank_3_prime)
            f3_end_1base = min_start - 1
            if f3_end_1base >= f3_start_1base:  # guard against TES at chr start
                records.append((chrom, f3_start_1base - 1, f3_end_1base, "downstream_3prime"))
    
    # Numbering Logic
    total_exons = len(exons)
    
    final_regions = []
    
    # We iterate genomic order (asc start)
    for i, ex in enumerate(exons):
        start = ex['start']
        end = ex['end'] # 1-based Ensembl
        
        # Determine Exon Number
        if strand == 1:
            exon_num = i + 1
        else:
            exon_num = total_exons - i
            
        region_base_name = f"exon{exon_num}"
        
        # Split logic into 5'UTR, CDS, 3'UTR parts
        # 1. Identify UTR content
        
        # If no CDS defined (non-coding), treat whole exon as a generic exon record.
        if cds_start is None:
            final_regions.append((start - 1, end, region_base_name))
            continue

        # Split Exon into potentially 3 parts: [Start, CDS_Start-1], [CDS_Start, CDS_End], [CDS_End+1, End]
        
        # Part 1: Pre-CDS (Genomic 5' if (+), Genomic 3' if (-))
        if start < cds_start:
            u_end = min(end, cds_start - 1)
            if u_end >= start:
                # This is UTR
                # Logic:
                # Strand (+): Start < CDS_Start -> 5' UTR
                # Strand (-): Start < CDS_Start -> this is physically 3' of the gene (higher coord is 5')
                # Wait:
                # (+) 5'----Start(100)....CDS(200)....End(300)----3'
                #     100-199 is 5'UTR.
                # (-) 3'----Start(100)....CDS(200)....End(300)----5'
                #     100-199 is 3'UTR (because gene runs 300->100).
                
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
                # Strand (+): > CDS_End -> 3' UTR
                # Strand (-): > CDS_End -> 5' UTR (lower coord is 3')
                # Wait:
                # (-) CDS(200)...End(300) -> 201-300 is 5' of CDS start (gene body), wait.
                # (-) Gene coords: 1000 (Start) to 500 (End). 
                # Ensembl always gives Start < End (Genomic).
                # ex['start'] = 500, ex['end'] = 600.
                # Translation['start'] = 550, ['end'] = 580.
                # 500-549: Genomic 'left' of CDS. For (-) strand, this is 3' side (Sequence is 600->500).
                # So 600 is 5', 500 is 3'.
                # 581-600 (Genomic Right): This is 5' side of CDS. -> 5' UTR.
                # 500-549 (Genomic Left): This is 3' side of CDS. -> 3' UTR.
                
                # Check previous logic:
                # Pre-CDS (Start < CDS_Start): 
                # (+) -> 5' UTR
                # (-) -> 3' UTR
                # Correct.
                
                # Post-CDS (End > CDS_End):
                # (+) -> 3' UTR
                # (-) -> 5' UTR
                
                type_ = '3UTR' if strand == 1 else '5UTR'
                if (type_ == '5UTR' and include_5utr) or (type_ == '3UTR' and include_3utr):
                    final_regions.append((u_start - 1, end, f"{region_base_name}_{type_}"))

    # Introns
    # Ensembl coords are 1-based inclusive.
    # Intron spans from (exon_i_end + 1) to (next_exon_start - 1), both 1-based inclusive.
    # BED 0-based half-open: start = exon_i_end (= 1-based intron start - 1),
    #                         end   = next_exon_start - 1 (= 1-based last intron base, BED end is exclusive so this is correct).
    if include_intron:
        for i in range(len(exons) - 1):
            bed_intron_start = exons[i]['end']          # 0-based: first intronic position
            bed_intron_end = exons[i+1]['start'] - 1    # 0-based exclusive: one past last intronic position

            if strand == 1:
                intron_num = i + 1
            else:
                intron_num = total_exons - i - 1

            if bed_intron_end > bed_intron_start:
                records.append((chrom, bed_intron_start, bed_intron_end, f'intron{intron_num}'))
                 
    # Add final regions
    for s, e, type_ in final_regions:
        records.append((chrom, s, e, type_))
        
    return records
