import streamlit as st
import pandas as pd
import mane_utils
import liftover_utils
import io
import re

def parse_exon_filter(value):
    """
    Parse an exon filter string into a set of ints.
    Returns:
        None          — if value is blank / "all" (include everything)
        set of ints   — the requested exon numbers (1-based)
    Raises ValueError with a descriptive message on bad syntax.

    Examples:
        "all"   -> None
        ""      -> None
        "3"     -> {3}
        "1,3,5" -> {1, 3, 5}
        "1-3,7" -> {1, 2, 3, 7}
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if str(value).strip().lower() in ("all", "", "nan"):
        return None
    result = set()
    for part in str(value).replace(" ", "").split(","):
        if "-" in part:
            bounds = part.split("-")
            if len(bounds) != 2 or not all(b.isdigit() for b in bounds):
                raise ValueError(f"Invalid range '{part}' — use format like '5-9'")
            lo, hi = int(bounds[0]), int(bounds[1])
            if lo > hi:
                raise ValueError(f"Invalid range '{part}' — start must be ≤ end")
            result.update(range(lo, hi + 1))
        elif part.isdigit():
            result.add(int(part))
        else:
            raise ValueError(f"Invalid token '{part}' — use integers, commas and ranges (e.g. '1,3,5-7')")
    if not result:
        raise ValueError("Exon filter is empty after parsing")
    return result


st.set_page_config(
    page_title="Gene Panel Designer",
    page_icon="🧬",
    layout="wide",
)

# ── Global CSS injected once ──────────────────────────────────────────────────
st.markdown("""
<style>
/* Step cards */
.step-card {
    background: #162032;
    border: 1px solid #1f3050;
    border-radius: 10px;
    padding: 18px 20px;
    text-align: center;
    height: 100%;
}
.step-number {
    display: inline-block;
    background: #00c9a7;
    color: #0f1923;
    font-weight: 700;
    font-size: 0.85rem;
    border-radius: 50%;
    width: 28px;
    height: 28px;
    line-height: 28px;
    margin-bottom: 8px;
}
.step-title {
    font-size: 0.95rem;
    font-weight: 600;
    color: #e8edf2;
    margin: 4px 0 2px 0;
}
.step-desc {
    font-size: 0.78rem;
    color: #8fa3b8;
}
/* Info card (liftover) */
.info-card {
    background: #162032;
    border-left: 4px solid #00c9a7;
    border-radius: 6px;
    padding: 16px 20px;
    margin-bottom: 16px;
    color: #b0c4d8;
    font-size: 0.9rem;
    line-height: 1.6;
}
/* Prominent download area */
.download-summary {
    background: #0d2b22;
    border: 1px solid #00c9a7;
    border-radius: 10px;
    padding: 20px 24px;
    margin-top: 12px;
    text-align: center;
}
.download-summary h3 {
    color: #00c9a7;
    margin: 0 0 4px 0;
    font-size: 1.15rem;
}
.download-summary p {
    color: #8fa3b8;
    font-size: 0.85rem;
    margin: 0 0 14px 0;
}
/* Subtle tag badges */
.badge {
    display: inline-block;
    background: #1f3050;
    color: #00c9a7;
    border-radius: 4px;
    padding: 2px 8px;
    font-size: 0.75rem;
    font-weight: 600;
    margin-right: 4px;
}
</style>
""", unsafe_allow_html=True)

# ── Hero header ───────────────────────────────────────────────────────────────
st.markdown("""
<div style="padding: 10px 0 6px 0;">
    <span style="font-size:2.1rem; font-weight:700; color:#e8edf2; letter-spacing:-0.5px;">
        🧬 Gene Panel Designer
    </span><br/>
    <span style="font-size:1rem; color:#8fa3b8; line-height:1.8;">
        Design and export clinical capture panels —
        <span style="color:#00c9a7; font-weight:600;">from gene list to ready-to-order BED file</span>
        in minutes.
    </span>
</div>
""", unsafe_allow_html=True)
st.divider()

# ── Load MANE Data ────────────────────────────────────────────────────────────
@st.cache_data
def load_data():
    return mane_utils.load_mane_data()

mane_df = load_data()

if mane_df is None:
    st.error("MANE Summary file not found. Please ensure 'MANE.summary.txt' is in the application directory.")
    st.stop()

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab1, tab2 = st.tabs(["🧬  Panel Designer", "🔄  LiftOver (hg19 → hg38)"])

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 1 — PANEL DESIGNER
# ═══════════════════════════════════════════════════════════════════════════════
with tab1:

    # ── Sidebar ───────────────────────────────────────────────────────────────
    st.sidebar.markdown("""
    <div style="padding:4px 0 10px 0;">
        <span style="font-size:1.1rem; font-weight:700; color:#00c9a7;">⚙️ Panel Settings</span><br/>
        <span style="font-size:0.78rem; color:#8fa3b8;">
            Configure defaults here.<br/>Per-gene overrides are set in the table below.
        </span>
    </div>
    """, unsafe_allow_html=True)

    genome_build = st.sidebar.radio(
        "Genome Build",
        ["hg19", "hg38"],
        index=1,
        help="Reference genome for BED coordinate output."
    )

    st.sidebar.divider()

    id_format = st.sidebar.radio(
        "Transcript ID Format",
        ["Ensembl (ENST)", "RefSeq (NM)"],
        help="Identifier used in the BED name column. RefSeq falls back to Ensembl when unavailable."
    )

    st.sidebar.divider()

    st.sidebar.markdown(
        "<span style='font-size:0.82rem; font-weight:600; color:#b0c4d8;'>TRANSCRIPT CONTENT</span>",
        unsafe_allow_html=True
    )
    default_5utr   = st.sidebar.checkbox("Include 5′ UTR", value=False, help="Can be toggled per gene in the table.")
    default_3utr   = st.sidebar.checkbox("Include 3′ UTR", value=False, help="Can be toggled per gene in the table.")
    default_intron = st.sidebar.checkbox("Include Introns", value=False, help="Can be toggled per gene in the table.")

    fetch_gc = st.sidebar.checkbox(
        "Calculate GC Content",
        value=True,
        help="Fetch cDNA GC % from Ensembl for each transcript. Adds ~10–30 s for large panels. "
             "Disable if you don't need it or the API is slow."
    )

    st.sidebar.divider()

    padding = st.sidebar.number_input(
        "Padding (bp)",
        min_value=0,
        value=20,
        help="Bases added to each side of every region. 20 bp is standard for clinical panels."
    )

    st.sidebar.divider()

    if st.sidebar.button(
        "🔄 Reset & Re-fetch",
        help="Clears the cached mapping and re-fetches transcripts from Ensembl. "
             "Use this if transcript data looks stale or if you changed the genome build.",
        use_container_width=True,
    ):
        for key in ['mapped_data', 'last_file_hash']:
            st.session_state.pop(key, None)
        st.rerun()

    # ── Step banner ───────────────────────────────────────────────────────────
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">1</div>
            <div class="step-title">Upload Gene List</div>
            <div class="step-desc">CSV, XLSX or plain TXT — one gene symbol per row</div>
        </div>""", unsafe_allow_html=True)
    with c2:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">2</div>
            <div class="step-title">Review & Customise</div>
            <div class="step-desc">Verify MANE transcripts · toggle UTRs, introns, promoter flanks per gene</div>
        </div>""", unsafe_allow_html=True)
    with c3:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">3</div>
            <div class="step-title">Generate &amp; Download</div>
            <div class="step-desc">Export a padded, strand-aware BED file ready to submit to your capture vendor</div>
        </div>""", unsafe_allow_html=True)

    st.markdown("<div style='margin-top:24px;'></div>", unsafe_allow_html=True)

    # ── File upload ───────────────────────────────────────────────────────────
    uploaded_file = st.file_uploader(
        "Upload your gene list",
        type=["csv", "xlsx", "txt"],
        key="designer_upload",
        help="One gene symbol per row. Aliases in parentheses (e.g. CD274 (PD-L1)) are stripped automatically."
    )

    if uploaded_file:
        # ── Parse ─────────────────────────────────────────────────────────────
        genes = []
        try:
            if uploaded_file.name.endswith('.csv'):
                df_raw = pd.read_csv(uploaded_file, header=None)
                genes = df_raw[0].astype(str).tolist()
            elif uploaded_file.name.endswith('.xlsx'):
                df_raw = pd.read_excel(uploaded_file, header=None)
                genes = df_raw[0].astype(str).tolist()
            else:
                content = uploaded_file.read().decode()
                genes = [line.strip() for line in content.splitlines() if line.strip()]

            # Clean: strip aliases in parentheses e.g. "CD274 (PD-L1)" → "CD274"
            cleaned = []
            for g in genes:
                g_clean = re.sub(r'\(.*?\)', '', g).strip()
                if g_clean:
                    cleaned.append(g_clean)
            genes = cleaned
        except Exception as e:
            st.error(f"Could not parse file: {e}")
            st.stop()

        st.success(f"✅ {len(genes)} gene symbols loaded.")

        # ── Map to MANE (cached by file content hash) ─────────────────────────
        # SCHEMA_VERSION must be bumped whenever new columns are added to mapped_df,
        # so stale session state from previous app versions is automatically invalidated.
        SCHEMA_VERSION = 3
        file_hash = hash(uploaded_file.getvalue())
        cache_key = (file_hash, SCHEMA_VERSION)
        if 'mapped_data' not in st.session_state or st.session_state.get('last_file_hash') != cache_key:

            mapped_df = mane_utils.get_mane_transcripts(genes, mane_df)

            # Fallback: Ensembl Canonical for genes absent from MANE
            found_symbols  = set(mapped_df['symbol'].str.upper())
            input_symbols  = set(g.upper() for g in genes)
            missing_genes  = input_symbols - found_symbols

            if missing_genes:
                status_ph = st.info(f"⏳ Fetching Ensembl Canonical for {len(missing_genes)} non-MANE gene(s)…")
                fallback_df = mane_utils.get_ensembl_canonical(missing_genes, build=genome_build)
                if not fallback_df.empty:
                    mapped_df = pd.concat([mapped_df, fallback_df], ignore_index=True)
                    status_ph.empty()

                    if 'original_input' in fallback_df.columns:
                        aliases = fallback_df[
                            fallback_df['original_input'].str.upper() != fallback_df['symbol'].str.upper()
                        ]
                        if not aliases.empty:
                            with st.expander(f"ℹ️ {len(aliases)} gene alias(es) resolved", expanded=True):
                                st.table(
                                    aliases[['original_input', 'symbol']].rename(
                                        columns={'original_input': 'Input', 'symbol': 'Official Symbol'}
                                    )
                                )

            # GC content — optional, always fetched from hg38 (MANE IDs are GRCh38-defined).
            if fetch_gc:
                gc_ph = st.info("⏳ Calculating GC content…")
                gc_map = mane_utils.fetch_transcript_gc_batch(mapped_df['Ensembl_nuc'].tolist(), build='hg38')
                mapped_df['Mean GC %'] = mapped_df['Ensembl_nuc'].map(gc_map)
                gc_ph.empty()

                n_gc_missing = mapped_df['Mean GC %'].isna().sum()
                if n_gc_missing == len(mapped_df):
                    st.warning("⚠️ GC content could not be fetched from Ensembl (API may be temporarily unavailable). "
                               "GC % will show as N/A — this does not affect BED generation.")
                elif n_gc_missing > 0:
                    st.caption(f"ℹ️ GC content unavailable for {n_gc_missing} transcript(s) — shown as N/A.")
            else:
                mapped_df['Mean GC %'] = None

            # Config columns
            mapped_df["Include 5' UTR"]  = default_5utr
            mapped_df["Include 3' UTR"]  = default_3utr
            mapped_df["Exon Filter"]     = "all"   # "all" or comma/range list e.g. "1,3,5-7"
            mapped_df["Include Intron"]  = default_intron
            mapped_df["Intron Filter"]   = "all"   # "all" or comma/range list e.g. "1,3,5-7"
            mapped_df["5' Flank"]        = 0
            mapped_df["3' Flank"]        = 0

            st.session_state['mapped_data']   = mapped_df
            st.session_state['last_file_hash'] = cache_key

        mapped_df = st.session_state['mapped_data']

        # ── Summary metric cards ───────────────────────────────────────────────
        n_total      = len(mapped_df)
        n_mane_sel   = (mapped_df['MANE_status'] == 'MANE Select').sum()
        n_other      = n_total - n_mane_sel
        n_not_found  = len(set(g.upper() for g in genes) - set(mapped_df['symbol'].str.upper()))

        mc1, mc2, mc3, mc4 = st.columns(4)
        mc1.metric("Genes submitted",   len(genes))
        mc2.metric("Mapped",            n_total)
        mc3.metric("MANE Select",       int(n_mane_sel))
        mc4.metric("Non-MANE / Alt",    int(n_other),
                   delta=f"{n_not_found} not found" if n_not_found else None,
                   delta_color="inverse" if n_not_found else "off")

        # ── Non-MANE transcript warning ────────────────────────────────────────
        non_mane = mapped_df[mapped_df['MANE_status'] != 'MANE Select']
        if not non_mane.empty:
            with st.expander(f"⚠️ {len(non_mane)} gene(s) using non-MANE Select transcripts", expanded=False):
                st.caption("These genes did not have a MANE Select transcript. An alternative was chosen automatically.")
                st.table(
                    non_mane[['symbol', 'Ensembl_nuc', 'MANE_status']].rename(
                        columns={'symbol': 'Gene', 'Ensembl_nuc': 'Transcript ID', 'MANE_status': 'Source'}
                    )
                )
        else:
            msg = ("✅ All genes mapped to MANE Select transcripts."
                   if genome_build == 'hg38'
                   else "✅ All genes mapped to MANE Select. (hg19 runtime switches may occur during BED generation.)")
            st.success(msg)

        st.divider()

        # ── Per-gene settings table ────────────────────────────────────────────
        st.markdown("### ⚙️ Customise Per-Gene Settings")
        st.caption("Edit UTR/intron toggles, exon/intron filters, and flanks per gene. "
                   "Filters accept blank/'all', numbers ('2,3'), or ranges ('19-21'). All other columns are read-only.")

        edited_df = st.data_editor(
            mapped_df,
            column_config={
                "symbol":         st.column_config.TextColumn("Gene",        disabled=True),
                "Ensembl_nuc":    st.column_config.TextColumn("Transcript",  disabled=True),
                "RefSeq_nuc":     st.column_config.TextColumn("RefSeq ID",   disabled=True),
                "MANE_status":    st.column_config.TextColumn("Type",        disabled=True),
                "Ensembl_Gene":   st.column_config.TextColumn("Gene ID",     disabled=True),
                "GRCh38_chr":     None,  # hidden
                "Mean GC %":      (
                                      st.column_config.NumberColumn(
                                          "GC %",
                                          format="%.1f",
                                          min_value=0, max_value=100,
                                          help="Mean GC content of the cDNA sequence (informational). N/A = could not fetch from Ensembl.",
                                          disabled=True,
                                      ) if fetch_gc else None   # hide column when GC fetch is disabled
                                  ),
                "Include 5' UTR": st.column_config.CheckboxColumn("5′ UTR",   help="Include 5′ UTR regions"),
                "Include 3' UTR": st.column_config.CheckboxColumn("3′ UTR",   help="Include 3′ UTR regions"),
                "Exon Filter":    st.column_config.TextColumn(
                                      "Exon Filter",
                                      help=(
                                          "Which exons to include. Leave blank or type 'all' for every exon. "
                                          "Otherwise enter exon numbers (transcript order, 1-based): "
                                          "e.g. '2,3' or '19-21' or '1,5-7,10'."
                                      ),
                                  ),
                "Include Intron": st.column_config.CheckboxColumn("Introns",  help="Enable intronic regions for this gene. Use Intron Filter to select specific introns."),
                "Intron Filter":  st.column_config.TextColumn(
                                      "Intron Filter",
                                      help=(
                                          "Which introns to include (only used when Introns is ticked). "
                                          "Leave blank or type 'all' for every intron. "
                                          "Otherwise enter intron numbers (transcript order, 1-based): "
                                          "e.g. '1,2' or '4-6' or '1,3,5-7'. "
                                          "Intron N lies between exon N and exon N+1."
                                      ),
                                  ),
                "5' Flank":       st.column_config.NumberColumn(
                                      "5′ Flank (bp)",
                                      help="Upstream promoter padding added upstream of TSS (e.g. 2000 for MYC).",
                                      min_value=0
                                  ),
                "3' Flank":       st.column_config.NumberColumn(
                                      "3′ Flank (bp)",
                                      help="Downstream padding added past the transcript end.",
                                      min_value=0
                                  ),
            },
            disabled=["symbol", "Ensembl_nuc", "RefSeq_nuc", "MANE_status", "Ensembl_Gene", "Mean GC %"],
            hide_index=True,
            use_container_width=True,
        )

        # Ensure filter columns are always present with safe defaults (defensive guard)
        for _col, _default in [("Exon Filter", "all"), ("Intron Filter", "all"),
                                ("Include 5' UTR", False), ("Include 3' UTR", False),
                                ("Include Intron", False), ("5' Flank", 0), ("3' Flank", 0)]:
            if _col not in edited_df.columns:
                edited_df[_col] = _default

        # Warn about genes that could not be mapped at all
        still_missing = set(g.upper() for g in genes) - set(edited_df['symbol'].str.upper())
        if still_missing:
            st.warning(f"⚠️ Could not find transcripts for: **{', '.join(sorted(still_missing))}**")

        st.divider()

        # ── Generate button ────────────────────────────────────────────────────
        _, btn_col, _ = st.columns([2, 1, 2])
        with btn_col:
            generate = st.button("⚡ Generate BED File", use_container_width=True, type="primary")

        if generate:
            st.sidebar.markdown(f"**Active build:** `{genome_build}`")

            # Quick connectivity check — warn early if Ensembl is unreachable
            import urllib.request as _ureq
            _api_server = "https://rest.ensembl.org" if genome_build == "hg38" else "https://grch37.rest.ensembl.org"
            try:
                with _ureq.urlopen(f"{_api_server}/info/ping?content-type=application/json", timeout=8) as _r:
                    pass  # just checking we can connect
            except Exception:
                st.warning(
                    "⚠️ Could not reach the Ensembl REST API. Check your internet connection, "
                    "or try again in a few minutes."
                )

            progress_bar = st.progress(0)
            status_ph    = st.empty()

            bed_lines      = []
            total          = len(edited_df)
            no_cds_genes   = []
            fallback_genes = []
            refseq_issues  = []
            exon_filter_errors = []   # {gene, requested, available}

            for i, (index, row) in enumerate(edited_df.iterrows()):
                progress_bar.progress((i + 1) / total)

                gene          = row['symbol']
                enst_req      = row['Ensembl_nuc']
                nm_req        = row.get('RefSeq_nuc')

                use_5utr      = row["Include 5' UTR"]
                use_3utr      = row["Include 3' UTR"]
                use_intron    = row.get('Include Intron', False)
                flank_5prime  = int(row.get("5' Flank", 0) or 0)
                flank_3prime  = int(row.get("3' Flank", 0) or 0)
                _ef_raw = row.get("Exon Filter", "all")
                exon_filter_raw   = "all" if pd.isna(_ef_raw) else str(_ef_raw).strip()
                _if_raw = row.get("Intron Filter", "all")
                intron_filter_raw = "all" if pd.isna(_if_raw) else str(_if_raw).strip()

                # Parse exon filter — report syntax errors immediately and skip gene
                try:
                    exon_filter_set = parse_exon_filter(exon_filter_raw)
                except ValueError as e:
                    st.error(f"❌ **{gene}** — invalid Exon Filter '{exon_filter_raw}': {e}")
                    continue

                # Parse intron filter — same logic
                try:
                    intron_filter_set = parse_exon_filter(intron_filter_raw)  # reuse same parser
                except ValueError as e:
                    st.error(f"❌ **{gene}** — invalid Intron Filter '{intron_filter_raw}': {e}")
                    continue

                status_ph.caption(f"⏳ Fetching {gene}  ({enst_req})…")

                transcript_info = mane_utils.fetch_ensembl_exons(enst_req, build=genome_build, gene_symbol=gene)

                if transcript_info:
                    fetched_id = transcript_info['transcript_id']

                    # Detect transcript switch (hg19 fallback)
                    is_switched = fetched_id.split('.')[0] != enst_req.split('.')[0]
                    if is_switched:
                        fallback_genes.append({
                            'Gene': gene,
                            'Requested (MANE)': enst_req,
                            'Used (hg19 canonical)': fetched_id,
                        })

                    # Determine output ID label
                    output_label = fetched_id
                    if "RefSeq" in id_format:
                        if is_switched:
                            refseq_issues.append({'Gene': gene, 'Reason': 'Transcript switched (hg19)'})
                        elif pd.isna(nm_req) or not str(nm_req).startswith("NM"):
                            refseq_issues.append({'Gene': gene, 'Reason': 'No RefSeq ID available'})
                        else:
                            output_label = nm_req

                    # CDS presence check (is None — avoid false positive at coord 0)
                    if transcript_info.get('cds_start') is None:
                        no_cds_genes.append(gene)

                    # Validate exon filter against actual exon count
                    n_exons = len(transcript_info.get('exons', []))
                    if exon_filter_set:
                        invalid_exons = sorted(e for e in exon_filter_set if e < 1 or e > n_exons)
                        if invalid_exons:
                            exon_filter_errors.append({
                                'Gene': gene,
                                'Requested': exon_filter_raw,
                                'Invalid exon(s)': ', '.join(str(e) for e in invalid_exons),
                                'Transcript exons': n_exons,
                            })
                            exon_filter_set = exon_filter_set - set(invalid_exons)
                            if not exon_filter_set:
                                st.error(f"❌ **{gene}** — all requested exons are out of range "
                                         f"(transcript has {n_exons} exon(s)). Gene skipped.")
                                continue

                    # Validate intron filter — a transcript with N exons has N-1 introns
                    n_introns = max(n_exons - 1, 0)
                    if intron_filter_set:
                        invalid_introns = sorted(n for n in intron_filter_set if n < 1 or n > n_introns)
                        if invalid_introns:
                            exon_filter_errors.append({
                                'Gene': gene,
                                'Requested': intron_filter_raw,
                                'Invalid exon(s)': ', '.join(str(n) for n in invalid_introns),
                                'Transcript exons': f"{n_introns} introns",
                            })
                            intron_filter_set = intron_filter_set - set(invalid_introns)
                            if not intron_filter_set and use_intron:
                                st.error(f"❌ **{gene}** — all requested introns are out of range "
                                         f"(transcript has {n_introns} intron(s)). Introns skipped for this gene.")
                                intron_filter_set = None  # fall back to no introns
                                use_intron = False

                    records = mane_utils.generate_bed_records(
                        transcript_info,
                        include_5utr=use_5utr,
                        include_3utr=use_3utr,
                        include_intron=use_intron,
                        flank_5_prime=flank_5prime,
                        flank_3_prime=flank_3prime,
                        exon_filter=exon_filter_set,
                        intron_filter=intron_filter_set,
                    )

                    for chrom, start, end, type_ in records:
                        p_start = max(0, start - padding)
                        p_end   = end + padding

                        # Column 8: exon/intron number (. for non-numbered regions)
                        m = re.search(r'(?:exon|intron)(\d+)', type_)
                        exon_num = m.group(1) if m else "."

                        full_name = f"{gene}_{output_label}_{type_}"
                        bed_lines.append(
                            f"{chrom}\t{p_start}\t{p_end}\t{full_name}\t.\t{transcript_info['strand']}\t{gene}\t{exon_num}"
                        )
                else:
                    st.error(
                        f"❌ Failed to fetch transcript for **{gene}** ({enst_req}). "
                        "Ensembl API may be temporarily unavailable or the transcript ID is not recognised on this build. "
                        "Try switching genome build or click **Reset & Re-fetch** in the sidebar."
                    )

            progress_bar.empty()
            status_ph.empty()

            # ── Generation Report ──────────────────────────────────────────────
            st.markdown("### 📋 Generation Report")

            rep1, rep2, rep3 = st.columns(3)
            with rep1:
                if fallback_genes:
                    with st.expander(f"⚠️ {len(fallback_genes)} transcript switch(es)", expanded=True):
                        st.caption(f"MANE Select IDs not available on {genome_build}; canonical transcript used instead.")
                        st.table(pd.DataFrame(fallback_genes))
                else:
                    st.success(f"✅ No transcript switches on {genome_build}.")

            with rep2:
                if "RefSeq" in id_format:
                    if refseq_issues:
                        with st.expander(f"⚠️ {len(refseq_issues)} gene(s) reverted to Ensembl ID", expanded=True):
                            st.caption("No valid RefSeq (NM) ID available; Ensembl ID used instead.")
                            st.table(pd.DataFrame(refseq_issues))
                    else:
                        st.success("✅ All regions labelled with RefSeq IDs.")
                else:
                    st.info("Ensembl IDs used (RefSeq mode off).")

            with rep3:
                if no_cds_genes:
                    st.warning(f"⚠️ No CDS data for: {', '.join(no_cds_genes)}\n\nFull exons included without CDS labels.")
                else:
                    st.success("✅ CDS data found for all transcripts.")

            # Exon filter validation report (full-width, below the 3-col block)
            if exon_filter_errors:
                with st.expander(f"⚠️ {len(exon_filter_errors)} gene(s) had out-of-range exon numbers", expanded=True):
                    st.caption(
                        "The exon numbers below do not exist in the selected transcript. "
                        "Valid exons were still included; only the invalid numbers were skipped."
                    )
                    st.table(pd.DataFrame(exon_filter_errors))

            # ── Download block ─────────────────────────────────────────────────
            if not bed_lines:
                st.error("No BED records were generated. All gene lookups may have failed — check errors above.")
            else:
                bed_content = "\n".join(bed_lines)
                n_genes_out = edited_df['symbol'].nunique()

                st.markdown(f"""
                <div class="download-summary">
                    <h3>✅ Panel ready</h3>
                    <p>{len(bed_lines):,} regions across {n_genes_out} gene(s) · build {genome_build} · {padding} bp padding</p>
                </div>
                """, unsafe_allow_html=True)

                st.download_button(
                    label="⬇️  Download BED File",
                    data=bed_content,
                    file_name=f"panel_{genome_build}.bed",
                    mime="text/plain",
                    use_container_width=True,
                    type="primary",
                )

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 2 — LIFTOVER
# ═══════════════════════════════════════════════════════════════════════════════
with tab2:

    st.markdown("""
    <div class="info-card">
        <strong style="color:#00c9a7;">🔄 Coordinate LiftOver — hg19 → hg38</strong><br/>
        Convert an existing BED file from GRCh37 (hg19) to GRCh38 (hg38) coordinates using
        UCSC chain files. Regions that cannot be mapped are collected separately so nothing
        is silently dropped.
    </div>
    """, unsafe_allow_html=True)

    lift_file = st.file_uploader(
        "Upload hg19 BED file",
        type=["bed", "txt"],
        key="lift_upload",
        help="Accepts tab-separated, comma-separated or space-separated BED files. Headers (track/browser/#) are preserved."
    )

    if lift_file:
        content = lift_file.read().decode("utf-8")
        n_lines = sum(1 for l in content.splitlines() if l.strip() and not l.startswith(('#', 'track', 'browser')))
        st.caption(f"📄 {lift_file.name} · {n_lines} data line(s) detected")

        _, btn_col2, _ = st.columns([2, 1, 2])
        with btn_col2:
            convert = st.button("🔄 Convert to hg38", use_container_width=True, type="primary")

        if convert:
            with st.spinner("Converting coordinates…"):
                converted, unmapped = liftover_utils.convert_bed_hg19_to_hg38(content)

            n_converted = sum(1 for l in converted.splitlines() if l.strip() and not l.startswith(('#', 'track', 'browser')))
            n_unmapped  = sum(1 for l in unmapped.splitlines() if l.strip())

            m1, m2 = st.columns(2)
            m1.metric("Regions converted", n_converted)
            m2.metric("Regions unmapped",  n_unmapped, delta_color="inverse" if n_unmapped else "off")

            st.divider()

            dl1, dl2 = st.columns(2)
            with dl1:
                st.download_button(
                    "⬇️  Download hg38 BED",
                    data=converted,
                    file_name="hg38_converted.bed",
                    mime="text/plain",
                    use_container_width=True,
                    type="primary",
                )
            with dl2:
                if unmapped.strip():
                    st.download_button(
                        "⬇️  Download Unmapped Regions",
                        data=unmapped,
                        file_name="unmapped_lines.txt",
                        mime="text/plain",
                        use_container_width=True,
                    )
                else:
                    st.success("✅ All regions mapped successfully — no unmapped lines.")
