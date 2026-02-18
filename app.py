import streamlit as st
import pandas as pd
import mane_utils
import liftover_utils
import io
import re

st.set_page_config(page_title="Gene Panel Designer", layout="wide")

st.title("Gene Panel Designer")

# Load MANE Data
@st.cache_data
def load_data():
    return mane_utils.load_mane_data()

mane_df = load_data()

if mane_df is None:
    st.error("MANE Summary file not found! Please ensure 'MANE.summary.txt' is in the application directory.")
    st.stop()

tab1, tab2 = st.tabs(["Panel Designer", "LiftOver Tool (hg19 -> hg38)"])

with tab1:
    st.markdown("""
    ### Identify clinically relevant MANE Select transcripts.
    1. Upload gene list.
    2. Review mapping and toggle gene-specific settings (**UTRs**, **Introns**).
    3. Generate BED file.
    """)

    # Initial Sidebar Config (Defaults)
    st.sidebar.header("Panel Settings (Defaults)")
    genome_build = st.sidebar.radio("Genome Build", ["hg19", "hg38"], index=1)
    
    # ID Format Selection
    id_format = st.sidebar.radio(
        "Output Transcript ID Format",
        ["Ensembl (ENST)", "RefSeq (NM)"],
        help="Choose the identifier used in the BED file name column."
    )
    
    st.sidebar.subheader("Transcript Content")
    default_5utr = st.sidebar.checkbox("Include 5' UTR (Default)", value=False, help="Can be customized per gene.")
    default_3utr = st.sidebar.checkbox("Include 3' UTR (Default)", value=False, help="Can be customized per gene.")
    default_intron = st.sidebar.checkbox("Include Introns (Default)", value=False, help="Can be customized per gene.")
    padding = st.sidebar.number_input("Padding (bp)", min_value=0, value=20, help="Recommended: 20bp for clinical panels.")
    
    # Removed sidebar flank inputs here in favor of gene-specific table settings

    uploaded_file = st.file_uploader("Upload Gene List", type=["csv", "xlsx", "txt"], key="designer_upload")

    if uploaded_file:
        # Parse file function...
        genes = []
        try:
            if uploaded_file.name.endswith('.csv'):
                df = pd.read_csv(uploaded_file, header=None)
                genes = df[0].astype(str).tolist()
            elif uploaded_file.name.endswith('.xlsx'):
                df = pd.read_excel(uploaded_file, header=None)
                genes = df[0].astype(str).tolist()
            else: # txt
                content = uploaded_file.read().decode()
                genes = [line.strip() for line in content.splitlines() if line.strip()]
            
            # Clean genes: remove parentheses and extra whitespace
            # e.g. "CD274 (PD-L1)" -> "CD274"
            cleaned_genes = []
            for g in genes:
                # Remove content in parens
                g_clean = re.sub(r'\(.*?\)', '', g)
                # Strip whitespace
                g_clean = g_clean.strip()
                if g_clean:
                    cleaned_genes.append(g_clean)
            genes = cleaned_genes
        except Exception as e:
            st.error(f"Error parsing file: {e}")
            st.stop()
            
        st.success(f"Loaded {len(genes)} genes. Processing...")
        
        # 1. Map to MANE
        # Use a content hash (not just filename) so re-uploading a same-named but different file
        # correctly triggers reprocessing.
        file_hash = hash(uploaded_file.getvalue())
        if 'mapped_data' not in st.session_state or st.session_state.get('last_file_hash') != file_hash:
             mapped_df = mane_utils.get_mane_transcripts(genes, mane_df)
             
             # Fallback for missing genes
             found_symbols = set(mapped_df['symbol'].str.upper())
             input_symbols = set(g.upper() for g in genes)
             missing_genes = input_symbols - found_symbols
             
             if missing_genes:
                 status_text = st.warning(f"Fetching Ensembl Canonical transcripts for {len(missing_genes)} non-MANE genes...")
                 fallback_df = mane_utils.get_ensembl_canonical(missing_genes, build=genome_build)
                 if not fallback_df.empty:
                     mapped_df = pd.concat([mapped_df, fallback_df], ignore_index=True)
                     status_text.empty() # Clear warning
                     
                     # Check for Aliases in Fallback Data
                     if 'original_input' in fallback_df.columns:
                         aliases = fallback_df[fallback_df['original_input'].str.upper() != fallback_df['symbol'].str.upper()]
                         if not aliases.empty:
                             with st.expander(f"ℹ️ {len(aliases)} Gene Aliases Recognized", expanded=True):
                                 st.info(f"The following genes were mapped to their official symbols for {genome_build}:")
                                 st.table(aliases[['original_input', 'symbol']].rename(columns={'original_input': 'Input Symbol', 'symbol': 'Official Symbol'}))
             
             # Fetch GC Content (Batch)
             status_text = st.warning(f"Calculating GC Content for {len(mapped_df)} transcripts...")
             transcript_ids = mapped_df['Ensembl_nuc'].tolist()
             # GC content always fetched from hg38 because MANE Select transcript IDs (ENST) are
             # defined on GRCh38. Ensembl Canonical fallback IDs may silently return None GC for
             # hg19-only transcripts — this is acceptable since GC is for informational display only.
             gc_map = mane_utils.fetch_transcript_gc_batch(transcript_ids, build='hg38')
             
             mapped_df['Mean GC %'] = mapped_df['Ensembl_nuc'].map(gc_map)
             status_text.empty()

             # Add Config Columns
             mapped_df['Include 5\' UTR'] = default_5utr
             mapped_df['Include 3\' UTR'] = default_3utr
             mapped_df['Include Intron'] = default_intron
             
             # NEW: Gene-specific flanks
             # Default to 0, but user can edit per gene
             mapped_df['5\' Flank'] = 0
             mapped_df['3\' Flank'] = 0
             
             st.session_state['mapped_data'] = mapped_df
             st.session_state['last_file_hash'] = file_hash
        
        # Report on Transcript Selection (User Request)
        non_mane_select = st.session_state['mapped_data'][
            st.session_state['mapped_data']['MANE_status'] != 'MANE Select'
        ]
        
        if not non_mane_select.empty:
            with st.expander(f"⚠️ {len(non_mane_select)} Genes used non-MANE Select transcripts", expanded=True):
                st.write("The following genes did not have a MANE Select transcript found. Alternate transcripts were chosen:")
                # Format for display
                display_warn = non_mane_select[['symbol', 'Ensembl_nuc', 'MANE_status']].rename(
                    columns={'symbol': 'Gene', 'Ensembl_nuc': 'Transcript ID', 'MANE_status': 'Reason / Source'}
                )
                st.table(display_warn)
        else:
            if genome_build == 'hg19':
                st.info("✅ All genes mapped to MANE Select. (Note: Runtime switches to hg19 Canonical may occur during generation if MANE versions are incompatible).")
            else:
                st.success("✅ All mapped genes are using MANE Select transcripts.")
        
        # Interactive Editor
        st.subheader("Customize Gene Settings")
        st.info("Tip: Set specific Promoter (5' Flank) regions here. (e.g. 2000 for MYC/BCL2, 500 for TERT/IG).")
        
        edited_df = st.data_editor(
            st.session_state['mapped_data'],
            column_config={
                "Include 5' UTR": st.column_config.CheckboxColumn("5' UTR", help="Toggle 5' UTR inclusion"),
                "Include 3' UTR": st.column_config.CheckboxColumn("3' UTR", help="Toggle 3' UTR inclusion"),
                "Include Intron": st.column_config.CheckboxColumn("Introns", help="Toggle Intron inclusion"),
                "5' Flank": st.column_config.NumberColumn("5' Flank (bp)", help="Upstream Promoter region", min_value=0),
                "3' Flank": st.column_config.NumberColumn("3' Flank (bp)", help="Downstream region", min_value=0),
                "Mean GC %": st.column_config.NumberColumn(
                    "GC %",
                    help="Mean GC content of the transcript (cDNA)",
                    format="%.1f%%",
                    min_value=0,
                    max_value=100
                ),
                "symbol": st.column_config.TextColumn("Gene", disabled=True),
                "Ensembl_nuc": st.column_config.TextColumn("Transcript", disabled=True),
                "MANE_status": st.column_config.TextColumn("Type", disabled=True),
                "GRCh38_chr": None, # Hide this column
            },
            disabled=["symbol", "Ensembl_nuc", "RefSeq_nuc", "MANE_status", "Ensembl_Gene", "Mean GC %"],
            hide_index=True,
            use_container_width=True
        )
        
        missing = set(g.upper() for g in genes) - set(edited_df['symbol'].str.upper())
        if missing:
            st.warning(f"Genes not found in MANE: {', '.join(missing)}")
            
        # Process Code
        if st.button("Generate Final BED"):
            st.sidebar.markdown(f"**Selected:** {genome_build}")
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            bed_lines = []
            total = len(edited_df)
            no_cds_genes = []
            fallback_genes = []
            refseq_issues = [] # Track genes that reverted to ENST in RefSeq mode
            
            for i, (index, row) in enumerate(edited_df.iterrows()):
                progress_value = (i + 1) / total
                count_progress = min(progress_value, 1.0)
                progress_bar.progress(count_progress)
                
                gene = row['symbol']
                enst_requested = row['Ensembl_nuc']
                nm_requested = row.get('RefSeq_nuc') # Might be NaN
                
                use_5utr = row["Include 5' UTR"]
                use_3utr = row["Include 3' UTR"]
                use_intron = row.get('Include Intron', False)
                # Read flanks from row
                flank_5_prime = row.get("5' Flank", 0)
                flank_3_prime = row.get("3' Flank", 0)
                
                status_text.text(f"Fetching {gene} ({enst_requested})...")
                
                # Pass gene_symbol to enable fallback lookup on specific build if ID fails
                transcript_info = mane_utils.fetch_ensembl_exons(enst_requested, build=genome_build, gene_symbol=gene)
                
                if transcript_info:
                    fetched_id = transcript_info['transcript_id']
                    
                    # Check for Switch
                    is_switched = False
                    if fetched_id.split('.')[0] != enst_requested.split('.')[0]:
                        fallback_genes.append({
                            'Gene': gene,
                            'Original (MANE)': enst_requested,
                            'Used (hg19)': fetched_id
                        })
                        is_switched = True
                    
                    # Determine Output ID Label
                    output_label = fetched_id # Default to Ensembl
                    
                    if "RefSeq" in id_format:
                        # Logic: Use RefSeq ONLY if:
                        # 1. No runtime switch (because original NM matches original ENST)
                        # 2. We actually have an NM ID
                        
                        if is_switched:
                            refseq_issues.append({'Gene': gene, 'Reason': 'Transcript Switched (hg19)'})
                        elif pd.isna(nm_requested) or not str(nm_requested).startswith("NM"):
                             refseq_issues.append({'Gene': gene, 'Reason': 'No RefSeq ID available'})
                        else:
                             output_label = nm_requested

                    # Check for CDS presence (use `is None` to avoid false positive when cds_start == 0)
                    if transcript_info.get('cds_start') is None:
                        no_cds_genes.append(gene)
                        
                    records = mane_utils.generate_bed_records(
                        transcript_info, 
                        include_5utr=use_5utr, 
                        include_3utr=use_3utr, 
                        include_intron=use_intron,
                        flank_5_prime=flank_5_prime,
                        flank_3_prime=flank_3_prime
                    )
                    
                    for chrom, start, end, type_ in records:
                        p_start = max(0, start - padding)
                        p_end = end + padding
                        
                        # Extract Exon Number for extra column
                        # type_ is usually 'exon7', 'exon7_CDS', 'intron1', '5UTR' (if not attached to exon?) 
                        # actually mane_utils usually outputs 'exonX_...'
                        exon_match = re.search(r'(?:exon|intron)(\d+)', type_)
                        exon_num = exon_match.group(1) if exon_match else "."
                        
                        full_name = f"{gene}_{output_label}_{type_}"
                        
                        # BED6 + 2 (Gene, Exon)
                        # Standard BED has 6 columns. We append Gene and Exon as 7 and 8.
                        bed_lines.append(f"{chrom}\t{p_start}\t{p_end}\t{full_name}\t.\t{transcript_info['strand']}\t{gene}\t{exon_num}")
                else:
                    st.error(f"Failed to fetch: {gene}")
               
            # ---------------------------
            # Final Generation Report
            # ---------------------------
            st.divider()
            st.subheader("Generation Report")
            
            # 1. Transcript Switches (hg19)
            if fallback_genes:
                with st.expander(f"⚠️ {len(fallback_genes)} Transcripts Changed (hg19 Fallback)", expanded=True):
                    st.warning(f"The original MANE Select transcript could not be found on {genome_build}. The system fell back to the canonical transcript for this build:")
                    st.table(pd.DataFrame(fallback_genes))
            else:
                 st.success(f"✅ No transcript switches occurred. All genes used the intended MANE Select transcript ID on {genome_build}.")

            # 2. RefSeq ID Status (if selected)
            if "RefSeq" in id_format:
                if refseq_issues:
                     with st.expander(f"⚠️ {len(refseq_issues)} Genes Reverted to Ensembl ID", expanded=True):
                        st.warning("The following genes did not have a valid RefSeq match suitable for output (or were switched at runtime), so we used the Ensembl ID (ENST) instead:")
                        st.table(pd.DataFrame(refseq_issues))
                else:
                     st.success("✅ All regions labeled with RefSeq (NM) IDs.")

            # 3. CDS Status
            if no_cds_genes:
                st.warning(f"⚠️ No Coding Sequence (CDS) data found for the following genes on {genome_build}: {', '.join(no_cds_genes)}.\n\nFull exons were included without '_CDS' labels.")
            else:
                st.success(f"✅ Coding Sequence (CDS) data found for all transcripts.")
                    
            if not bed_lines:
                status_text.error("No BED records were generated. All gene lookups may have failed.")
            else:
                status_text.success(f"Generation Complete! {len(bed_lines)} records generated.")
                st.download_button(
                    "Download BED",
                    "\n".join(bed_lines),
                    file_name=f"panel_{genome_build}.bed",
                    mime="text/plain"
                )

with tab2:
    st.header("LiftOver Assistant (hg19 -> hg38)")
    st.markdown("Convert an existing BED file from hg19 (GRCh37) to hg38 (GRCh38) coordinates.")
    
    lift_file = st.file_uploader("Upload hg19 BED File", type=["bed", "txt"], key="lift_upload")
    
    if lift_file:
        content = lift_file.read().decode("utf-8")
        if st.button("Convert to hg38"):
            with st.spinner("Converting..."):
                converted, unmapped = liftover_utils.convert_bed_hg19_to_hg38(content)
            
            st.success("Conversion finished.")
            
            col1, col2 = st.columns(2)
            with col1:
                st.download_button("Download hg38 BED", converted, "hg38_converted.bed")
            with col2:
                if unmapped.strip():
                     st.download_button("Download Unmapped Lines", unmapped, "unmapped_lines.txt")
                else:
                     st.info("All regions mapped successfully.")
