import streamlit as st
import pandas as pd
import mane_utils
import liftover_utils
import io

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
    
    default_utr = st.sidebar.checkbox("Include UTRs (Default)", value=False, help="Can be customized per gene.")
    default_intron = st.sidebar.checkbox("Include Introns (Default)", value=False, help="Can be customized per gene.")
    padding = st.sidebar.number_input("Padding (bp)", min_value=0, value=20, help="Recommended: 20bp for clinical panels.")

    uploaded_file = st.file_uploader("Upload Gene List", type=["csv", "xlsx", "txt"], key="designer_upload")

    if uploaded_file:
        # Parse file
        genes = []
        try:
import re

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
        if 'mapped_data' not in st.session_state or st.session_state.get('last_uploaded') != uploaded_file.name:
             mapped_df = mane_utils.get_mane_transcripts(genes, mane_df)
             # Add Config Columns
             mapped_df['Include UTR'] = default_utr
             mapped_df['Include Intron'] = default_intron
             
             st.session_state['mapped_data'] = mapped_df
             st.session_state['last_uploaded'] = uploaded_file.name
        
        # Interactive Editor
        st.subheader("Customize Gene Settings")
        edited_df = st.data_editor(
            st.session_state['mapped_data'],
            column_config={
                "Include UTR": st.column_config.CheckboxColumn("Include UTR", help="Toggle 5'/3' UTR inclusion"),
                "Include Intron": st.column_config.CheckboxColumn("Include Intron", help="Toggle Intron inclusion"),
                "symbol": st.column_config.TextColumn("Gene", disabled=True),
                "Ensembl_nuc": st.column_config.TextColumn("Transcript", disabled=True),
                "MANE_status": st.column_config.TextColumn("Type", disabled=True),
                "GRCh38_chr": None, # Hide this column
            },
            disabled=["symbol", "Ensembl_nuc", "RefSeq_nuc", "MANE_status", "Ensembl_Gene"],
            hide_index=True,
            use_container_width=True
        )
        
        missing = set(g.upper() for g in genes) - set(edited_df['symbol'].str.upper())
        if missing:
            st.warning(f"Genes not found in MANE: {', '.join(missing)}")
            
        # Process Code
        if st.button("Generate Final BED"):
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            bed_lines = []
            total = len(edited_df)
            
            for i, row in edited_df.iterrows():
                progress_bar.progress((i + 1) / total)
                
                gene = row['symbol']
                enst = row['Ensembl_nuc']
                use_utr = row['Include UTR']
                use_intron = row['Include Intron']
                
                status_text.text(f"Fetching {gene} ({enst})...")
                
                transcript_info = mane_utils.fetch_ensembl_exons(enst, build=genome_build)
                
                if transcript_info:
                    records = mane_utils.generate_bed_records(transcript_info, include_utr=use_utr, include_intron=use_intron)
                    
                    for chrom, start, end, type_ in records:
                        p_start = max(0, start - padding)
                        p_end = end + padding
                        bed_lines.append(f"{chrom}\t{p_start}\t{p_end}\t{gene}_{enst}_{type_}\t.\t{transcript_info['strand']}")
                else:
                    st.error(f"Failed to fetch: {gene}")
                    
            status_text.success("Generation Complete!")
            
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
