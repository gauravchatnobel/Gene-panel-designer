from pyliftover import LiftOver
import os
import pandas as pd
import io

def convert_bed_hg19_to_hg38(bed_content):
    """
    Convert hg19 BED content (str) to hg38.
    Returns: (converted_bed_str, unmapped_bed_str)
    """
    # Initialize LiftOver (downloads chain file to /tmp or cache on first run)
    # Using 'hg19' to 'hg38'
    lo = LiftOver('hg19', 'hg38')
    
    converted_lines = []
    unmapped_lines = []
    
    # Process line by line
    for line in bed_content.splitlines():
        # Strip BOM if present (often at start of first line)
        line = line.lstrip('\ufeff')
        
        if not line.strip() or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
            converted_lines.append(line) # Keep headers
            continue
            
        # Detect delimiter: try tab first, then comma
        # Detect delimiter: try tab first, then comma, then whitespace
        parts = line.strip().split('\t')
        if len(parts) < 3:
            parts = line.strip().split(',')
            parts = [p.strip('"') for p in parts] # Strip quotes for CSV
            
        if len(parts) < 3:
             # Fallback to general whitespace (handles spaces)
             parts = line.strip().split()
            
        if len(parts) < 3:
            unmapped_lines.append(line + "\t# Invalid Format")
            continue
            
        chrom = parts[0]
        # Normalize chromosome (Ensembl '1' -> UCSC 'chr1')
        if not chrom.startswith('chr'):
            chrom = f"chr{chrom}"
            
        # Handle MT -> chrM special case if needed, usually M -> chrM
        if chrom == 'chrM': 
            pass # Standard
        elif chrom == 'chrMT':
             chrom = 'chrM'
             
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError:
             unmapped_lines.append(line + "\t# Invalid Coordinates")
             continue
             
        rest = parts[3:] if len(parts) > 3 else []
        
        # LiftOver start and end
        # LiftOver returns list of matches. Usually one.
        # Check start
        new_start_list = lo.convert_coordinate(chrom, start)
        new_end_list = lo.convert_coordinate(chrom, end)
        
        if not new_start_list or not new_end_list:
            unmapped_lines.append(line + "\t# Mapping Failed")
            continue
            
        # Take first match
        new_chrom = new_start_list[0][0]
        new_start = new_start_list[0][1]
        new_end = new_end_list[0][1]
        
        # Check integrity (same chromosome, minimal inversion check though liftover handles this)
        if new_chrom != new_end_list[0][0]:
             unmapped_lines.append(line + "\t# Split Mapping")
             continue
             
        # Reconstruct line
        # Standardize chrom name if needed (e.g. if source didn't have chr)
        # pyliftover target usually follows UCSC convention (chrN)
        new_line = f"{new_chrom}\t{new_start}\t{new_end}"
        if rest:
            new_line += "\t" + "\t".join(rest)
            
        converted_lines.append(new_line)
        
    return "\n".join(converted_lines), "\n".join(unmapped_lines)
