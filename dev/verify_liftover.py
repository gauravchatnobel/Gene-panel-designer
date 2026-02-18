from liftover_utils import convert_bed_hg19_to_hg38

# Test hg19 coordinates (chr7:140453136 for BRAF V600E) -> hg38 should be chr7:140753336 approx
# Note: LiftOver constructor downloads heavy chain files, so first run might be slow.
bed_content = "chr7\t140453130\t140453140\tBRAF_TEST\n"
print("Converting test coordinate...")

try:
    conv, unmap = convert_bed_hg19_to_hg38(bed_content)
    print("Converted:\n", conv)
    print("Unmapped:\n", unmap)
    
    if "14075333" in conv: # Rough check
        print("LiftOver Success!")
    else:
        print("LiftOver Result Unexpected.")

except Exception as e:
    print("LiftOver Failed:", e)
