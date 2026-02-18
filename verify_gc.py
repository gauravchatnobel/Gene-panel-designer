import mane_utils

def text_gc_batch():
    # Test IDs from different builds to be safe, but we primarily target hg38 for now
    ids = [
        "ENST00000344686", # MAP3K14
        "ENST00000641515", # GAPDH (hg38 MANE)
        "ENST00000621411"  # H3C2
    ]
    
    print("Fetching GC for batch (hg38)...")
    results = mane_utils.fetch_transcript_gc_batch(ids, build='hg38')
    
    for tid, gc in results.items():
        if gc is not None:
            print(f"{tid}: {gc:.2f}%")
        else:
            print(f"{tid}: Failed")

text_gc_batch()
