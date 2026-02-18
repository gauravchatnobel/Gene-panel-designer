import requests
import json

def get_canonical(gene_symbol):
    server = "https://rest.ensembl.org"
    # Search for the gene symbol to get the Ensembl Gene ID and Canonical Transcript
    ext = f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1;content-type=application/json"
    
    try:
        r = requests.get(server+ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            print(f"{gene_symbol}: Failed {r.status_code}")
            return
            
        data = r.json()
        
        # In the gene object, 'Transcript' is a list.
        # We look for 'is_canonical': 1
        
        transcripts = data.get('Transcript', [])
        canonical = None
        for t in transcripts:
            if t.get('is_canonical') == 1:
                canonical = t
                break
        
        if canonical:
            print(f"Gene: {gene_symbol}")
            print(f"Canonical ID: {canonical['id']}")
            print(f"Version: {canonical.get('version')}")
            print(f"Canonical Status: Found")
        else:
            print(f"{gene_symbol}: No canonical transcript found in {len(transcripts)} transcripts.")

    except Exception as e:
        print(f"Error: {e}")

# Test with confirmed non-MANE gene from user list if possible, or just a known one.
# User mentioned "HIST1H3B" etc. Let's try one of those.
get_canonical("HIST1H3B")
get_canonical("H2AC11")
