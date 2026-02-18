import urllib.request
import json
import ssl

def check_gene(gene_symbol):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1;content-type=application/json"
    
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    req = urllib.request.Request(server+ext, headers={ "Content-Type" : "application/json"})
    
    try:
        with urllib.request.urlopen(req, context=ctx) as r:
            data = json.loads(r.read().decode())
            
        print(f"Gene: {data.get('display_name')}")
        
        transcripts = data.get('Transcript', [])
        print(f"Found {len(transcripts)} transcripts.")
        
        for t in transcripts:
            display_name = t.get('display_name')
            # Inspect structure
            if t['id'] == transcripts[0]['id']:
                 print("Sample Transcript Keys:", t.keys())
            
            # Look for MANE
            if 'MANE' in str(t):
                print(f"Found MANE potential in {t['id']} ({display_name})")
                
    except Exception as e:
        print(f"Failed: {e}")

check_gene("BRCA1")
