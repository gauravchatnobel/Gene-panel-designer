import requests
import json

def check_bulk(genes):
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens"
    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    
    # We need expand=1 to get transcripts? 
    # The POST endpoint documentation says: "dict of lists of strings" for symbols.
    # To get transcripts, we might need to add "expand" param to the URL even for POST?
    # Or maybe it returns minimal info.
    # Let's try adding expand=1 to URL.
    
    r = requests.post(server+ext + "?expand=1", headers=headers, data=json.dumps({ "symbols" : genes }))
    
    if not r.ok:
        print(f"Failed: {r.status_code} {r.text}")
        return
        
    data = r.json()
    for gene, info in data.items():
        if info:
            print(f"Gene: {gene}")
            transcripts = info.get('Transcript', [])
            found = False
            for t in transcripts:
                if t.get('is_canonical'):
                    print(f"  Canonical: {t['id']} ({t.get('version')})")
                    found = True
                    break
            if not found:
                print("  No canonical found.")
        else:
            print(f"Gene: {gene} -> Not Found")

check_bulk(["H3C2"])
