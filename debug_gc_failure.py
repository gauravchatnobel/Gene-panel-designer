import requests
import json
import requests
import json

def text_gc_manual():
    # ID from screenshot: ACTB = ENST00000646664.1
    # 1. Test Versioned vs Unversioned on hg38
    id_ver = "ENST00000646664.1"
    id_base = "ENST00000646664"
    
    server_38 = "https://rest.ensembl.org"
    server_19 = "https://grch37.rest.ensembl.org"
    ext = "/sequence/id"
    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}

    print("--- TEST 1: hg38 (Correct Server) ---")
    payload = { "ids": [id_ver], "type": "cdna" }
    r = requests.post(server_38+ext, headers=headers, data=json.dumps(payload))
    print(f"Versioned ({id_ver}): {r.status_code}")
    if r.ok: print("   -> Success") 
    else: print(f"   -> Fail: {r.text}")

    payload = { "ids": [id_base], "type": "cdna" }
    r = requests.post(server_38+ext, headers=headers, data=json.dumps(payload))
    print(f"Base ({id_base}): {r.status_code}")
    if r.ok: print("   -> Success")

    print("\n--- TEST 2: hg19 (Incorrect Server for MANE IDs) ---")
    payload = { "ids": [id_ver], "type": "cdna" }
    r = requests.post(server_19+ext, headers=headers, data=json.dumps(payload))
    print(f"Versioned ({id_ver}): {r.status_code}")
    
    payload = { "ids": [id_base], "type": "cdna" }
    r = requests.post(server_19+ext, headers=headers, data=json.dumps(payload))
    print(f"Base ({id_base}): {r.status_code}")

text_gc_manual()
