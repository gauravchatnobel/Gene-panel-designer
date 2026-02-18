import requests
import json
import time

def calculate_gc(seq):
    if not seq: return 0
    g = seq.count('G') + seq.count('g')
    c = seq.count('C') + seq.count('c')
    return (g + c) / len(seq) * 100

def fetch_batch_gc(transcript_ids):
    server = "https://rest.ensembl.org"
    ext = "/sequence/id"
    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    
    # Payload for POST
    # Note: Ensembl allows up to 50 IDs per call usually
    payload = {
        "ids" : transcript_ids,
        "type": "cdna" 
    }
    
    print(f"Fetching sequences for {len(transcript_ids)} transcripts...")
    t0 = time.time()
    r = requests.post(server+ext, headers=headers, data=json.dumps(payload))
    t1 = time.time()
    
    if not r.ok:
        print(f"Failed: {r.status_code}")
        print(r.text)
        return

    data = r.json()
    print(f"Success! Took {t1-t0:.2f}s")
    
    for item in data[:3]: # Show first 3 result stats
        tid = item['id']
        seq = item['seq']
        gc = calculate_gc(seq)
        print(f"{tid}: Len={len(seq)}, GC={gc:.2f}%")

# Test with a few mixed IDs
test_ids = ["ENST00000641515", "ENST00000361390", "ENST00000361453", "ENST00000361624", "ENST00000361739"]
fetch_batch_gc(test_ids)
