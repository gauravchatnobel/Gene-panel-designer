import mane_utils
import pandas as pd
import unittest
from unittest.mock import MagicMock

# Mock MANE Data
mock_df = pd.DataFrame({
    'symbol': ['BRCA1', 'BRCA1', 'TP53', 'GENEX', 'GENEY'],
    'Ensembl_nuc': ['ENST001', 'ENST002', 'ENST003', 'ENST004', 'ENST005'],
    'RefSeq_nuc': ['NM_001', 'NM_002', 'NM_003', 'NM_004', 'NM_005'],
    'MANE_status': ['MANE Select', 'MANE Plus Clinical', 'MANE Select', 'MANE Plus Clinical', 'MANE Plus Clinical'],
    'Ensembl_Gene': ['ENSG1', 'ENSG1', 'ENSG2', 'ENSG4', 'ENSG5'],
    'GRCh38_chr': ['17', '17', '17', '1', '2']
})

class TestManeLogic(unittest.TestCase):
    def test_prioritization(self):
        # BRCA1 has both. Should pick Select.
        # TP53 has Select.
        # GENEX has Plus only.
        # GENEY has Plus only.
        
        genes = ['BRCA1', 'TP53', 'GENEX', 'MISSING']
        
        result = mane_utils.get_mane_transcripts(genes, mock_df)
        
        # Check BRCA1
        brca1 = result[result['symbol'] == 'BRCA1']
        self.assertEqual(len(brca1), 1)
        self.assertEqual(brca1.iloc[0]['MANE_status'], 'MANE Select')
        
        # Check TP53
        tp53 = result[result['symbol'] == 'TP53']
        self.assertEqual(len(tp53), 1)
        self.assertEqual(tp53.iloc[0]['MANE_status'], 'MANE Select')
        
        # Check GENEX (only has Plus)
        genex = result[result['symbol'] == 'GENEX']
        self.assertEqual(len(genex), 1)
        self.assertEqual(genex.iloc[0]['MANE_status'], 'MANE Plus Clinical')
        
        print("Logic Verification Passed!")

if __name__ == '__main__':
    unittest.main()
