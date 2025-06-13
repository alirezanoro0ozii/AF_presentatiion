from utils.utils import tm_align, esm_fold

query = "T1025_ESM.pdb"
reference = "casp14_targets/T1025.pdb"

pdb_file, plddt = esm_fold(query, 'T1025')
print(tm_align(pdb_file, reference, 'T1025'))