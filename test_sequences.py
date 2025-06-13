import os
from utils.utils import read_fasta_sequences

def check_pdb_exists(seq_id, pdb_dir):
    pdb_path = os.path.join(pdb_dir, f"{seq_id}.pdb")
    return os.path.exists(pdb_path)

# Read sequences from casp-14.txt
sequences = read_fasta_sequences('casp-14.txt')

# Check each sequence against PDB files
pdb_dir = 'casp14_targets'
missing_pdbs = []

for seq_id in sequences:
    if not check_pdb_exists(seq_id, pdb_dir):
        missing_pdbs.append(seq_id)

# Print results
print(f"Total sequences found: {len(sequences)}")
print(f"Missing PDB files: {len(missing_pdbs)}")
if missing_pdbs:
    print("\nMissing PDB files for:")
    for seq_id in missing_pdbs:
        print(f"- {seq_id}")
