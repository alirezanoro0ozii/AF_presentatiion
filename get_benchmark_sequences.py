import os
import json
from utils.utils import get_sequence_from_pdb


seq2pdb = {}

for uid in os.listdir('CASP/casp14_targets'):
    pdb_name = f"CASP/casp14_targets/{uid}"
    sequence = get_sequence_from_pdb(pdb_name)
    seq2pdb[sequence] = pdb_name

for uid in os.listdir('CASP/casp14_targets_domain'):
    pdb_name = f"CASP/casp14_targets_domain/{uid}"
    sequence = get_sequence_from_pdb(pdb_name)
    seq2pdb[sequence] = pdb_name

with open('seq2pdb.json', 'w') as f:
    json.dump(seq2pdb, f)
    

Af2_seq2pdb = {}

for uid in os.listdir('pdbs/AF2_paper'):
    pdb_name = f"pdbs/AF2_paper/{uid}"
    sequence = get_sequence_from_pdb(pdb_name)
    Af2_seq2pdb[sequence] = pdb_name

with open('Af2_seq2pdb.json', 'w') as f:
    json.dump(Af2_seq2pdb, f)


# Read the JSON files
with open('seq2pdb.json', 'r') as f:
    seq2pdb = json.load(f)

with open('Af2_seq2pdb.json', 'r') as f:
    Af2_seq2pdb = json.load(f)

# Find sequences that exist in both dictionaries
common_sequences = set(seq2pdb.keys()) & set(Af2_seq2pdb.keys())

print(f"Found {len(common_sequences)} sequences that exist in both datasets")

# Create a dictionary with common sequences and their corresponding PDB files
common_seq2pdb = {}
for sequence in common_sequences:
    common_seq2pdb[sequence] = {
        'casp14_pdb': seq2pdb[sequence],
        'af2_pdb': Af2_seq2pdb[sequence]
    }

# Save the common sequences to a new JSON file
with open('common_sequences.json', 'w') as f:
    json.dump(common_seq2pdb, f, indent=2)

print(f"Common sequences saved to 'common_sequences.json'")
print(f"Sample of common sequences:")
for i, (seq, pdbs) in enumerate(list(common_seq2pdb.items())[:5]):
    print(f"Sequence {i+1}: {seq[:50]}...")
    print(f"  CASP14: {pdbs['casp14_pdb']}")
    print(f"  AF2: {pdbs['af2_pdb']}")
    print()


with open('common_sequences.json', 'r') as f:
    common_seq2pdb = json.load(f)

print(len(common_seq2pdb))