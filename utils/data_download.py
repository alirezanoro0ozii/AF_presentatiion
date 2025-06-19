import requests
import os

os.makedirs("pdbs/alphafold", exist_ok=True)
os.makedirs("pdbs/experimental", exist_ok=True)

def download_alphafold_pdb(uniprot_id):
    """Download the AlphaFold DB prediction for a UniProt ID."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    r = requests.get(url)
    r.raise_for_status()
    with open(f"pdbs/alphafold/{uniprot_id}_alphafold.pdb", "wb") as f:
        f.write(r.content)
    print(f"✅ AlphaFold PDB written to pdbs/alphafold/{uniprot_id}_alphafold.pdb")
    return f"pdbs/alphafold/{uniprot_id}_alphafold.pdb"

def fetch_best_pdb_id(uniprot_id):
    """Use PDBe's mapping API to find PDB entries for this UniProt."""
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{uniprot_id}"
    r = requests.get(url)
    r.raise_for_status()
    data = r.json()
    mappings = data.get(uniprot_id, {}).get("PDB", {})
    if not mappings:
        raise ValueError(f"No PDB entries found for UniProt {uniprot_id}")
    # pick the first entry (usually highest coverage)
    pdb_id = next(iter(mappings))
    print(f"🔍 Found PDB ID {pdb_id} for UniProt {uniprot_id}")
    return pdb_id

def download_experimental_pdb(pdb_id):
    """Download the experimental PDB by its PDB code."""
    url = "https://files.rcsb.org/download/7VCF.cif"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url)
    r.raise_for_status()
    with open(f"pdbs/experimental/{pdb_id}.pdb", "wb") as f:
        f.write(r.content)
    print(f"✅ Experimental PDB written to pdbs/experimental/{pdb_id}.pdb")
    return f"pdbs/experimental/{pdb_id}.pdb"
