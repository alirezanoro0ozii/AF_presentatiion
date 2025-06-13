import requests
import sys

def download_alphafold_pdb(uniprot_id, out_path=None):
    """Download the AlphaFold DB prediction for a UniProt ID."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    # url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}.pdb"
    r = requests.get(url)
    r.raise_for_status()
    fn = out_path or f"data/{uniprot_id}_alphafold.pdb"
    with open(fn, "wb") as f:
        f.write(r.content)
    print(f"✅ AlphaFold PDB written to {fn}")
    return fn

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

def download_experimental_pdb(pdb_id, out_path=None):
    """Download the experimental PDB by its PDB code."""
    url = "https://files.rcsb.org/download/7VCF.cif"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url)
    r.raise_for_status()
    fn = out_path or f"data/{pdb_id}.pdb"
    with open(fn, "wb") as f:
        f.write(r.content)
    print(f"✅ Experimental PDB written to {fn}")
    return fn

if __name__ == "__main__":
    uni_id = sys.argv[1] if len(sys.argv) > 1 else "A0A000"  # Default to SARS-CoV-2 spike protein
    # 1) Download AF2 prediction
    af_pdb = download_alphafold_pdb(uni_id)

    # 2) Find and download an experimental structure
    try:
        pdb_code = fetch_best_pdb_id(uni_id)
    except Exception as e:
        pdb_code = "MoeA5"
    exp_pdb = download_experimental_pdb(pdb_code)


import requests
import os

# 1) Find all PDBs mapped to the UniProt accession
def get_pdb_ids_for_uniprot(uniprot_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{uniprot_id}"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json().get(uniprot_id, {}).get("PDB", {})
    return list(data.keys())

# 2) Filter for truly experimental entries
def is_experimental(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/method/{pdb_id.lower()}"
    resp = requests.get(url)
    resp.raise_for_status()
    methods = resp.json().get(pdb_id.lower(), {}).get("experimental_method", [])
    # keep only if one of these methods appears
    allowed = {"X-ray diffraction", "NMR", "Electron microscopy"}
    return any(m in allowed for m in methods)

# 3) Download the PDB from RCSB
def download_pdb_from_rcsb(pdb_id, out_dir="."):
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(url)
    resp.raise_for_status()
    path = os.path.join(out_dir, f"{pdb_id}.pdb")
    with open(path, "wb") as f:
        f.write(resp.content)
    print(f"Downloaded experimental PDB: {path}")
    return path

if __name__ == "__main__":
    uni = "P69905"  # replace with your UniProt accession
    outdir = "structures"
    os.makedirs(outdir, exist_ok=True)

    pdbs = get_pdb_ids_for_uniprot(uni)
    print(f"Found PDB IDs for {uni}: {pdbs}")

    for pdb in pdbs:
        if is_experimental(pdb):
            download_pdb_from_rcsb(pdb, outdir)
        else:
            print(f"Skipping {pdb}: not experimental (likely AFDB or other model)")
