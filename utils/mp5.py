import requests
import json
import matplotlib.pyplot as plt
import os
os.makedirs("pdbs/mp5", exist_ok=True)

def save_pd_file(name: str, distance_matrix: list):
    os.makedirs(f'PD/{name}', exist_ok=True)
    # Plot the distance matrix
    plt.figure(figsize=(10, 8))
    plt.imshow(distance_matrix, cmap='viridis')
    plt.colorbar(label='Distance (Å)')
    plt.title('Pairwise Distance Matrix')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.savefig(f'PD/{name}/MP5_distance_matrix.png')
    plt.close()
        
def mp5_fold(name: str, sequence: str) -> dict:
    """
    Call the 310.ai API to predict PDB structure from protein sequence
    
    Args:
        sequence: Protein amino acid sequence
        
    Returns:
        API response as dictionary
    """
    url = "https://310.ai/api/mp5/v1/predict_pdb"
    
    headers = {
        "Content-Type": "application/json"
    }
    
    data = {
        "sequence": sequence,
        "user_id": "1",
        "method": "hybrid"
    }
    
    # if os.path.exists(f"pdbs/mp5/{name}.pdb"):
    #     return f"pdbs/mp5/{name}.pdb"
    
    os.makedirs(f"MP5_jsons/{name}", exist_ok=True)
    try:
        response = requests.post(url, headers=headers, json=data)
        response.raise_for_status()
        
        save_pd_file(name, response.json()["pd_matrix"])
        
        with open(f"MP5_jsons/{name}/predicted_pdb.json", "w") as f:
            json.dump(response.json(), f, indent=2)
            
        # Download and save PDB file
        pdb_url = response.json()["pdb_file_path"]
        pdb_response = requests.get(pdb_url)
        pdb_response.raise_for_status()
        
        with open(f"pdbs/mp5/{name}.pdb", "w") as f:
            f.write(pdb_response.text)
        return f"pdbs/mp5/{name}.pdb"
    
    except requests.exceptions.RequestException as e:
        raise Exception(f"API request failed: {e}")
    
if __name__ == "__main__":
    from utils import read_fasta_sequences
    sequences = read_fasta_sequences("CASP/casp-14.txt")
    for name, sequence in sequences.items():
        if os.path.exists(f"api_responses/{name}"):
            print(f"Skipping {name} - folder already exists")
            continue
        result = mp5_fold(name, sequence)
        print(result)
        print(result.keys())
        # break