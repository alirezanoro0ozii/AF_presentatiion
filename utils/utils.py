import requests
import json
import re
import os
import subprocess
from Bio.PDB import PDBParser, Superimposer

def esm_fold(sequence, num_cycle = 6, name = 'esm'):
    os.makedirs(f'esm_fold/{name}', exist_ok=True)
    # url = "http://54.214.144.188/omega/get_pdb_url_from_seq_plddt/"
    url = "http://35.192.181.98:8000/omega/local_nimfold/"
    data = {
        "sequence": sequence.replace(" ", "").replace("\"", ""),
        "num_cycle": num_cycle
    }
    response = requests.post(url, json=data)
    if response.status_code == 200:
        res = response.content.decode("utf-8")
        res_json = json.loads(res)
        plddt = res_json["result"]["pLDDT_mean"]
        pdb_link = res_json["pdb"]
        pdb_file = requests.get(pdb_link)
        with open(f'esm_fold/{name}.pdb', 'w') as f:
            f.write(pdb_file.text)
        return pdb_file.text, plddt
    else:
        return None, None
    
def compute_rmsd(query, reference):
    # Parse structures
    parser = PDBParser(QUIET=True)
    predicted = parser.get_structure('pred', query)
    experimental = parser.get_structure('exp', reference)

    # Select C-alpha atoms
    pred_atoms = [atom for atom in predicted.get_atoms() if atom.get_name() == 'CA']
    exp_atoms = [atom for atom in experimental.get_atoms() if atom.get_name() == 'CA']

    # Ensure equal length
    if len(pred_atoms) != len(exp_atoms):
        raise ValueError('The number of CA atoms in predicted and experimental structures must match.')

    # Superimpose and compute RMSD
    sup = Superimposer()
    sup.set_atoms(exp_atoms, pred_atoms)

    rmsd = sup.rms
    return rmsd

def tm_align(query, reference, name=''):    
    # if not os.path.isfile('TMalign'):
    #     print(subprocess.getoutput('g++ -O3 -ffast-math -lm -o TMalign {s}'.format(s='TMalign.cpp')))
    
    output_dir = f'TMalign/{name}'
    os.makedirs(output_dir, exist_ok=True)

    batcmd = f'./TMalign/TMalign {query} {reference} -o {output_dir}/TM_sup'
    try:
        result = subprocess.check_output(batcmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running TMalign: {e}")
    
    score_lines = []
    with open(f'{output_dir}/results.txt', 'w') as f:
        for line in result.decode().split("\n"):
            f.write(line + "\n")
            if line.startswith("TM-score"):
                score_lines.append(line)

    # Fetch the chain number
    key_getter = lambda s: re.findall(r"Chain_[12]{1}", s)[0]
    score_getter = lambda s: float(re.findall(r"=\s+([0-9.]+)", s)[0])
    results_dict = {key_getter(s): score_getter(s) for s in score_lines}
    
    filename = output_dir + '/TM_sup_all_atm'

    with open(filename, "r") as f:
        pdb_lines = f.readlines()

    pdbfilename = output_dir + f'/result_impose.pdb'
    icounter = 1
    with open(pdbfilename, "w") as f:
        for line in pdb_lines:
            if  icounter >13:
                f.write(line)
            icounter+=1

    return results_dict

"""lDDT protein distance score."""
# import jax.numpy as jnp
import numpy as jnp


def lddt(predicted_points,
         true_points,
         true_points_mask,
         cutoff=15.,
         per_residue=False):
  """Measure (approximate) lDDT for a batch of coordinates.

  lDDT reference:
  Mariani, V., Biasini, M., Barbato, A. & Schwede, T. lDDT: A local
  superposition-free score for comparing protein structures and models using
  distance difference tests. Bioinformatics 29, 2722–2728 (2013).

  lDDT is a measure of the difference between the true distance matrix and the
  distance matrix of the predicted points.  The difference is computed only on
  points closer than cutoff *in the true structure*.

  This function does not compute the exact lDDT value that the original paper
  describes because it does not include terms for physical feasibility
  (e.g. bond length violations). Therefore this is only an approximate
  lDDT score.

  Args:
    predicted_points: (batch, length, 3) array of predicted 3D points
    true_points: (batch, length, 3) array of true 3D points
    true_points_mask: (batch, length, 1) binary-valued float array.  This mask
      should be 1 for points that exist in the true points.
    cutoff: Maximum distance for a pair of points to be included
    per_residue: If true, return score for each residue.  Note that the overall
      lDDT is not exactly the mean of the per_residue lDDT's because some
      residues have more contacts than others.

  Returns:
    An (approximate, see above) lDDT score in the range 0-1.
  """

  assert len(predicted_points.shape) == 3
  assert predicted_points.shape[-1] == 3
  assert true_points_mask.shape[-1] == 1
  assert len(true_points_mask.shape) == 3

  # Compute true and predicted distance matrices.
  dmat_true = jnp.sqrt(1e-10 + jnp.sum(
      (true_points[:, :, None] - true_points[:, None, :])**2, axis=-1))

  dmat_predicted = jnp.sqrt(1e-10 + jnp.sum(
      (predicted_points[:, :, None] -
       predicted_points[:, None, :])**2, axis=-1))

  dists_to_score = (
      (dmat_true < cutoff).astype(jnp.float32) * true_points_mask *
      jnp.transpose(true_points_mask, [0, 2, 1]) *
      (1. - jnp.eye(dmat_true.shape[1]))  # Exclude self-interaction.
  )

  # Shift unscored distances to be far away.
  dist_l1 = jnp.abs(dmat_true - dmat_predicted)

  # True lDDT uses a number of fixed bins.
  # We ignore the physical plausibility correction to lDDT, though.
  score = 0.25 * ((dist_l1 < 0.5).astype(jnp.float32) +
                  (dist_l1 < 1.0).astype(jnp.float32) +
                  (dist_l1 < 2.0).astype(jnp.float32) +
                  (dist_l1 < 4.0).astype(jnp.float32))

  # Normalize over the appropriate axes.
  reduce_axes = (-1,) if per_residue else (-2, -1)
  norm = 1. / (1e-10 + jnp.sum(dists_to_score, axis=reduce_axes))
  score = norm * (1e-10 + jnp.sum(dists_to_score * score, axis=reduce_axes))

  return score[0]

import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import *

def get_coordinates(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    coordinates = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == 'CA':  # Only consider alpha carbons
                        coordinates.append(atom.get_coord())
    
    return np.array([coordinates])

def compute_lddt(ref_pdb, pred_pdb, cutoff=15.0):
    true_points = get_coordinates(ref_pdb)
    predicted_points = get_coordinates(pred_pdb)
    true_points_mask = np.ones(shape=[1, len(true_points), 1])
    return lddt(predicted_points, true_points, true_points_mask)


from Bio.PDB import PDBParser
import numpy as np
from scipy.spatial import cKDTree

def parse_ca_coords(pdb_path):
    """Return sorted list of residue IDs and corresponding Cα coords."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('S', pdb_path)
    coords = []
    res_ids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
                    res_ids.append((chain.id, residue.id[1]))
    return np.array(coords), res_ids

def compute_lddt2(ref_pdb, pred_pdb, cutoff=15.0):
    # load coordinates and residue mappings
    ref_coords, ref_ids = parse_ca_coords(ref_pdb)
    pred_coords, pred_ids = parse_ca_coords(pred_pdb)
    assert ref_ids == pred_ids, "Residue ordering/mapping must match"

    # build a neighbor-search tree on the reference coords
    tree = cKDTree(ref_coords)
    thresholds = [0.5, 1.0, 2.0, 4.0]
    per_residue_scores = []

    for i, r_coord in enumerate(ref_coords):
        # neighbors within cutoff (excluding self)
        neighbors = [j for j in tree.query_ball_point(r_coord, cutoff) if j != i]
        if not neighbors:
            per_residue_scores.append(np.nan)
            continue

        # distance arrays
        d_ref  = np.linalg.norm(ref_coords[neighbors]  - r_coord, axis=1)
        d_pred = np.linalg.norm(pred_coords[neighbors] - pred_coords[i], axis=1)
        diffs  = np.abs(d_pred - d_ref)

        # fraction under each threshold
        fracs = [(diffs < t).sum() / len(diffs) for t in thresholds]
        per_residue_scores.append(np.mean(fracs))

    # global LDDT
    global_lddt = np.nanmean(per_residue_scores)
    return per_residue_scores, global_lddt
    
aa_heavy_atoms = {
    "ALA": ["N", "CA", "C", "O", "CB"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
    "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"]
}

Three_letter_to_one_letter = {
    "ALA": "A", "ARG": "R", "ASN": "N",
    "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V"
}

one_letter_to_three_letter = {v: k for k, v in Three_letter_to_one_letter.items()}

def element(atom):
    if "N" in atom:
        return "N"
    if "O" in atom:
        return "O"
    else:
        return "C"

import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import *

def get_coordinates(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    coordinates = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == 'CA':  # Only consider alpha carbons
                        coordinates.append(atom.get_coord())
    
    return np.array(coordinates)

def calculate_distance_matrix(coordinates):
    n = len(coordinates)
    distance_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            distance_matrix[i, j] = np.linalg.norm(coordinates[i] - coordinates[j])
    
    return distance_matrix

def plot_distance_matrix(pdb_file):
    # Get coordinates from PDB file
    coordinates = get_coordinates(pdb_file)

    # Calculate distance matrix
    distance_matrix = calculate_distance_matrix(coordinates)

    # Plot the distance matrix
    plt.figure(figsize=(10, 8))
    plt.imshow(distance_matrix, cmap='viridis')
    plt.colorbar(label='Distance (Å)')
    plt.title('Pairwise Distance Matrix')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.savefig(f'PD/{pdb_file.split("/")[-1]}_distance_matrix.png')

# clean unnecessary rosetta information from PDB
def clean_pdb(pdb_file):
    # Read the pdb file and filter relevant lines
    with open(pdb_file, 'r') as f_in:
        relevant_lines = [line for line in f_in if line.startswith(('ATOM', 'HETATM', 'MODEL', 'TER', 'END', 'LINK'))]

    # Write the cleaned lines back to the original pdb file
    with open(pdb_file, 'w') as f_out:
        f_out.writelines(relevant_lines)
        
def compute_rmsd_95(query, reference):
    """
    Calculate RMSD95 as implemented in AlphaFold2.
    This calculates RMSD after removing the worst 5% of residues.
    
    Args:
        query (str): Path to query PDB file
        reference (str): Path to reference PDB file
        
    Returns:
        float: RMSD95 value
    """
    # Parse structures
    parser = PDBParser(QUIET=True)
    predicted = parser.get_structure('pred', query)
    experimental = parser.get_structure('exp', reference)

    # Get CA atoms and their coordinates
    pred_atoms = []
    exp_atoms = []
    pred_coords = []
    exp_coords = []
    
    for atom in predicted.get_atoms():
        if atom.get_name() == 'CA':
            pred_atoms.append(atom)
            pred_coords.append(atom.get_coord())
            
    for atom in experimental.get_atoms():
        if atom.get_name() == 'CA':
            exp_atoms.append(atom)
            exp_coords.append(atom.get_coord())
    
    pred_coords = np.array(pred_coords)
    exp_coords = np.array(exp_coords)
    
    # Calculate initial RMSD
    sup = Superimposer()
    sup.set_atoms(exp_atoms, pred_atoms)
    sup.apply(pred_atoms)
    
    # Calculate distances for each residue
    distances = np.linalg.norm(pred_coords - exp_coords, axis=1)
    
    # Sort distances and remove worst 5%
    sorted_indices = np.argsort(distances)
    num_to_remove = int(len(distances) * 0.05)
    keep_indices = sorted_indices[:-num_to_remove]
    
    # Calculate RMSD95
    rmsd_95 = np.sqrt(np.mean(distances[keep_indices]**2))
    
    return rmsd_95


def read_sequence_from_pdb(pdb_file):
    """
    Read a PDB file and return the protein sequence without using any external libraries.
    
    Args:
        pdb_file (str): Path to the PDB file
        
    Returns:
        str: Protein sequence
    """
    # Dictionary mapping three-letter amino acid codes to one-letter codes
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
    }
    
    sequence = ""
    current_chain = None
    current_residue = None
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # Extract chain ID and residue number
                chain = line[21]
                res_num = int(line[22:26])
                res_name = line[17:20].strip()
                
                # Only process if it's a new residue
                if chain != current_chain or res_num != current_residue:
                    if res_name in three_to_one:
                        sequence += three_to_one[res_name]
                    current_chain = chain
                    current_residue = res_num
                    
    return sequence

from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import load

def calculate_similarity(seq1, seq2, matrix_name="BLOSUM62", gap_open=-10, gap_extend=-0.5):
    """
    Calculate percent identity and percent similarity between two protein sequences using a global alignment.

    Parameters:
    - seq1 (str): First amino acid sequence.
    - seq2 (str): Second amino acid sequence.
    - matrix_name (str): Name of substitution matrix to use (default: BLOSUM62).
    - gap_open (float): Gap opening penalty (default: -10).
    - gap_extend (float): Gap extension penalty (default: -0.5).

    Returns:
    - percent_identity (float): Percentage of identical matches in the alignment.
    - percent_similarity (float): Percentage of identical or positively scored substitutions.
    - alignment (tuple): The best alignment (aligned_seq1, aligned_seq2, score, begin, end).
    """
    # Load substitution matrix
    matrix = load(matrix_name)
    
    # Create aligner
    aligner = PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    aligner.mode = 'global'
    
    # Perform alignment
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    
    # Get aligned sequences
    aligned_seq1 = best_alignment[0]
    aligned_seq2 = best_alignment[1]
    
    # Calculate metrics
    identical = 0
    similar = 0
    aligned_length = 0

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == '-' or b == '-':
            continue
        aligned_length += 1
        if a == b:
            identical += 1
        else:
            # Check if substitution is positively scored in the matrix
            if matrix[a][b] > 0:
                similar += 1

    percent_identity = (identical / aligned_length) * 100 if aligned_length > 0 else 0
    percent_similarity = ((identical + similar) / aligned_length) * 100 if aligned_length > 0 else 0

    return percent_identity, percent_similarity, best_alignment


from Bio.PDB import *
from Bio.PDB.Polypeptide import is_aa

Three_letter_to_one_letter = {
    "ALA": "A", "ARG": "R", "ASN": "N",
    "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V"
}

def get_sequence_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    sequence = ""
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue):
                    sequence += Three_letter_to_one_letter[residue.get_resname()]
                    
    return sequence

def read_fasta_sequences(file_path):
    sequences = {}
    current_seq = ""
    current_id = ""
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = current_seq
                current_id = line[1:].split(' ')[0]
                current_seq = ""
            else:
                current_seq += line
        if current_id:
            sequences[current_id] = current_seq
            
    return sequences