import requests
import re
import os
import subprocess
from Bio.PDB import PDBParser, Superimposer
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import load
from Bio.PDB.Polypeptide import is_aa

# detect C alpha clashes for deformed trajectories
def calculate_clash_score(pdb_file, threshold=2.4, only_ca=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    atoms = []
    atom_info = []  # Detailed atom info for debugging and processing

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'H':  # Skip hydrogen atoms
                        continue
                    if only_ca and atom.get_name() != 'CA':
                        continue
                    atoms.append(atom.coord)
                    atom_info.append((chain.id, residue.id[1], atom.get_name(), atom.coord))

    tree = cKDTree(atoms)
    pairs = tree.query_pairs(threshold)

    valid_pairs = set()
    for (i, j) in pairs:
        chain_i, res_i, name_i, coord_i = atom_info[i]
        chain_j, res_j, name_j, coord_j = atom_info[j]

        # Exclude clashes within the same residue
        if chain_i == chain_j and res_i == res_j:
            continue

        # Exclude directly sequential residues in the same chain for all atoms
        if chain_i == chain_j and abs(res_i - res_j) == 1:
            continue

        # If calculating sidechain clashes, only consider clashes between different chains
        if not only_ca and chain_i == chain_j:
            continue

        valid_pairs.add((i, j))

    return len(valid_pairs)

def calculate_pdb_average(pdb_file):
    """
    Calculate average of the last column (B-factor/occupancy) in PDB file
    
    Args:
        pdb_file (str): Path to PDB file
        
    Returns:
        float: Average value of the last column
    """
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    total = 0
    count = 0
    for line in lines:
        if line.startswith('ATOM'):
            try:
                # The last column is typically at position 60-66
                value = float(line[60:66].strip())
                total += value
                count += 1
            except ValueError:
                continue
    
    return total / count if count > 0 else 0

def fold(seq, name):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    payload = seq
    headers = {
        'Content-Type': 'text/plain'
    }

    response = requests.request("POST", url, headers=headers, data=payload)
    with open(f'pdbs/esm_fold/{name}.pdb', 'w') as f:
        f.write(response.text)
    
    avg_value = calculate_pdb_average(f'pdbs/esm_fold/{name}.pdb')
    return f'pdbs/esm_fold/{name}.pdb', avg_value
    
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
    
    # Extract aligned length, RMSD, and sequence identity from score_lines
    aligned_info = None
    for line in result.decode().split("\n"):
        if line.startswith("Aligned length="):
            aligned_info = line
            break
    
    if aligned_info:
        # Parse aligned length, RMSD, and sequence identity
        rmsd_value = float(re.findall(r"RMSD=\s+([0-9.]+)", aligned_info)[0])
        seq_id = float(re.findall(r"Seq_ID=n_identical/n_aligned=\s+([0-9.]+)", aligned_info)[0])
    
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

    return results_dict['Chain_1'], rmsd_value, seq_id

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
  dmat_true = np.sqrt(1e-10 + np.sum(
      (true_points[:, :, None] - true_points[:, None, :])**2, axis=-1))

  dmat_predicted = np.sqrt(1e-10 + np.sum(
      (predicted_points[:, :, None] -
       predicted_points[:, None, :])**2, axis=-1))

  dists_to_score = (
      (dmat_true < cutoff).astype(np.float32) * true_points_mask *
      np.transpose(true_points_mask, [0, 2, 1]) *
      (1. - np.eye(dmat_true.shape[1]))  # Exclude self-interaction.
  )

  # Shift unscored distances to be far away.
  dist_l1 = np.abs(dmat_true - dmat_predicted)

  # True lDDT uses a number of fixed bins.
  # We ignore the physical plausibility correction to lDDT, though.
  score = 0.25 * ((dist_l1 < 0.5).astype(np.float32) +
                  (dist_l1 < 1.0).astype(np.float32) +
                  (dist_l1 < 2.0).astype(np.float32) +
                  (dist_l1 < 4.0).astype(np.float32))

  # Normalize over the appropriate axes.
  reduce_axes = (-1,) if per_residue else (-2, -1)
  norm = 1. / (1e-10 + np.sum(dists_to_score, axis=reduce_axes))
  score = norm * (1e-10 + np.sum(dists_to_score * score, axis=reduce_axes))

  return score[0]

def compute_lddt(ref_pdb, pred_pdb):
    true_points = get_coordinates(ref_pdb)
    predicted_points = get_coordinates(pred_pdb)
    true_points_mask = np.ones(shape=[1, len(true_points), 1])
    return lddt(predicted_points, true_points, true_points_mask)

def compute_lddt2(ref_pdb, pred_pdb, cutoff=15.0):
    ref_coords = get_coordinates(ref_pdb).squeeze(0)
    pred_coords = get_coordinates(pred_pdb).squeeze(0)

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

def calculate_distance_matrix(coordinates):
    coordinates = coordinates.squeeze(0)
    n = len(coordinates)
    distance_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            distance_matrix[i, j] = np.linalg.norm(coordinates[i] - coordinates[j])
    
    return distance_matrix

def plot_distance_matrix(pdb_file, name, remove_domain=False):
    if remove_domain:
        os.makedirs(f'PD/{pdb_file.split("/")[-1].split(".")[0].split("-")[0]}', exist_ok=True)
    else:
        os.makedirs(f'PD/{pdb_file.split("/")[-1].split(".")[0]}', exist_ok=True)
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
    if remove_domain:
        plt.savefig(f'PD/{pdb_file.split("/")[-1].split(".")[0].split("-")[0]}/{name}_distance_matrix.png')
    else:
        plt.savefig(f'PD/{pdb_file.split("/")[-1].split(".")[0]}/{name}_distance_matrix.png')
    plt.close()

def plot_distance_matrix_difference(pdb_file_ref, pdb_file_query, model_name):
    coordinates_ref = get_coordinates(pdb_file_ref)
    coordinates_query = get_coordinates(pdb_file_query)
    distance_matrix_ref = calculate_distance_matrix(coordinates_ref)
    distance_matrix_query = calculate_distance_matrix(coordinates_query)
    distance_matrix_difference = np.abs(distance_matrix_ref - distance_matrix_query)
    plt.figure(figsize=(10, 8))
    plt.imshow(distance_matrix_difference, cmap='viridis')
    plt.colorbar(label='Distance (Å)')
    plt.title('Pairwise Distance Matrix Difference')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.savefig(f'PD/{pdb_file_ref.split("/")[-1].split(".")[0]}/{model_name}_distance_matrix_difference.png')
    plt.close()

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

def parse_ca_coordinates(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    coordinates = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == 'CA':  # Only consider alpha carbons
                        coordinates.append((residue.get_id()[1], atom.get_coord()))
    
    return coordinates

def kabsch(P, Q):
    """
    Kabsch algorithm: given two (N×3) point sets P and Q (already centered),
    returns rotation matrix U that minimizes RMSD(P, Q U).
    """
    C = P.T @ Q
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1,1,d])
    U = V @ D @ Wt
    return U

def compute_gdt_ts(pdb_model, pdb_target, cutoffs=(1,2,4,8), max_iter=10):
    """
    Compute GDT-TS between model and target PDBs.
    
    Returns: GDT_TS float in [0,100].
    """
    # 1) Parse CA coords and build matching lists
    m_ca = parse_ca_coordinates(pdb_model)
    t_ca = parse_ca_coordinates(pdb_target)
    # build dict for target by residue_id
    t_dict = {rid: coord for rid, coord in t_ca}
    # keep only residues present in both
    pairs = [(m_coord, t_dict[rid]) for rid, m_coord in m_ca if rid in t_dict]
    if not pairs:
        raise ValueError("No matching Cα residues found.")
    P = np.array([p for p, _ in pairs])
    Q = np.array([q for _, q in pairs])
    L = len(P)

    # 2) center both on their centroids
    P_centroid = P.mean(axis=0)
    Q_centroid = Q.mean(axis=0)
    P -= P_centroid
    Q -= Q_centroid

    scores = []
    for d in cutoffs:
        P_sub, Q_sub = P.copy(), Q.copy()
        # iterative selection & superposition
        for _ in range(max_iter):
            # find rotation on current subset
            U = kabsch(P_sub, Q_sub)
            P_rot = P @ U
            # compute distances to target
            dists = np.linalg.norm(P_rot - Q, axis=1)
            # select only those within cutoff
            mask = dists <= d
            newP = P[mask]
            newQ = Q[mask]
            # if no change or empty, stop
            if len(newP)==len(P_sub) or len(newP)==0:
                break
            P_sub, Q_sub = newP, newQ

        Nd = len(P_sub)
        scores.append(100.0 * Nd / L)

    # 3) average the four percentages
    gdt_ts = float(np.mean(scores))
    return gdt_ts

# import os

# def check_files_in_folder(folder1, folder2):
#     """
#     Check if all files in folder1 are present in folder2.
#     Returns a list of files that are missing in folder2.
#     """
#     files1 = set(os.listdir(folder1))
#     files2 = set(os.listdir(folder2))
#     # Only consider files, not directories
#     files1 = {f for f in files1 if os.path.isfile(os.path.join(folder1, f))}
#     files2 = {f for f in files2 if os.path.isfile(os.path.join(folder2, f))}
#     missing_files = files1 - files2
#     common_files = files1 & files2
#     if not missing_files:
#         print("All files in '{}' are present in '{}'.".format(folder1, folder2))
#     else:
#         print("Missing files in '{}':".format(folder2))
#         for f in missing_files:
#             print(f)
#     print("Common files:", len(common_files))
#     return missing_files, common_files

# missing_files, common_files = check_files_in_folder('pdbs/casp14_alphafold', 'CASP/casp14_targets')
# print(missing_files)
# print(common_files)