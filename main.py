from utils.utils import esm_fold, get_sequence_from_pdb, tm_align, compute_rmsd, compute_lddt, compute_lddt2, plot_distance_matrix, compute_rmsd_95

def calculate_metrics(pdb_name):
    # sequence = get_sequence_from_pdb(pdb_name)
    # print("Sequence: ", sequence)
    # pdb_name_esm, plddt = esm_fold(sequence, pdb_name)
    # print("esm_fold: ", pdb_name_esm, "plddt: ", plddt)
    pdb_name_esm = "esm_fold/T1025_esm.pdb"
    pdb_name_tm, tm_score = tm_align(pdb_name, pdb_name_esm, '1025')
    print("tm_align: ", pdb_name_tm, "tm_score: ", tm_score)
    rmsd = compute_rmsd(pdb_name, pdb_name_esm)
    print("rmsd: ", rmsd)
    lddt = compute_lddt(pdb_name, pdb_name_esm)
    print("lddt: ", lddt)
    per_res, g_lddt = compute_lddt2(pdb_name, pdb_name_esm)
    print("per_res: ", per_res)
    print("g_lddt: ", g_lddt)
    plot_distance_matrix(pdb_name)
    plot_distance_matrix(pdb_name_esm)
    print("plot_distance_matrix: ", pdb_name, pdb_name_esm)
    rmsd_95 = compute_rmsd_95(pdb_name, pdb_name_esm)
    return plddt, tm_score, rmsd, lddt, per_res, g_lddt, rmsd_95


pdb_name = "CASP/casp14_targets/T1025.pdb"
plddt, tm_score, rmsd, lddt, per_res, g_lddt, rmsd_95 = calculate_metrics(pdb_name)
print("plddt: ", plddt)
print("tm_score: ", tm_score)
print("rmsd: ", rmsd)
print("lddt: ", lddt)
print(f"Global LDDT = {g_lddt*100:.2f}")
print(f"RMSD95 = {rmsd_95:.2f}")
for idx, score in enumerate(per_res):
    print(f"Res {idx+1:4d}: LDDT = {score*100:5.1f}")