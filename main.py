from utils.utils import *

def calculate_metrics(pdb_name):
    sequence = get_sequence_from_pdb(pdb_name)
    pdb_name_query, plddt = fold(sequence, pdb_name.split('/')[-1].split('.')[0])
    results_dict = tm_align(pdb_name, pdb_name_query, pdb_name.split('/')[-1].split('.')[0])
    rmsd = compute_rmsd(pdb_name, pdb_name_query)
    lddt = compute_lddt(pdb_name, pdb_name_query)
    per_res, g_lddt = compute_lddt2(pdb_name, pdb_name_query)
    plot_distance_matrix(pdb_name, 'ref')
    plot_distance_matrix(pdb_name_query, 'query')
    plot_distance_matrix_difference(pdb_name, pdb_name_query)
    rmsd_95 = compute_rmsd_95(pdb_name, pdb_name_query)
    ref_clash_score = calculate_clash_score(pdb_name)
    query_clash_score = calculate_clash_score(pdb_name_query)
    return sequence, plddt, results_dict, rmsd, lddt, per_res, g_lddt, rmsd_95, ref_clash_score, query_clash_score

pdb_name = "CASP/casp14_targets/T1025.pdb"
sequence, plddt, results_dict, rmsd, lddt, per_res, g_lddt, rmsd_95, ref_clash_score, query_clash_score = calculate_metrics(pdb_name)
print("sequence: ", sequence)
print("plddt: ", plddt)
print("tm_score: ", results_dict)
print("rmsd: ", rmsd)
print("lddt: ", lddt)
print("Global LDDT: ", g_lddt)
print("RMSD95: ", rmsd_95)
print("ref_clash_score: ", ref_clash_score)
print("query_clash_score: ", query_clash_score)
# for idx, score in enumerate(per_res):
#     print("Res ", idx+1, ": LDDT = ", score)