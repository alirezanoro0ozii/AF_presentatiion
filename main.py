from utils.utils import *
import warnings
warnings.filterwarnings('ignore')
import csv

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
    gdt = compute_gdt_ts(pdb_name, pdb_name_query)
    return sequence, plddt, results_dict, rmsd, lddt, per_res, g_lddt, rmsd_95, ref_clash_score, query_clash_score, gdt

uid = "T1025"
pdb_name = f"CASP/casp14_targets/{uid}.pdb"
sequence, plddt, results_dict, rmsd, lddt, per_res, g_lddt, rmsd_95, ref_clash_score, query_clash_score, gdt = calculate_metrics(pdb_name)

# Write main metrics to CSV
with open(f'results/{uid}_metrics.csv', 'w', newline='') as csvfile:
    fieldnames = ['metric', 'value']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    writer.writerow({'metric': 'sequence', 'value': sequence})
    writer.writerow({'metric': 'length', 'value': len(sequence)})
    writer.writerow({'metric': 'plddt', 'value': plddt})
    writer.writerow({'metric': 'tm_score', 'value': results_dict})
    writer.writerow({'metric': 'rmsd', 'value': rmsd})
    writer.writerow({'metric': 'lddt', 'value': lddt})
    writer.writerow({'metric': 'Global_LDDT', 'value': g_lddt})
    writer.writerow({'metric': 'RMSD95', 'value': rmsd_95})
    writer.writerow({'metric': 'ref_clash_score', 'value': ref_clash_score})
    writer.writerow({'metric': 'query_clash_score', 'value': query_clash_score})
    writer.writerow({'metric': 'gdt', 'value': gdt})

# Write per-residue LDDT scores to separate CSV
with open(f'results/{uid}_per_residue_lddt.csv', 'w', newline='') as csvfile:
    fieldnames = ['residue', 'lddt_score']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    for idx, score in enumerate(per_res):
        writer.writerow({'residue': idx+1, 'lddt_score': score})