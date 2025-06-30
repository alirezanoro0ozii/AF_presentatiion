from utils.utils import *
from utils.mp5 import mp5_fold
import warnings
warnings.filterwarnings('ignore')
import csv
import json

with open('common_sequences.json', 'r') as f:
    common_seq2pdb = json.load(f)

def write_metrics_csv(metrics, sequence, pdb_name_ref, model_name):
    os.makedirs('results', exist_ok=True)
    os.makedirs(f'results/{model_name}', exist_ok=True)
    with open(f'results/{model_name}/{pdb_name_ref.split("/")[-1].split(".")[0]}_metrics.csv', 'w', newline='') as csvfile:
        fieldnames = ['metric', 'value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        writer.writerow({'metric': 'sequence', 'value': sequence})
        writer.writerow({'metric': 'length', 'value': len(sequence)})
        writer.writerow({'metric': 'tm_score', 'value': metrics['tm_score']})
        writer.writerow({'metric': 'rmsd_value', 'value': metrics['rmsd_value']})
        writer.writerow({'metric': 'seq_id', 'value': metrics['seq_id']})
        writer.writerow({'metric': 'CA_rmsd', 'value': metrics['CA_rmsd']})
        writer.writerow({'metric': 'lddt', 'value': metrics['lddt']})
        writer.writerow({'metric': 'gdt', 'value': metrics['gdt']})
        
def calculate_metrics(pdb_name_ref, pdb_name_query, model_name):
    tm_score, rmsd_value, seq_id = tm_align(pdb_name_ref, pdb_name_query, model_name + '_' + pdb_name_ref.split('/')[-1].split('.')[0])
    CA_rmsd = compute_rmsd(pdb_name_ref, pdb_name_query)
    lddt = compute_lddt(pdb_name_ref, pdb_name_query)
    plot_distance_matrix(pdb_name_ref, 'ref')
    plot_distance_matrix(pdb_name_query, model_name, remove_domain=False if '-' in pdb_name_ref.split('/')[-1].split('.')[0] else True)
    plot_distance_matrix_difference(pdb_name_ref, pdb_name_query, model_name)
    gdt = compute_gdt_ts(pdb_name_ref, pdb_name_query)
    return {
        'tm_score': tm_score,
        'rmsd_value': rmsd_value, 
        'seq_id': seq_id,
        'CA_rmsd': CA_rmsd,
        'lddt': lddt,
        'gdt': gdt
    }

for i, (sequence, path_dict) in enumerate(common_seq2pdb.items()):
    try:
        print(f'{i+1}/{len(common_seq2pdb)} {sequence}')
        pdb_name_ref = path_dict['casp14_pdb']

        pdb_name_query_AF2 = path_dict['af2_pdb']
        print(pdb_name_ref, " ", pdb_name_query_AF2)
        pdb_name_query_mp5 = mp5_fold(pdb_name_ref.split('/')[-1].split('.')[0], sequence)

        metrics_mp5 = calculate_metrics(pdb_name_ref, pdb_name_query_mp5, 'MP5')
        metrics_AF2 = calculate_metrics(pdb_name_ref, pdb_name_query_AF2, 'AF2')

        write_metrics_csv(metrics_mp5, sequence, pdb_name_ref, 'MP5')
        write_metrics_csv(metrics_AF2, sequence, pdb_name_ref, 'AF2')

    except Exception as e:
        print(e)
        continue