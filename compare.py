from utils.utils import read_sequence_from_pdb, calculate_similarity

source_seq = read_sequence_from_pdb("source_pdb.pdb")
source_seq = "QVQLVESGGGLVQPGGSLRLSCAASGGSEYSYSTFSLGWFRQAPGQGLEAVAAIASMGGLTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAAVRGYFMRLPSSHNFRYWGQGTLVTVS"

import pandas as pd

df = pd.read_csv("RFAntibody.csv")
df['compare'] = df['sequence'].apply(lambda x: calculate_similarity(x, source_seq)[0])

print(df['compare'])

df.to_csv("RFAntibody_compare.csv", index=False)