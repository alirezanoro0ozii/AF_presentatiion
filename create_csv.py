import os
import pandas as pd
import os
import boto3
from utils.utils import get_sequence_from_pdb

FOLDER_NAME = "results/example_outputs"
CSV_FILE_NAME = "RFAntibody.csv"

AWS_ACCESS_KEY_ID="AKIA5KLWPDHX6CINHWXU"
AWS_SECRET_ACCESS_KEY="xjzpbxRpVeeE94+xQqSHhKFO7vXNvuQO+dVoKza3"
S3_BUCKET_NAME="model-server-data"
AWS_REGION="us-west-2"

# Initialize S3 client
s3 = boto3.client(
    's3',
    region_name=AWS_REGION,
    aws_access_key_id=AWS_ACCESS_KEY_ID,
    aws_secret_access_key=AWS_SECRET_ACCESS_KEY
)

def add_to_csv(file_name):
    sequence = get_sequence_from_pdb(f"{FOLDER_NAME}/{file_name}")

    with open(f"{FOLDER_NAME}/{file_name}", 'r') as f:
        lines = f.readlines()
        # SCORE interaction_pae: 3.99
        # SCORE pae: 7.14
        # SCORE pred_lddt: 0.92
        # SCORE target_aligned_antibody_rmsd: 15.11
        # SCORE target_aligned_cdr_rmsd: 8.16
        # SCORE framework_aligned_antibody_rmsd: 0.84
        # SCORE framework_aligned_cdr_rmsd: 1.67
        # SCORE framework_aligned_H1_rmsd: 1.12
        # SCORE framework_aligned_H2_rmsd: 0.51
        # SCORE framework_aligned_H3_rmsd: 2.32
        # SCORE framework_aligned_L1_rmsd: nan
        # SCORE framework_aligned_L2_rmsd: nan
        # SCORE framework_aligned_L3_rmsd: nan
        last_12_lines = lines[-12:]
    
    # Read existing CSV
    try:
        df = pd.read_csv(CSV_FILE_NAME)
    except:
        df = pd.DataFrame(columns=['file_name', 'url', 'sequence'] + [line.split(":")[0].split(" ")[1].strip() for line in last_12_lines])
    
    # Upload to S3
    s3.put_object(
        Bucket=S3_BUCKET_NAME,
        Key=file_name,
        Body=open(f"{FOLDER_NAME}/{file_name}", 'rb'),
        ACL='public-read'
    )

    # Generate public URL
    url = f"https://{S3_BUCKET_NAME}.s3.{AWS_REGION}.amazonaws.com/{file_name}"
    
    new_row = [file_name.split("_")[3], url, sequence] + [float(line.split(":")[1].strip()) for line in last_12_lines]
    df.loc[len(df)] = new_row

    # Save back to CSV
    df.to_csv(CSV_FILE_NAME, index=False)

for files in sorted(os.listdir(FOLDER_NAME)):
    if files.endswith(".pdb"):
        add_to_csv(files)

df = pd.read_csv(CSV_FILE_NAME)
print(df.iloc[0])