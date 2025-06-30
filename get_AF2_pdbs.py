domains = """T1024-D1, T1024-D2, T1026-D1, T1027-D1, T1029-D1, T1030-D1, T1030-D2, T1031-D1, T1032-D1,
T1033-D1, T1034-D1, T1035-D1, T1037-D1, T1038-D1, T1038-D2, T1039-D1, T1040-D1, T1041-D1,
T1042-D1, T1043-D1, T1045s2-D1, T1046s1-D1, T1046s2-D1, T1047s1-D1, T1047s2-D1, T1047s2-D2,
T1047s2-D3, T1049-D1, T1050-D1, T1050-D2, T1050-D3, T1052-D1, T1052-D2, T1052-D3, T1053-D1,
T1053-D2, T1054-D1, T1055-D1, T1056-D1, T1057-D1, T1058-D1, T1058-D2, T1060s2-D1, T1060s3-D1,
T1061-D1, T1061-D2, T1061-D3, T1064-D1, T1065s1-D1, T1065s2-D1, T1067-D1, T1068-D1, T1070-D1,
T1070-D2, T1070-D3, T1070-D4, T1073-D1, T1074-D1, T1076-D1, T1078-D1, T1079-D1, T1080-D1,
T1082-D1, T1083-D1, T1084-D1, T1087-D1, T1089-D1, T1090-D1, T1091-D1, T1091-D2, T1091-D3,
T1091-D4, T1092-D1, T1092-D2, T1093-D1, T1093-D2, T1093-D3, T1094-D1, T1094-D2, T1095-D1,
T1096-D1, T1096-D2, T1099-D1, T1100-D1, T1100-D2, T1101-D1, T1101-D2"""
domains =domains.replace("\n", " ")
domains_list = domains.split(", ") 
domains_list = [domain.strip() for domain in domains_list]
print(len(domains_list))
import pandas as pd
df = pd.read_csv("CASP/casp_targets_filtered.csv")
casp_domains = [d.split(":")[0] for d in df['Domain_Info'].values]
print(len(casp_domains))
for domain in casp_domains:
    if domain not in domains_list:
        print(domain)
        
# import os

# os.makedirs("pdbs/AF2_paper", exist_ok=True)
# # Download each model's prediction page using wget
# for domain in domains_list:
#     # https://predictioncenter.org/casp14/results.cgi?target=T1099-D1&model=T1099TS427_1-D1&view=prediction
#     output_file = f"pdbs/AF2_paper/{domain}.pdb"
#     if not os.path.exists(output_file):
#         cmd = f'wget "https://predictioncenter.org/casp14/results.cgi?target={domain}&model={domain.split("-")[0]}TS427_1-{domain.split("-")[1]}&view=prediction" -O "{output_file}"'
#         print(f"Downloading: {cmd}")
#         os.system(cmd)
#     else:
#         print(f"Skipping {domain}.pdb - file already exists")