#   <tr>
# 	<td  style="width: 35px;">1.</td>
# 	<td>T1104</a></td>
# 	<td>117</td>
#         <td>All groups</td>
#         <td>T1104-D1: 1-117 </td>
# 	<td><font >117</font></td>
# 	<td>FM</td>
# 	<!---td>2022-05-02</td--->
# 	<td><a  href="http://www.pdb.org/pdb/explore/explore.do?structureId=7roa">7roa</a></td>

# </tr>

import csv
import re
from bs4 import BeautifulSoup

def extract_casp_data_from_html(html_file_path, csv_file_path):
    """
    Extract CASP target data from HTML file and save to CSV
    """
    # Read the HTML file
    with open(html_file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()
    
    # Parse HTML with BeautifulSoup
    soup = BeautifulSoup(html_content, 'html.parser')
    
    # Find all table rows that contain CASP target data
    # Looking for rows with target IDs like T1104, T1026, etc.
    rows = soup.find_all('tr')
    
    casp_data = []
    
    for row in rows:
        # Find all cells in the row
        cells = row.find_all('td')
        
        if len(cells) >= 8:  # Ensure we have enough cells
            # Extract data from each cell
            try:
                rank = cells[0].get_text(strip=True)
                target_id = cells[1].get_text(strip=True)
                length = cells[2].get_text(strip=True)
                groups = cells[3].get_text(strip=True)
                domain_info = cells[4].get_text(strip=True)
                domain_length = cells[5].get_text(strip=True)
                category = cells[6].get_text(strip=True)
                
                # Extract PDB ID from the link if present
                pdb_link = cells[7].find('a')
                pdb_id = pdb_link.get_text(strip=True) if pdb_link else ""
                
                # Only add if we have a valid target ID (starts with T and has numbers)
                if re.match(r'^T\d+', target_id):
                    casp_data.append([
                        rank,
                        target_id,
                        length,
                        groups,
                        domain_info,
                        domain_length,
                        category,
                        pdb_id
                    ])
                    
            except (IndexError, AttributeError):
                continue
    
    # Write to CSV file
    with open(csv_file_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header
        writer.writerow([
            'Rank',
            'Target_ID',
            'Length',
            'Groups',
            'Domain_Info',
            'Domain_Length',
            'Category',
            'PDB_ID'
        ])
        
        # Write data
        writer.writerows(casp_data)
    
    print(f"Extracted {len(casp_data)} CASP targets to {csv_file_path}")
    return casp_data

# Example usage
if __name__ == "__main__":
    html_file = "CASP_page.html"
    csv_file = "casp_targets.csv"
    
    data = extract_casp_data_from_html(html_file, csv_file)
    
    # Print first few entries as preview
    print("\nFirst 5 entries:")
    for entry in data[:5]:
        print(entry)
        
        # Read the CSV file and filter out unwanted entries
        import pandas as pd
        
        # Read the CSV file
        df = pd.read_csv(csv_file)
        
        # Remove entries where Groups column contains "Server only"
        df = df[~df['Groups'].str.contains('Server only', na=False)]
        
        # Remove entries where Category column contains "MultiDom"
        df = df[~df['Category'].str.contains('MultiDom', na=False)]
        
        # Write the filtered data back to CSV
        df.to_csv("casp_targets_filtered.csv", index=False)
        
        print(f"Filtered CSV saved. Remaining entries: {len(df)}")
