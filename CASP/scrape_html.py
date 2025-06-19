import re

#                <td title="Model name"  >
                #         <a href="results.cgi?target=T1076-D1&model=T1076TS427_1-D1&view=prediction"> 
                #         <nobr>T1076TS427_1-D1</nobr>
                #         </a>
                # </td>           

def find_model_names(html):
    """
    Find all model name entries in the Results html.
    Returns a list of (model_name, model_href) tuples.
    """
    pattern = re.compile(
        r'<td\s+title="Model name"\s*>\s*'
        r'<a\s+href="([^"]+)">\s*'
        r'<nobr>([^<]+)</nobr>\s*'
        r'</a>\s*'
        r'</td>',
        re.IGNORECASE
    )
    return [(match.group(2), match.group(1)) for match in pattern.finditer(html)]

# Read html and find model names and hrefs
with open("Results - CASP14.html", "r", encoding="utf-8") as f:
    html_content = f.read()
results = [item for item in find_model_names(html_content) if "TS427" in item[0]]
for name, href in results:
    print(f"Model: {name}, Href: {href}")
print(len(results))

import os
os.makedirs("../pdbs/casp14_alphafold", exist_ok=True)
# Download each model's prediction page using wget
for name, href in results:
    # Only download if the href looks like a results.cgi link
    if href.startswith("results.cgi"):
        # Compose the wget command
        # https://predictioncenter.org/casp14/results.cgi?target=T1099-D1&model=T1099TS427_1-D1&view=prediction
        cmd = f'wget "https://predictioncenter.org/casp14/{href}" -O "../pdbs/casp14_alphafold/{name.split("TS")[0]}.pdb"'
        print(f"Downloading: {cmd}")
        os.system(cmd)
