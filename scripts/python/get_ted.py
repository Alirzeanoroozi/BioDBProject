import pandas as pd
import requests
import csv

# Input and output file paths
input_file = "data/part4/idmapping.tsv"
output_file = "data/part6/alphafold_cath.csv"

# Read the tsv file
df = pd.read_csv(input_file, sep="\t")

results = []

for idx, row in df.iterrows():
    alphafold_id = row['AlphaFoldDB']
    if pd.isnull(alphafold_id):
        print(f"Warning: No data for {row['Entry']}")
        continue
    alphafold_id = alphafold_id.replace(";", "")

    # Construct API URL for the Alphafold (should be Uniprot ID)
    url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{alphafold_id}?skip=0&limit=100"
    try:
        response = requests.get(url, headers={"accept": "application/json"})
        if response.status_code == 200:
            data = response.json()
            # This assumes response contains a 'cath' or similar field; adjust if the key is different
            cath_info = data.get("data", [])
            # If multiple CATH domains are present, store as list or semicolon joined

            if isinstance(cath_info, list):
                for cath_dom in cath_info:
                    result = {"Alphafold": alphafold_id}
                    result.update(cath_dom) 
                    results.append(result)
        else:
            print(f"Warning: No data for {alphafold_id} (status code {response.status_code})")
    except Exception as e:
        print(f"Error fetching CATH for {alphafold_id}: {e}")

# Write out the CSV file
with open(output_file, "w", newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=results[0].keys())
    writer.writeheader()
    for row in results:
        writer.writerow(row)