from rcsbapi.data import DataQuery as Query
import csv

# List of PDB IDs
pdb_ids = [
    "1OXR", "1TGM"
]

# Define CSV file headers
headers = ['rcsb_id', 'method_details', 'details', 'method', 'crystals_number']

# List to store query results
results = []

# Query each PDB ID
for pdb_id in pdb_ids:
    print(f"\n--- Querying PDB ID: {pdb_id} ---")
    try:
        query = Query(
            input_type="entries",
            input_ids=[pdb_id],
            return_data_list=["exptl"]
        )
        result_dict = query.exec()
        
        # Extract entry data
        entry = result_dict.get('data', {}).get('entries', [{}])[0]
        exptl = entry.get('exptl', [{}])[0]
        
        # Prepare row data
        row = {
            'rcsb_id': entry.get('rcsb_id', pdb_id),
            'method_details': exptl.get('method_details'),
            'details': exptl.get('details'),
            'method': exptl.get('method'),
            'crystals_number': exptl.get('crystals_number')
        }
        results.append(row)
        print(row)
    except Exception as e:
        print(f"Error querying {pdb_id}: {e}")
        # Add row with error indication
        results.append({
            'rcsb_id': pdb_id,
            'method_details': None,
            'details': f"Error: {str(e)}",
            'method': None,
            'crystals_number': None
        })

# Write results to CSV, replacing None with '-'
with open('pdb_daata.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=headers)
    writer.writeheader()
    for row in results:
        # Create a new row where None is replaced with '-'
        csv_row = {key: '-' if value is None else value for key, value in row.items()}
        writer.writerow(csv_row)

print("\nCSV file 'pdb_data.csv' has been created successfully.")
