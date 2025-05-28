from rcsbapi.search import TextQuery
from rcsbapi.data import DataQuery as Query
import csv

# Step 1: Get user input for search query
search_term = input("Enter search query for PDB structures (e.g., Hemoglobin): ").strip()
if not search_term:
    print("Error: Search query cannot be empty.")
    exit(1)

# Step 2: Search for PDB IDs
print(f"\n--- Searching for structures associated with '{search_term}' ---")
try:
    query = TextQuery(value=search_term)
    pdb_ids = list(query())  # Convert iterator to list of PDB IDs
    if not pdb_ids:
        print("No PDB IDs found for the search query.")
        exit(1)
    print(f"Found {len(pdb_ids)} PDB IDs:")
    for i, pdb_id in enumerate(pdb_ids, 1):
        print(f"{i}. {pdb_id}")
except Exception as e:
    print(f"Error during search: {e}")
    exit(1)

# Step 3: Ask user to limit the number of PDB IDs
try:
    limit = input(f"\nEnter the number of PDB IDs to process (max {len(pdb_ids)}): ").strip()
    limit = int(limit)
    if limit <= 0 or limit > len(pdb_ids):
        print(f"Invalid input. Must be between 1 and {len(pdb_ids)}.")
        exit(1)
except ValueError:
    print("Invalid input. Please enter a valid number.")
    exit(1)

# Limit the PDB IDs list
limited_pdb_ids = pdb_ids[:limit]
print(f"\nProcessing {limit} PDB IDs: {limited_pdb_ids}")

# Step 4: Query experimental data for each PDB ID
headers = ['rcsb_id', 'method_details', 'details', 'method', 'crystals_number']
results = []

for pdb_id in limited_pdb_ids:
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

# Step 5: Write results to CSV, replacing None with '-'
if results:
    with open('pdb_data.csv', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for row in results:
            # Replace None with '-' in CSV output
            csv_row = {key: '-' if value is None else value for key, value in row.items()}
            writer.writerow(csv_row)
    print("\nCSV file 'pdb_data.csv' has been created successfully.")
else:
    print("\nNo results to write to CSV. Check API connectivity or PDB IDs.")
