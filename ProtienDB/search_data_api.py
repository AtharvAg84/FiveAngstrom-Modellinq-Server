from rcsbapi.search import TextQuery
from rcsbapi.data import DataQuery as Query
import csv
import io

def search_pdb_ids(search_term):
    """
    Search for PDB IDs using a query term.

    Args:
        search_term (str): The search query (e.g., 'Hemoglobin').

    Returns:
        list: List of PDB IDs, or empty list if search fails or no results.
    """
    if not search_term or not search_term.strip():
        return []
    
    try:
        query = TextQuery(value=search_term)
        pdb_ids = list(query())  # Convert iterator to list of PDB IDs
        return pdb_ids
    except Exception as e:
        return []

def query_experimental_data(pdb_ids, limit):
    """
    Query experimental data for a list of PDB IDs, up to a specified limit.

    Args:
        pdb_ids (list): List of PDB IDs to query.
        limit (int): Maximum number of PDB IDs to process.

    Returns:
        tuple: (headers, results) where headers is the list of CSV columns and
               results is a list of dictionaries with experimental data.
    """
    headers = ['rcsb_id', 'method_details', 'details', 'method', 'crystals_number']
    results = []
    
    limited_pdb_ids = pdb_ids[:min(limit, len(pdb_ids))]
    
    for pdb_id in limited_pdb_ids:
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
        except Exception as e:
            results.append({
                'rcsb_id': pdb_id,
                'method_details': None,
                'details': f"Error: {str(e)}",
                'method': None,
                'crystals_number': None
            })
    
    return headers, results

def generate_csv_data(headers, results):
    """
    Generate CSV data as a string from query results, replacing None with '-'.

    Args:
        headers (list): List of CSV column headers.
        results (list): List of dictionaries with row data.

    Returns:
        str: CSV data as a string, or empty string if no results.
    """
    if not results:
        return ""
    
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=headers)
    writer.writeheader()
    for row in results:
        csv_row = {key: '-' if value is None else value for key, value in row.items()}
        writer.writerow(csv_row)
    csv_data = output.getvalue()
    output.close()
    return csv_data