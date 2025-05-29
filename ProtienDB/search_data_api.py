import csv
import io

class PDBSearchAPI:
    def search_pdb_ids(self, search_term, limit=None):
        from rcsbapi.search import TextQuery
        if not search_term or not search_term.strip():
            return []
        try:
            # Pass the limit directly to TextQuery if supported
            query = TextQuery(value=search_term, max_results=limit) if limit is not None else TextQuery(value=search_term)
            pdb_ids = list(query())
            return pdb_ids
        except Exception:
            return []

    def query_experimental_data(self, pdb_ids, limit):
        from rcsbapi.data import DataQuery as Query
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
                entry = result_dict.get('data', {}).get('entries', [{}])[0]
                exptl = entry.get('exptl', [{}])[0]
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

    def generate_csv_data(self, headers, results):
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