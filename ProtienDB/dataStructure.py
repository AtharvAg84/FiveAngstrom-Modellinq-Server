from rcsbapi.data import DataQuery as Query

# List of PDB IDs
pdb_ids = [
    "1OXR", "1TGM", "1PTH", "6MQF", "6UX1", "2QQT", "3GCL", "3KK6", "8J3W",
    "3IAZ", "4NSB", "3N8Y", "5F19", "1CQE", "1PGE", "2I2Z", "1PGF", "1PGG",
    "3JUT", "3K1X", "5FDQ", "6OFY", "5F1A"
]

# Traverse the list and query "exptl" data for each entry
for pdb_id in pdb_ids:
    print(f"Querying PDB ID: {pdb_id}")
    query = Query(
        input_type="entries",
        input_ids=[pdb_id],
        return_data_list=["exptl"]
    )
    result_dict = query.exec()
    exptl_data = result_dict.get("result_set", [{}])[0].get("exptl", [])
    print(f"{pdb_id} - exptl: {exptl_data}\n")
