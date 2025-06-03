import psycopg2
import json
import numpy as np
from fairchem.core.datasets import AseDBDataset

# Database connection parameters
db_params = {
    "dbname": "omol4m",
    "user": "atharvag",  # PostgreSQL username
    "password": "atharv8484#",  # PostgreSQL password
    "host": "localhost",
    "port": "5432"
}

# Dataset path
dataset_path = "./train_4M/"
dataset = AseDBDataset({"src": dataset_path})

# Connect to PostgreSQL
conn = psycopg2.connect(**db_params)
cur = conn.cursor()

# SQL insert query
insert_query = """
INSERT INTO info (
    index, source, reference_source, data_id, charge, spin, num_atoms, num_electrons,
    num_ecp_electrons, n_scf_steps, n_basis, unrestricted, nl_energy, integrated_densities,
    homo_energy, homo_lumo_gap, s_squared, s_squared_dev, warnings, mulliken_charges,
    lowdin_charges, nbo_charges, composition, mulliken_spins, lowdin_spins, nbo_spins,
    chemical_formula, chemical_symbol, coordinates
) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
"""

# Iterate over dataset
for i in range(len(dataset)):
    atoms = dataset.get_atoms(i)
    info = atoms.info

    # Helper function to convert NumPy arrays to lists
    def to_json_serializable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (list, dict)):
            return obj
        if isinstance(obj, (int, float, str, bool)) or obj is None:
            return obj
        return str(obj)  # Fallback for unexpected types

    # Prepare data for insertion
    data = (
        i,  # index
        info.get("source"),
        info.get("reference_source"),
        info.get("data_id"),
        info.get("charge"),
        info.get("spin"),
        info.get("num_atoms"),
        info.get("num_electrons"),
        info.get("num_ecp_electrons"),
        info.get("n_scf_steps"),
        info.get("n_basis"),
        info.get("unrestricted"),
        info.get("nl_energy"),
        json.dumps(to_json_serializable(info.get("integrated_densities"))),
        json.dumps(to_json_serializable(info.get("homo_energy"))),
        json.dumps(to_json_serializable(info.get("homo_lumo_gap"))),
        info.get("s_squared"),
        info.get("s_squared_dev"),
        json.dumps(to_json_serializable(info.get("warnings"))),
        json.dumps(to_json_serializable(info.get("mulliken_charges"))),
        json.dumps(to_json_serializable(info.get("lowdin_charges"))),
        json.dumps(to_json_serializable(info.get("nbo_charges"))),
        info.get("composition"),
        json.dumps(to_json_serializable(info.get("mulliken_spins"))) if info.get("unrestricted") else None,
        json.dumps(to_json_serializable(info.get("lowdin_spins"))) if info.get("unrestricted") else None,
        json.dumps(to_json_serializable(info.get("nbo_spins"))) if info.get("unrestricted") else None,
        atoms.get_chemical_formula(),  # chemical_formula
        json.dumps(atoms.get_chemical_symbols()),  # chemical_symbol
        json.dumps(to_json_serializable(atoms.positions))  # coordinates
    )

    # Execute insert
    try:
        cur.execute(insert_query, data)
        print(f"Processed {i+1} records", end="\r")  # Erase print for record count
    except Exception as e:
        print(f"\nError inserting record at index {i}: {e}")  # Print error on new line
        conn.rollback()
        continue

# Commit and close
conn.commit()
cur.close()
conn.close()

print("\nData inserted successfully!")  # Final message on new line