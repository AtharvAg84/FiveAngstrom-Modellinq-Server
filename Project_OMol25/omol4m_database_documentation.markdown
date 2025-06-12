# OMol4M PostgreSQL Database Documentation

This document provides a complete guide to setting up the `omol4m` PostgreSQL database, creating the `info` table, managing the `atharvag` user, and migrating the database to another PostgreSQL instance. It includes all necessary commands, file names, and steps to replicate the setup, intended for handover to the project manager.

## Project Overview
The `omol4m` database stores chemical structure data from the `AseDBDataset` (part of the `fairchem.core.datasets` library), specifically for a dataset located at `./train_4M/`. The `info` table contains metadata and properties of molecular structures, such as composition, atomic positions, and DFT-calculated properties (e.g., HOMO energy, charges). A `composition_hash` column stores a numeric (64-bit integer) SHA-256 hash of the `composition` column for fast searching.

The setup includes:
- A PostgreSQL user (`atharvag`) to own and manage the database.
- The `omol4m` database.
- The `info` table with 30 columns, including `id`, `index`, `composition`, and `composition_hash`.
- Python scripts to populate the table and compute hashes.
- Migration instructions to transfer the database to another PostgreSQL instance.

## Files
The following scripts are provided (all should be saved in the project directory):
- `setup_omol4m.sql`: Creates the user, database, and table.
- `populate_omol4m_info.py`: Populates the `info` table from the dataset.
- `add_composition_hash_numeric.py`: Adds and populates the `composition_hash` column.
- `migrate_omol4m.sh`: Shell script for database migration.

## Prerequisites
- **PostgreSQL**: Version 9.6 or later (for `SHA256` function support, if used in triggers).
- **Python**: Version 3.11 or later, with required packages:
  ```bash
  pip install psycopg2-binary fairchem-core ase numpy
  ```
- **Dataset**: The `AseDBDataset` files at `./train_4M/` (update path if different).
- **Access**: A superuser account (e.g., `postgres`) to create users and databases.

## Database Setup

### Step 1: Download the train_4M dataset
https://huggingface.co/facebook/OMol25/blob/main/DATASET.md

### Step 2: Create the PostgreSQL User
Create the `atharvag` user to own the database and table.

**File**: `setup_omol4m.sql` (partial)

**Command**:
```sql
-- Run as postgres user
CREATE USER atharvag WITH PASSWORD 'atharv**';
```
- **Note**: Replace `'atharv**'` with a secure password. Update `db_params` in Python scripts if changed.

**Execution**:
```bash
psql -U postgres -c "CREATE USER atharvag WITH PASSWORD 'atharv**';"
```

### Step 3: Create the Database
Create the `omol4m` database and grant ownership to `atharvag`.

**File**: `setup_omol4m.sql` (partial)

**Command**:
```sql
-- Run as postgres user
CREATE DATABASE omol4m OWNER atharvag;
```
- **Note**: This creates the database with `atharvag` as the owner, avoiding permission issues.

**Execution**:
```bash
psql -U postgres -c "CREATE DATABASE omol4m OWNER atharvag;"
```

### Step 4: Switch to the `atharvag` User
Connect to PostgreSQL as `atharvag` to perform subsequent operations.

**Command**:
```bash
psql -U atharvag -d omol4m
```
- **Alternative**: Run SQL scripts as `atharvag`:
  ```bash
  psql -U atharvag -d omol4m -f setup_omol4m.sql
  ```

### Step 5: Create the `info` Table
Create the `info` table in the `omol4m` database with 30 columns to store molecular data and the numeric hash.

**File**: `setup_omol4m.sql` (complete)

**Content**:
```sql
-- setup_omol4m.sql
-- Run as postgres user to create user and database
CREATE USER atharvag WITH PASSWORD 'atharv**';
CREATE DATABASE omol4m OWNER atharvag;

-- Connect to omol4m as atharvag
\c omol4m atharvag

-- Create info table
CREATE TABLE info (
    id SERIAL PRIMARY KEY,
    index INTEGER NOT NULL,
    source TEXT NOT NULL,
    reference_source TEXT,
    data_id VARCHAR(50) NOT NULL,
    charge INTEGER NOT NULL,
    spin INTEGER NOT NULL,
    num_atoms INTEGER NOT NULL,
    num_electrons INTEGER NOT NULL,
    num_ecp_electrons INTEGER NOT NULL,
    n_scf_steps INTEGER NOT NULL,
    n_basis INTEGER NOT NULL,
    unrestricted BOOLEAN NOT NULL,
    nl_energy FLOAT,
    integrated_densities JSONB,
    homo_energy JSONB,
    homo_lumo_gap JSONB,
    s_squared FLOAT,
    s_squared_dev FLOAT,
    warnings JSONB,
    mulliken_charges JSONB,
    lowdin_charges JSONB,
    nbo_charges JSONB,
    composition TEXT NOT NULL,
    mulliken_spins JSONB,
    lowdin_spins JSONB,
    nbo_spins JSONB,
    chemical_formula TEXT NOT NULL,
    chemical_symbol JSONB NOT NULL,
    coordinates JSONB NOT NULL,
    composition_hash BIGINT
);

-- Grant privileges to atharvag
GRANT ALL ON TABLE info TO atharvag;
GRANT USAGE, SELECT ON SEQUENCE info_id_seq TO atharvag;
```

**Table Schema**:
| Column Name           | Type       | Description                                                                 |
|-----------------------|------------|-----------------------------------------------------------------------------|
| `id`                  | SERIAL     | Auto-incrementing primary key (1, 2, 3, ...).                               |
| `index`               | INTEGER    | Dataset index from `dataset.get_atoms(i)`.                                  |
| `source`              | TEXT       | Unique identifier (e.g., file path like `omol/electrolytes/...`).           |
| `reference_source`    | TEXT       | Internal identifier (can be NULL).                                          |
| `data_id`             | VARCHAR(50)| Dataset domain (e.g., `elytes`).                                            |
| `charge`              | INTEGER    | Total charge (e.g., `0`).                                                   |
| `spin`                | INTEGER    | Total spin (e.g., `1`).                                                     |
| `num_atoms`           | INTEGER    | Number of atoms (e.g., `84`).                                               |
| `num_electrons`       | INTEGER    | Number of electrons (e.g., `396`).                                           |
| `num_ecp_electrons`   | INTEGER    | Effective core potential electrons (e.g., `0`).                              |
| `n_scf_steps`         | INTEGER    | Number of DFT steps (e.g., `16`).                                           |
| `n_basis`             | INTEGER    | Number of basis functions (e.g., `2177`).                                    |
| `unrestricted`        | BOOLEAN    | Restricted/unrestricted flag (e.g., `False`).                                |
| `nl_energy`           | FLOAT      | Dispersion energy (e.g., `41.255288292679985`).                             |
| `integrated_densities`| JSONB      | Electron density integrals (e.g., `[197.99992624, 197.99992624, ...]`).     |
| `homo_energy`         | JSONB      | HOMO energy (e.g., `[-8.86718388]`).                                        |
| `homo_lumo_gap`       | JSONB      | HOMO-LUMO gap (e.g., `[8.31819417]`).                                       |
| `s_squared`           | FLOAT      | Net magnetization (e.g., `0.0`).                                            |
| `s_squared_dev`       | FLOAT      | S^2 deviation (e.g., `0.0`).                                                |
| `warnings`            | JSONB      | DFT warning messages (e.g., `[]`).                                           |
| `mulliken_charges`    | JSONB      | Partial Mulliken charges (array, `natoms x 1`).                              |
| `lowdin_charges`      | JSONB      | Partial Löwdin charges (array, `natoms x 1`).                                |
| `nbo_charges`         | JSONB      | Partial NBO charges (array, `natoms x 1`, can be NULL).                      |
| `composition`         | TEXT       | Chemical composition (e.g., `B1Br1C27H36N2O16S1`).                           |
| `mulliken_spins`      | JSONB      | Partial Mulliken spins (array, `natoms x 1`, NULL if `unrestricted` False).  |
| `lowdin_spins`        | JSONB      | Partial Löwdin spins (array, `natoms x 1`, NULL if `unrestricted` False).    |
| `nbo_spins`           | JSONB      | Partial NBO spins (array, `natoms x 1`, NULL if `unrestricted` False).       |
| `chemical_formula`    | TEXT       | Formula from `atoms.get_chemical_formula()` (e.g., `C27H36N2O16S1B1Br1`).    |
| `chemical_symbol`     | JSONB      | Chemical symbols array (e.g., `["C", "H", ...]`).                           |
| `coordinates`         | JSONB      | Atomic positions (e.g., `[[x1, y1, z1], [x2, y2, z2], ...]`).               |
| `composition_hash`    | BIGINT     | 64-bit integer hash of `composition` (first 8 bytes of SHA-256, masked).     |

**Execution**:
```bash
psql -U atharvag -d omol4m -f setup_omol4m.sql
```
- **Note**: If run as `postgres`, the `CREATE USER` and `CREATE DATABASE` commands will execute, but switch to `atharvag` for the table creation.

### Step 6: Populate the `info` Table
Use the Python script to populate the `info` table from the `AseDBDataset`.

**File**: `populate_omol4m_info.py`

**Content**:
```python
import psycopg2
import json
import numpy as np
from fairchem.core.datasets import AseDBDataset

# Database connection parameters
db_params = {
    "dbname": "omol4m",
    "user": "atharvag",
    "password": "atharv**",
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
        atoms.get_chemical_formula(),
        json.dumps(atoms.get_chemical_symbols()),
        json.dumps(to_json_serializable(atoms.positions))
    )

    # Execute insert
    try:
        cur.execute(insert_query, data)
        print(f"Processed {i+1} records", end="\r")
    except Exception as e:
        print(f"\nError inserting record at index {i}: {e}")
        conn.rollback()
        continue

# Commit and close
conn.commit()
cur.close()
conn.close()

print("\nData inserted successfully!")
```

**Execution**:
```bash
python populate_omol4m_info.py
```
- **Notes**:
  - Ensure `./train_4M/` exists and is accessible. Update `dataset_path` if needed (e.g., `/home/atharvag/train_4M/`).
  - The script inserts rows without `composition_hash` (added in the next step).
  - For a 4M-row dataset, consider batch inserts for performance (contact developer for modification).

### Step 7: Add Numeric Composition Hash
Add the `composition_hash` column (`BIGINT`) and populate it with a 64-bit integer derived from the SHA-256 hash of `composition`.

**File**: `add_composition_hash_numeric.py`

**Content**:
```python
import psycopg2
import hashlib

# Database connection parameters
db_params = {
    "dbname": "omol4m",
    "user": "atharvag",
    "password": "atharv**",
    "host": "localhost",
    "port": "5432"
}

try:
    # Connect to PostgreSQL
    conn = psycopg2.connect(**db_params)
    cur = conn.cursor()

    # Add composition_hash column as BIGINT (if not exists)
    cur.execute("""
        ALTER TABLE info
        ADD COLUMN IF NOT EXISTS composition_hash BIGINT;
    """)

    # Fetch all composition values
    cur.execute("SELECT id, composition FROM info")
    rows = cur.fetchall()

    # Update each row with numeric SHA-256 hash (first 8 bytes as 64-bit integer)
    for i, (id_val, composition) in enumerate(rows, 1):
        if composition:
            # Compute SHA-256 hash and take first 8 bytes as integer
            hash_bytes = hashlib.sha256(composition.encode('utf-8')).digest()
            hash_int = int.from_bytes(hash_bytes[:8], 'big') & 0x7FFFFFFFFFFFFFFF  # Mask to fit signed BIGINT
            # Update the row
            cur.execute("""
                UPDATE info
                SET composition_hash = %s
                WHERE id = %s
            """, (hash_int, id_val))
        else:
            # Handle NULL or empty composition
            cur.execute("""
                UPDATE info
                SET composition_hash = NULL
                WHERE id = %s
            """, (id_val,))

        # Erase print for progress
        print(f"Processed {i} records", end="\r")

    # Commit changes
    conn.commit()
    print("\nNumeric composition hashes added successfully!")

except Exception as e:
    print(f"\nError: {e}")
    conn.rollback()

finally:
    # Close cursor and connection
    cur.close()
    conn.close()
```

**Execution**:
```bash
python add_composition_hash_numeric.py
```
- **Notes**:
  - Adds `composition_hash` and computes a 64-bit integer hash (first 8 bytes of SHA-256, masked to fit `BIGINT`’s signed range).
  - Uses erase print for progress (`Processed {i} records`).
  - For large datasets, batch updates can improve performance (contact developer).

### Step 8: Optimize Searching
Create an index on `composition_hash` for fast searches.

**Command**:
```sql
\c omol4m
CREATE INDEX idx_composition_hash ON info (composition_hash);
```

**Execution**:
```bash
psql -U atharvag -d omol4m -c "CREATE INDEX idx_composition_hash ON info (composition_hash);"
```
- **Note**: This speeds up queries like `SELECT * FROM info WHERE composition_hash = 1234567890;`.

### Step 9: Verify the Setup
Check the table and data:
```sql
\c omol4m
SELECT COUNT(*) FROM info;
SELECT composition, composition_hash FROM info LIMIT 5;
```
- **Expected**: Row count matches dataset size; `composition_hash` shows 64-bit integers.

## Database Migration
To migrate the `omol4m` database to another PostgreSQL instance (e.g., from a development server to production), use `pg_dump` and `pg_restore`.

### Prerequisites
- **Source Server**: The current PostgreSQL instance with `omol4m`.
- **Destination Server**: A PostgreSQL instance (same or newer version).
- **Access**: Superuser credentials (`postgres`) on both servers.
- **Network**: Connectivity between servers (or file transfer for dump file).

### Migration Steps
1. **Dump the Database (Source Server)**:
   Create a dump file of `omol4m`.

   **File**: `migrate_omol4m.sh` (partial)

   **Command**:
   ```bash
   pg_dump -U postgres -Fc omol4m > omol4m_backup.dump
   ```
   - **Options**:
     - `-Fc`: Custom format for efficient restore.
     - `omol4m_backup.dump`: Output file (transfer to destination server).
   - **Execution**:
     ```bash
     pg_dump -U postgres -Fc omol4m > omol4m_backup.dump
     ```

2. **Transfer the Dump File**:
   Copy `omol4m_backup.dump` to the destination server (e.g., via `scp`):
   ```bash
   scp omol4m_backup.dump user@destination_host:/path/to/destination/
   ```

3. **Create User and Database (Destination Server)**:
   Set up `atharvag` and `omol4m` on the destination server.

   **File**: `migrate_omol4m.sh` (partial)

   **Commands**:
   ```bash
   psql -U postgres -c "CREATE USER atharvag WITH PASSWORD 'atharv**';"
   psql -U postgres -c "CREATE DATABASE omol4m OWNER atharvag;"
   ```

4. **Restore the Database**:
   Restore the dump file to the new `omol4m` database.

   **Command**:
   ```bash
   pg_restore -U postgres -d omol4m --verbose omol4m_backup.dump
   ```
   - **Options**:
     - `-d omol4m`: Target database.
     - `--verbose`: Detailed output for debugging.
   - **Execution**:
     ```bash
     pg_restore -U postgres -d omol4m --verbose /path/to/omol4m_backup.dump
     ```

5. **Verify Migration**:
   Connect to the new database and check:
   ```sql
   \c omol4m
   SELECT COUNT(*) FROM info;
   SELECT composition, composition_hash FROM info LIMIT 5;
   ```

**Complete Migration Script**:
**File**: `migrate_omol4m.sh`
```bash
#!/bin/bash

# Source server: Dump the database
pg_dump -U postgres -Fc omol4m > omol4m_backup.dump

# Transfer dump file (update user@destination_host and paths)
scp omol4m_backup.dump user@destination_host:/path/to/destination/

# Destination server: Create user and database
psql -U postgres -c "CREATE USER atharvag WITH PASSWORD 'atharv**';"
psql -U postgres -c "CREATE DATABASE omol4m OWNER atharvag;"

# Destination server: Restore the database
pg_restore -U postgres -d omol4m --verbose /path/to/omol4m_backup.dump

# Clean up
rm omol4m_backup.dump
```
**Execution**:
```bash
chmod +x migrate_omol4m.sh
./migrate_omol4m.sh
```
- **Notes**:
  - Update `user@destination_host` and paths in `scp`.
  - Ensure `postgres` credentials are valid on both servers.
  - If `atharvag` already exists on the destination, skip the `CREATE USER` command.

### Troubleshooting
- **Permission Denied**:
  - Error: `permission denied for table info`.
  - Fix:
    ```sql
    \c omol4m
    GRANT ALL ON TABLE info TO atharvag;
    GRANT USAGE, SELECT ON SEQUENCE info_id_seq TO atharvag;
    ```
  - Check owner:
    ```sql
    SELECT tableowner FROM pg_tables WHERE tablename = 'info';
    ALTER TABLE info OWNER TO atharvag;
    ```
- **BIGINT Out of Range**:
  - Fixed in `add_composition_hash_numeric.py` by masking (`& 0x7FFFFFFFFFFFFFFF`).
- **Dataset Path**:
  - Ensure `./train_4M/` exists. Update `dataset_path` in `populate_omol4m_info.py` if needed.
- **Slow Performance**:
  - For 4M rows, batch inserts/updates can be implemented. Contact developer for modified scripts.
- **Migration Issues**:
  - Check PostgreSQL versions (`SELECT version();`). Destination must be same or newer.
  - Verify dump file integrity:
    ```bash
    pg_restore -l omol4m_backup.dump
    ```

## Handover Checklist
- **Files**:
  - `setup_omol4m.sql`
  - `populate_omol4m_info.py`
  - `add_composition_hash_numeric.py`
  - `migrate_omol4m.sh`
- **Dependencies**:
  - PostgreSQL, Python, dataset libraries, `./train_4M/`.
- **Credentials**:
  - User: `atharvag`, password: `atharv**` (update securely).
  - Superuser: `postgres` (provide password separately).
- **Steps**:
  1. Run `setup_omol4m.sql`.
  2. Run `populate_omol4m_info.py`.
  3. Run `add_composition_hash_numeric.py`.
  4. Create index on `composition_hash`.
  5. Use `migrate_omol4m.sh` for migration.
- **Verification**:
  - Check row count and sample data.
  - Test search performance: `SELECT * FROM info WHERE composition_hash = ...;`.

## Contact
For issues, contact the developer or database administrator with:
- PostgreSQL logs: `/var/log/postgresql/`.
- Error messages and script outputs.
- Dataset and server details.

**Prepared by**: Atharv Agarwal
**Date**: June 12, 2025