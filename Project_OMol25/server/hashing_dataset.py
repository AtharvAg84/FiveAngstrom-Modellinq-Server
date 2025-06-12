import psycopg2
import hashlib

# Database connection parameters
db_params = {
    "dbname": "omol4m",
    "user": "atharvag",
    "password": "atharv8484#",
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