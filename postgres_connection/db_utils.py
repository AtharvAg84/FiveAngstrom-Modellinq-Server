import psycopg2
from psycopg2.extras import RealDictCursor
import logging
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Database connection parameters
DB_PARAMS = {
    "dbname": "something",
    "user": "BhaRak",
    "password": "BhaRak2506",
    "host": "localhost",
    "port": "5432"
}

def get_db_connection():
    """
    Establish a connection to the PostgreSQL database.

    Returns:
        connection: psycopg2 connection object, or None if connection fails.
    """
    try:
        conn = psycopg2.connect(**DB_PARAMS)
        logger.info("Connected to PostgreSQL database 'something'")
        return conn
    except psycopg2.Error as e:
        logger.error(f"Failed to connect to database: {str(e)}")
        return None

def get_table_names():
    """
    Fetch the list of table names in the 'public' schema of the database.

    Returns:
        list: List of table names, or empty list if failed.
    """
    conn = get_db_connection()
    if not conn:
        return []

    try:
        with conn.cursor() as cur:
            cur.execute("""
                SELECT table_name 
                FROM information_schema.tables 
                WHERE table_schema = 'public'
            """)
            tables = [row[0] for row in cur.fetchall()]
            logger.info(f"Retrieved {len(tables)} table names")
            return tables
    except psycopg2.Error as e:
        logger.error(f"Failed to fetch table names: {str(e)}")
        return []
    finally:
        conn.close()

def get_table_schema(table_name):
    """
    Fetch the schema (columns and data types) of a specific table.

    Args:
        table_name (str): Name of the table.

    Returns:
        list: List of dictionaries with column_name and data_type, or None if failed.
    """
    conn = get_db_connection()
    if not conn:
        return None

    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute("""
                SELECT column_name, data_type 
                FROM information_schema.columns 
                WHERE table_schema = 'public' AND table_name = %s
            """, (table_name,))
            schema = cur.fetchall()
            if not schema:
                logger.warning(f"No schema found for table: {table_name}")
                return None
            logger.info(f"Retrieved schema for table: {table_name}")
            return schema
    except psycopg2.Error as e:
        logger.error(f"Failed to fetch schema for table {table_name}: {str(e)}")
        return None
    finally:
        conn.close()

def get_row_count(table_name):
    """
    Fetch the number of rows in a specific table.

    Args:
        table_name (str): Name of the table.

    Returns:
        int: Number of rows, or None if failed.
    """
    conn = get_db_connection()
    if not conn:
        return None

    try:
        with conn.cursor() as cur:
            cur.execute(f"SELECT COUNT(*) FROM {table_name}")
            count = cur.fetchone()[0]
            logger.info(f"Retrieved row count for table {table_name}: {count}")
            return count
    except psycopg2.Error as e:
        logger.error(f"Failed to fetch row count for table {table_name}: {str(e)}")
        return None
    finally:
        conn.close()

def get_table_data(table_name, limit=100):
    """
    Fetch data from a specific table, with an optional row limit.

    Args:
        table_name (str): Name of the table.
        limit (int): Maximum number of rows to fetch (default: 100).

    Returns:
        list: List of dictionaries representing rows, or None if failed.
    """
    conn = get_db_connection()
    if not conn:
        return None

    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(f"SELECT * FROM {table_name} LIMIT %s", (limit,))
            data = cur.fetchall()
            logger.info(f"Retrieved {len(data)} rows from table {table_name}")
            return data
    except psycopg2.Error as e:
        logger.error(f"Failed to fetch data from table {table_name}: {str(e)}")
        return None
    finally:
        conn.close()

def get_table_data_filtered(table_name, filter_condition, limit=100):
    """
    Fetch filtered data from a table using a WHERE clause.

    Args:
        table_name (str): Name of the table.
        filter_condition (str): SQL WHERE clause (e.g., "name = 'Alice'").
        limit (int): Maximum number of rows to fetch (default: 100).

    Returns:
        list: List of dictionaries representing rows, or None if failed.
    """
    conn = get_db_connection()
    if not conn:
        return None

    try:
        # Basic validation to prevent SQL injection
        if re.search(r';\s*(drop|alter|delete|update|insert)\s', filter_condition, re.IGNORECASE):
            logger.warning("Unsafe filter condition detected")
            return None

        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            query = f"SELECT * FROM {table_name} WHERE {filter_condition} LIMIT %s"
            cur.execute(query, (limit,))
            data = cur.fetchall()
            logger.info(f"Retrieved {len(data)} rows from table {table_name} with filter")
            return data
    except psycopg2.Error as e:
        logger.error(f"Failed to fetch filtered data from table {table_name}: {str(e)}")
        return None
    finally:
        conn.close()

def get_column_info(table_name, column_name):
    """
    Fetch details about a specific column in a table.

    Args:
        table_name (str): Name of the table.
        column_name (str): Name of the column.

    Returns:
        dict: Column details (column_name, table_name, data_type, is_nullable), or None if failed.
    """
    conn = get_db_connection()
    if not conn:
        return None

    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute("""
                SELECT column_name, table_name, data_type, is_nullable
                FROM information_schema.columns
                WHERE table_schema = 'public' AND table_name = %s AND column_name = %s
            """, (table_name, column_name))
            info = cur.fetchone()
            if not info:
                logger.warning(f"Column {column_name} not found in table {table_name}")
                return None
            logger.info(f"Retrieved info for column {column_name} in table {table_name}")
            return info
    except psycopg2.Error as e:
        logger.error(f"Failed to fetch column info for {column_name} in table {table_name}: {str(e)}")
        return None
    finally:
        conn.close()

def execute_custom_query(sql_query, limit=100):
    """
    Execute a user-provided SELECT query, restricted to SELECT statements.

    Args:
        sql_query (str): SQL SELECT query.
        limit (int): Maximum number of rows to fetch (default: 100).

    Returns:
        list: List of dictionaries representing rows, or None if failed or unsafe.
    """
    conn = get_db_connection()
    if not conn:
        return None

    try:
        # Restrict to SELECT queries and prevent dangerous keywords
        if not sql_query.strip().upper().startswith("SELECT") or re.search(r';\s*(drop|alter|delete|update|insert)\s', sql_query, re.IGNORECASE):
            logger.warning("Unsafe or non-SELECT query detected")
            return None

        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            query = f"{sql_query} LIMIT %s"
            cur.execute(query, (limit,))
            data = cur.fetchall()
            logger.info(f"Executed custom query, retrieved {len(data)} rows")
            return data
    except psycopg2.Error as e:
        logger.error(f"Failed to execute custom query: {str(e)}")
        return None
    finally:
        conn.close()

def get_database_stats():
    """
    Fetch general statistics about the database (e.g., total tables, total rows).

    Returns:
        dict: Statistics dictionary, or empty dict if failed.
    """
    conn = get_db_connection()
    if not conn:
        return {}

    try:
        stats = {}
        with conn.cursor() as cur:
            # Total tables
            cur.execute("""
                SELECT COUNT(*) 
                FROM information_schema.tables 
                WHERE table_schema = 'public'
            """)
            stats["total_tables"] = cur.fetchone()[0]

            # Total rows across all tables
            cur.execute("""
                SELECT table_name 
                FROM information_schema.tables 
                WHERE table_schema = 'public'
            """)
            tables = [row[0] for row in cur.fetchall()]
            total_rows = 0
            for table in tables:
                cur.execute(f"SELECT COUNT(*) FROM {table}")
                total_rows += cur.fetchone()[0]
            stats["total_rows"] = total_rows

            logger.info(f"Retrieved database statistics: {stats}")
            return stats
    except psycopg2.Error as e:
        logger.error(f"Failed to fetch database stats: {str(e)}")
        return {}
    finally:
        conn.close()

def search_tables_columns(keyword):
    """
    Search for tables or columns matching a keyword.

    Args:
        keyword (str): Keyword to search for.

    Returns:
        list: List of dictionaries with type (table/column), table_name, and optional column_name.
    """
    conn = get_db_connection()
    if not conn:
        return []

    try:
        results = []
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            # Search tables
            cur.execute("""
                SELECT table_name 
                FROM information_schema.tables 
                WHERE table_schema = 'public' AND table_name ILIKE %s
            """, (f"%{keyword}%",))
            for row in cur.fetchall():
                results.append({"type": "table", "table_name": row["table_name"]})

            # Search columns
            cur.execute("""
                SELECT table_name, column_name 
                FROM information_schema.columns 
                WHERE table_schema = 'public' AND column_name ILIKE %s
            """, (f"%{keyword}%",))
            for row in cur.fetchall():
                results.append({"type": "column", "table_name": row["table_name"], "column_name": row["column_name"]})

            logger.info(f"Found {len(results)} matches for keyword: {keyword}")
            return results
    except psycopg2.Error as e:
        logger.error(f"Failed to search tables/columns: {str(e)}")
        return []
    finally:
        conn.close()