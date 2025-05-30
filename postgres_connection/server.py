import logging
import json
import os
from fastmcp import FastMCP
import google.generativeai as genai
from db_utils import (
    get_table_names, get_table_schema, get_row_count, get_table_data,
    get_column_info, get_table_data_filtered, execute_custom_query,
    get_database_stats, search_tables_columns
)
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

path="../.env"
# Load environment variables
load_dotenv(path)
logger.info("Loading environment variables")

# Configure Gemini API
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if not GEMINI_API_KEY:
    logger.error("GEMINI_API_KEY not found in .env file")
    raise ValueError("GEMINI_API_KEY not found in .env file")
genai.configure(api_key=GEMINI_API_KEY)
logger.info("Gemini API configured")

# Initialize Gemini model
try:
    model = genai.GenerativeModel("gemini-1.5-flash")
    logger.info("Gemini model initialized")
except Exception as e:
    logger.error(f"Failed to initialize Gemini model: {str(e)}")
    raise

# Initialize FastMCP
try:
    mcp = FastMCP("Database Information Server")
    logger.info("FastMCP server initialized with name: Database Information Server")
except Exception as e:
    logger.error(f"Failed to initialize FastMCP: {str(e)}")
    raise

# Helper function to parse Gemini response
def parse_gemini_response(gemini_response):
    """Parse and clean Gemini JSON response."""
    if not hasattr(gemini_response, 'text') or not gemini_response.text:
        logger.error("Invalid Gemini response: No text content")
        return None
    gemini_data = gemini_response.text.strip()
    if gemini_data.startswith("```json"):
        gemini_data = gemini_data[7:].strip()
    if gemini_data.endswith("```"):
        gemini_data = gemini_data[:-3].strip()
    try:
        return json.loads(gemini_data)
    except json.JSONDecodeError as e:
        logger.error(f"Failed to parse Gemini response: {str(e)}")
        return None

# FastMCP Tools
@mcp.tool()
def get_database_info(query: str) -> dict:
    """
    Route database queries to appropriate tools using Gemini to interpret the query.

    Args:
        query (str): Query about the database (e.g., "list tables", "data in table users where name is Alice").

    Returns:
        dict: Result from the appropriate tool or error message.
    """
    logger.info(f"Received query: {query}")
    try:
        if len(query) > 1000:
            logger.warning(f"Query too long: {query[:50]}...")
            return {
                "status": "error",
                "error": "Query exceeds maximum length of 1000 characters",
                "gemini_message": ""
            }

        # Use Gemini to interpret the query
        prompt = f"""
        You are a database expert. Interpret the user's query about a PostgreSQL database named 'something' and determine which action to take.
        Available actions and their tools:
        - list_tables: List all table names (get_table_names)
        - table_schema: Get schema of a specific table (get_table_schema)
        - row_count: Get number of rows in a specific table (get_row_count)
        - table_data: Get all data from a specific table (get_table_data)
        - table_data_filtered: Get filtered data from a table with conditions (get_table_data_filtered)
        - column_info: Get details of a specific column in a table (get_column_info)
        - custom_query: Execute a user-provided SELECT query (execute_custom_query)
        - database_stats: Get database statistics (get_database_stats)
        - search_tables_columns: Search for tables/columns by keyword (search_tables_columns)
        Input: {query}
        Output format: {{
            "tool": "tool_name",
            "table_name": "table_name or empty",
            "column_name": "column_name or empty",
            "filter_condition": "WHERE clause or empty",
            "sql_query": "SELECT query or empty",
            "keyword": "search keyword or empty",
            "is_valid": true/false,
            "message": "explanation"
        }}
        """
        logger.debug("Sending prompt to Gemini")
        gemini_response = model.generate_content(prompt)
        gemini_result = parse_gemini_response(gemini_response)
        if not gemini_result:
            return {
                "status": "error",
                "error": "Failed to parse Gemini response",
                "gemini_message": ""
            }

        if not gemini_result.get("is_valid"):
            logger.warning(f"Invalid query from Gemini: {gemini_result.get('message')}")
            return {
                "status": "error",
                "error": gemini_result.get("message", "Invalid query from Gemini"),
                "gemini_message": gemini_result.get("message", "")
            }

        # Route to appropriate tool
        tool = gemini_result.get("tool")
        table_name = gemini_result.get("table_name", "")
        column_name = gemini_result.get("column_name", "")
        filter_condition = gemini_result.get("filter_condition", "")
        sql_query = gemini_result.get("sql_query", "")
        keyword = gemini_result.get("keyword", "")
        result = {"status": "success", "gemini_message": gemini_result.get("message", "")}

        if tool == "list_tables":
            logger.info("Fetching list of table names")
            result["tables"] = get_table_names()
        elif tool == "table_schema":
            if not table_name:
                return {"status": "error", "error": "Table name required", "gemini_message": gemini_result.get("message", "")}
            logger.info(f"Fetching schema for table: {table_name}")
            schema = get_table_schema(table_name)
            if not schema:
                return {"status": "error", "error": f"Table {table_name} not found", "gemini_message": gemini_result.get("message", "")}
            result["schema"] = schema
        elif tool == "row_count":
            if not table_name:
                return {"status": "error", "error": "Table name required", "gemini_message": gemini_result.get("message", "")}
            logger.info(f"Fetching row count for table: {table_name}")
            count = get_row_count(table_name)
            if count is None:
                return {"status": "error", "error": f"Table {table_name} not found", "gemini_message": gemini_result.get("message", "")}
            result["row_count"] = count
        elif tool == "table_data":
            if not table_name:
                return {"status": "error", "error": "Table name required", "gemini_message": gemini_result.get("message", "")}
            logger.info(f"Fetching data for table: {table_name}")
            data = get_table_data(table_name)
            if data is None:
                return {"status": "error", "error": f"Table {table_name} not found", "gemini_message": gemini_result.get("message", "")}
            result["data"] = data
        elif tool == "table_data_filtered":
            if not table_name or not filter_condition:
                return {"status": "error", "error": "Table name and filter condition required", "gemini_message": gemini_result.get("message", "")}
            logger.info(f"Fetching filtered data for table: {table_name} with condition: {filter_condition}")
            data = get_table_data_filtered(table_name, filter_condition)
            if data is None:
                return {"status": "error", "error": f"Table {table_name} not found or invalid filter", "gemini_message": gemini_result.get("message", "")}
            result["data"] = data
        elif tool == "column_info":
            if not table_name or not column_name:
                return {"status": "error", "error": "Table and column names required", "gemini_message": gemini_result.get("message", "")}
            logger.info(f"Fetching info for column {column_name} in table: {table_name}")
            info = get_column_info(table_name, column_name)
            if not info:
                return {"status": "error", "error": f"Column {column_name} not found in table {table_name}", "gemini_message": gemini_result.get("message", "")}
            result["column_info"] = info
        elif tool == "custom_query":
            if not sql_query:
                return {"status": "error", "error": "SQL query required", "gemini_message": gemini_result.get("message", "")}
            logger.info(f"Executing custom query: {sql_query}")
            data = execute_custom_query(sql_query)
            if data is None:
                return {"status": "error", "error": "Invalid or unsafe query", "gemini_message": gemini_result.get("message", "")}
            result["data"] = data
        elif tool == "database_stats":
            logger.info("Fetching database statistics")
            result["stats"] = get_database_stats()
        elif tool == "search_tables_columns":
            if not keyword:
                return {"status": "error", "error": "Search keyword required", "gemini_message": gemini_result.get("message", "")}
            logger.info(f"Searching tables/columns with keyword: {keyword}")
            result["search_results"] = search_tables_columns(keyword)
        else:
            logger.warning(f"Unknown tool: {tool}")
            return {
                "status": "error",
                "error": f"Unknown tool: {tool}",
                "gemini_message": gemini_result.get("message", "")
            }

        return result

    except Exception as e:
        logger.error(f"Server error processing query '{query}': {str(e)}")
        return {
            "status": "error",
            "error": f"Server error: {str(e)}",
            "gemini_message": ""
        }

@mcp.tool()
def get_table_data_filtered(table_name: str, filter_condition: str) -> dict:
    """
    Fetch filtered data from a specific table using a WHERE clause.

    Args:
        table_name (str): Name of the table.
        filter_condition (str): SQL WHERE clause (e.g., "name = 'Alice'").

    Returns:
        dict: Filtered data or error message.
    """
    logger.info(f"Fetching filtered data for table: {table_name} with condition: {filter_condition}")
    try:
        data = get_table_data_filtered(table_name, filter_condition)
        if data is None:
            return {"status": "error", "error": f"Table {table_name} not found or invalid filter"}
        return {"status": "success", "data": data}
    except Exception as e:
        logger.error(f"Error fetching filtered data: {str(e)}")
        return {"status": "error", "error": f"Error: {str(e)}"}

@mcp.tool()
def get_column_info(table_name: str, column_name: str) -> dict:
    """
    Fetch details about a specific column in a table.

    Args:
        table_name (str): Name of the table.
        column_name (str): Name of the column.

    Returns:
        dict: Column details or error message.
    """
    logger.info(f"Fetching info for column {column_name} in table: {table_name}")
    try:
        info = get_column_info(table_name, column_name)
        if not info:
            return {"status": "error", "error": f"Column {column_name} not found in table {table_name}"}
        return {"status": "success", "column_info": info}
    except Exception as e:
        logger.error(f"Error fetching column info: {str(e)}")
        return {"status": "error", "error": f"Error: {str(e)}"}

@mcp.tool()
def execute_custom_query(sql_query: str) -> dict:
    """
    Execute a user-provided SELECT query (restricted to SELECT statements).

    Args:
        sql_query (str): SQL SELECT query.

    Returns:
        dict: Query results or error message.
    """
    logger.info(f"Executing custom query: {sql_query}")
    try:
        data = execute_custom_query(sql_query)
        if data is None:
            return {"status": "error", "error": "Invalid or unsafe query"}
        return {"status": "success", "data": data}
    except Exception as e:
        logger.error(f"Error executing custom query: {str(e)}")
        return {"status": "error", "error": f"Error: {str(e)}"}

@mcp.tool()
def get_database_stats() -> dict:
    """
    Fetch general statistics about the database.

    Returns:
        dict: Database statistics or error message.
    """
    logger.info("Fetching database statistics")
    try:
        stats = get_database_stats()
        return {"status": "success", "stats": stats}
    except Exception as e:
        logger.error(f"Error fetching database stats: {str(e)}")
        return {"status": "error", "error": f"Error: {str(e)}"}

@mcp.tool()
def search_tables_columns(keyword: str) -> dict:
    """
    Search for tables or columns matching a keyword.

    Args:
        keyword (str): Keyword to search for.

    Returns:
        dict: Search results or error message.
    """
    logger.info(f"Searching tables/columns with keyword: {keyword}")
    try:
        results = search_tables_columns(keyword)
        return {"status": "success", "search_results": results}
    except Exception as e:
        logger.error(f"Error searching tables/columns: {str(e)}")
        return {"status": "error", "error": f"Error: {str(e)}"}

# Run the server
if __name__ == "__main__":
    try:
        logger.info("Starting FastMCP server with stdio transport")
        mcp.run(transport="stdio")
        logger.info("FastMCP server started successfully")
    except Exception as e:
        logger.error(f"Failed to start server: {str(e)}")
        raise