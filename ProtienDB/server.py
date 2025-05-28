import logging
import os
from fastmcp import FastMCP
from search_data_api import search_pdb_ids, query_experimental_data, generate_csv_data
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()
logger.info("Loading environment variables")

# Initialize FastMCP
try:
    mcp = FastMCP("PDB Search Server")
    logger.info("FastMCP server initialized with name: PDB Search Server")
except Exception as e:
    logger.error(f"Failed to initialize FastMCP: {str(e)}")
    raise

# Define MCP tool
@mcp.tool()
def search_pdb_data(query: str, limit: int = 10) -> dict:
    """
    Search for PDB IDs by query term and retrieve experimental data as CSV.

    Args:
        query (str): Search term (e.g., 'Hemoglobin').
        limit (int): Maximum number of PDB IDs to query (default: 10).

    Returns:
        dict: Contains status, CSV data, or error message.
    """
    logger.info(f"Received query: {query}, limit: {limit}")
    try:
        # Validate query
        if not query or len(query.strip()) == 0:
            logger.warning("Empty or invalid query")
            return {
                "status": "error",
                "error": "Query cannot be empty"
            }

        if len(query) > 500:
            logger.warning(f"Query too long: {query[:50]}...")
            return {
                "status": "error",
                "error": "Query exceeds maximum length of 500 characters"
            }

        if limit < 1 or limit > 100:
            logger.warning(f"Invalid limit: {limit}")
            return {
                "status": "error",
                "error": "Limit must be between 1 and 100"
            }

        # Search for PDB IDs
        logger.info(f"Searching PDB IDs for query: {query}")
        pdb_ids = search_pdb_ids(query)
        logger.debug(f"Found {len(pdb_ids)} PDB IDs")

        if not pdb_ids:
            logger.info("No PDB IDs found for query")
            return {
                "status": "success",
                "message": "No PDB IDs found for the query",
                "csv_data": ""
            }

        # Query experimental data
        logger.info(f"Querying experimental data for {len(pdb_ids)} PDB IDs with limit {limit}")
        headers, results = query_experimental_data(pdb_ids, limit)
        logger.debug(f"Retrieved {len(results)} results")

        # Generate CSV data
        logger.info("Generating CSV data")
        csv_data = generate_csv_data(headers, results)
        logger.debug(f"CSV data length: {len(csv_data)} characters")

        return {
            "status": "success",
            "message": f"Found {len(pdb_ids)} PDB IDs, processed {len(results)}",
            "csv_data": csv_data
        }
    except Exception as e:
        logger.error(f"Server error processing query '{query}': {str(e)}")
        return {
            "status": "error",
            "error": f"Server error: {str(e)}"
        }

# Run the server
if __name__ == "__main__":
    try:
        logger.info("Starting FastMCP server with stdio transport")
        mcp.run(transport="stdio")
        logger.info("FastMCP server started successfully")
    except Exception as e:
        logger.error(f"Failed to start server: {str(e)}")
        raise
