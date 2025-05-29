import asyncio
import json
import time
import logging
import os

# Configure logging to file and console for debugging
logging.basicConfig(
    level=logging.DEBUG,  # Increased to DEBUG for more detail
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)

# Server details
SERVER_SCRIPT = "pdb_search_server.py"
TOOL_NAME = "search_pdb_data"

async def call_mcp_server(query: str, limit: int = 10, retries: int = 3) -> dict:
    """
    Send a query to the FastMCP server.

    Args:
        query (str): Search term for PDB database (e.g., 'Hemoglobin').
        limit (int): Maximum number of PDB IDs to process (default: 10).
        retries (int): Number of retry attempts for connection failures.

    Returns:
        dict: Server response containing PDB IDs, experimental data, CSV, or error.
    """
    logger.debug(f"Attempting to connect to server for query: {query}, limit: {limit}")
    for attempt in range(retries):
        try:
            # Initialize FastMCP client
            from fastmcp.client import Client
            client = Client(SERVER_SCRIPT)
            logger.info(f"Connected to server via stdio transport: {client.transport}")
            # Use async context manager to connect to the server
            async with client:
                # Call the search_pdb_data tool with query as a dict
                logger.debug(f"Calling tool {TOOL_NAME} with query: {query}")
                response = await client.call_tool(TOOL_NAME, {"query": query})
                logger.info("Received response from server")
                logger.debug(f"Raw response: {response}")
                # Handle response type
                if isinstance(response, list):
                    logger.debug(f"Response is a list: {response}")
                    if len(response) == 0:
                        return {
                            "status": "error",
                            "error": "Empty response list from server"
                        }
                    response = response[0]  # Take the first element
                    logger.info("Extracted first element from response list")
                # Handle TextContent object
                if hasattr(response, 'text'):
                    logger.debug("Response is a TextContent object")
                    try:
                        response_text = response.text.strip()
                        # Remove code block markers if present
                        if response_text.startswith("```json"):
                            response_text = response_text[7:].strip()
                        if response_text.endswith("```"):
                            response_text = response_text[:-3].strip()
                        response = json.loads(response_text)
                        logger.info("Parsed TextContent text to dictionary")
                    except json.JSONDecodeError as e:
                        logger.error(f"Failed to parse TextContent text: {str(e)}")
                        return {
                            "status": "error",
                            "error": f"Failed to parse TextContent text: {str(e)}"
                        }
                # Ensure response is a dictionary
                if not isinstance(response, dict):
                    logger.error(f"Unexpected response type: {type(response)}")
                    return {
                        "status": "error",
                        "error": f"Unexpected response type: {type(response)}"
                    }
                return response
        except Exception as e:
            logger.error(f"Connection failed: {str(e)}")
            if attempt == retries - 1:
                return {
                    "status": "error",
                    "error": f"Failed to connect to server via stdio after {retries} attempts. Ensure server.py is accessible and running. Error: {str(e)}"
                }
            logger.info(f"Retrying ({attempt + 1}/{retries})...")
            time.sleep(1)

async def main():
    print("Welcome to the FastMCP PDB Search Client!")
    print("Enter a search term (e.g., 'Hemoglobin') to find PDB IDs and experimental data.")
    print(f"Connecting to server: {SERVER_SCRIPT} (Transport: {os.environ.get('MCP_TRANSPORT', 'stdio')})")
    print("Optionally specify a limit (e.g., 'Hemoglobin 5' for max 5 PDB IDs).")
    print("Type 'exit' to quit.")

    while True:
        user_input = input("\nYou: ").strip()
        if user_input.lower() == "exit":
            print("Goodbye!")
            break
        # Parse query and limit
        parts = user_input.split()
        query = " ".join(parts[:-1]) if parts and parts[-1].isdigit() else user_input
        limit = int(parts[-1]) if parts and parts[-1].isdigit() else 10
        logger.debug(f"Processing user query: {query}, limit: {limit}")
        result = await call_mcp_server(query, limit)
        print("Response from server:")
        logger.debug(f"Result: {result}")
        if result.get("status") == "success":
            # Extract response data
            search_term = result.get("gemini_search_term")
            pdb_ids = result.get("pdb_ids")
            experimental_data = result.get("experimental_data")
            csv_data = result.get("csv_data")
            # Save CSV to file
            query_clean = "".join(c for c in query if c.isalnum() or c in "_-") or "search"
            output_file = f"{query_clean}_output.csv"
            try:
                with open(output_file, "w", encoding="utf-8") as f:
                    f.write(csv_data)
                print(f"CSV file saved: {output_file}")
            except Exception as e:
                print(f"Error saving CSV file: {str(e)}")
            # Display results
            print(f"Search Term: {search_term}")
            print(f"PDB IDs ({len(pdb_ids)}): {', '.join(pdb_ids)}")
            print(f"Experimental Data ({len(experimental_data)} entries):")
            for entry in experimental_data:
                print(f"- {entry['rcsb_id']}: Method={entry['method'] or '-'}, Details={entry['details'] or '-'}")
            print(f"Message: {result.get('gemini_message')}")
        else:
            print(f"Error: {result.get('error')}")
            if result.get("gemini_message"):
                print(f"Message: {result.get('gemini_message')}")

if __name__ == "__main__":
    asyncio.run(main())