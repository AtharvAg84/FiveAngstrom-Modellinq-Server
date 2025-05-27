import asyncio
import json
import time
import logging
import os
from fetch_coordinates import coords_to_pdb

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

GEMINI_API_KEY = os.environ.get("GOOGLE_API_KEY")

# Server details
SERVER_SCRIPT = "server.py"
TOOL_NAME = "get_coordinates"

async def call_mcp_server(query: str, retries=3) -> dict:
    """
    Send a query to the FastMCP server using the stdio transport.

    Args:
        query (str): SMILES string, chemical name, or description.
        retries (int): Number of retry attempts for connection failures.

    Returns:
        dict: Server response containing coordinates, Gemini data, or error.
    """
    logger.info(f"Attempting to connect to server via stdio for query: {query}")
    for attempt in range(retries):
        try:
            # Initialize FastMCP client
            from fastmcp.client import Client
            client = Client(SERVER_SCRIPT)
            logger.info(f"Connected to server via stdio transport: {client.transport}")

            # Use async context manager to connect to the server
            async with client:
                # Call the get_coordinates tool with query as a dict
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
                if hasattr(response, 'text'):  # Check for TextContent-like object
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
    print("Welcome to the FastMCP Chemistry Coordinates Client!")
    print("Enter a SMILES string, chemical name, or description to get 2D coordinates in PDB format.")
    print(f"Connecting to server via stdio: {SERVER_SCRIPT}")
    print("Type 'exit' to quit.")

    while True:
        query = input("\nYou: ").strip()
        if query.lower() == "exit":
            print("Goodbye!")
            break

        logger.info(f"Processing user query: {query}")
        result = await call_mcp_server(query)
        print("Response from server:")
        logger.debug(f"Result: {result}")

        if result.get("status") == "success":
            # Generate PDB content
            smiles = result.get("gemini_smiles")
            coords = result.get("coordinates")
            source = result.get("source")
            # Get RDKit molecule for bonds
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print("Error: Failed to parse SMILES for PDB generation.")
                continue

            # Generate PDB
            query_clean = "".join(c for c in query if c.isalnum() or c in "_-")
            name = query_clean or "MOLECULE"
            pdb_content = coords_to_pdb(smiles, coords, mol, name)

            # Save to file
            output_file = f"{query_clean or 'molecule'}_output.pdb"
            try:
                with open(output_file, "w") as f:
                    f.write(pdb_content)
                print(f"PDB file saved: {output_file}")
                print(f"SMILES: {smiles}")
                print(f"Source: {source}")
                print(f"Message: {result.get('gemini_message')}")
            except Exception as e:
                print(f"Error saving PDB file: {str(e)}")
        else:
            print(f"Error: {result.get('error')}")
            if result.get("gemini_message"):
                print(f"Gemini Message: {result.get('gemini_message')}")

if __name__ == "__main__":
    asyncio.run(main())