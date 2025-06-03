import asyncio
import json
import time
import logging
import os
from fastmcp.client import Client
from typing import Dict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MoleculeClient:
    """FastMCP client for searching molecules by composition."""
    
    def __init__(self, server_script: str, tool_name: str, retries: int = 3):
        """
        Initialize client with server and tool details.
        
        Args:
            server_script (str): Path to server script (e.g., 'server.py').
            tool_name (str): Name of the server tool (e.g., 'search_molecule').
            retries (int): Number of retry attempts.
        """
        self.server_script = server_script
        self.tool_name = tool_name
        self.retries = retries

    async def call_server(self, composition: str) -> Dict:
        """
        Call the server with a composition query.
        
        Args:
            composition (str): Chemical composition or query.
            
        Returns:
            Dict: Server response or error.
        """
        logger.info(f"Sending query: {composition[:50]}{'...' if len(composition) > 50 else ''}")
        for attempt in range(self.retries):
            try:
                client = Client(self.server_script)
                logger.info(f"Connected to server: {client.transport}")
                async with client:
                    response = await client.call_tool(self.tool_name, {"composition": composition})
                    logger.info("Received response")

                    if isinstance(response, list):
                        if len(response) == 0:
                            return {"status": "error", "error": "Empty response list"}
                        response = response[0]

                    if hasattr(response, 'text'):
                        text = response.text.strip()
                        if text.startswith("```json"):
                            text = text[7:].strip()
                        if text.endswith("```"):
                            text = text[:-3].strip()
                        response = json.loads(text)

                    if not isinstance(response, dict):
                        logger.error(f"Unexpected response type: {type(response)}")
                        return {"status": "error", "error": f"Unexpected response type: {type(response)}"}
                    return response

            except Exception as e:
                logger.error(f"Connection failed: {str(e)}")
                if attempt == self.retries - 1:
                    return {
                        "status": "error",
                        "error": f"Failed to connect after {self.retries} attempts. Ensure server is running."
                    }
                logger.info(f"Retrying ({attempt + 1}/{self.retries})...")
                await asyncio.sleep(1)

    async def run(self) -> None:
        """Run the interactive client loop."""
        print("Welcome to the Molecule Search Client!")
        print("Enter a chemical composition (e.g., 'C9H8O4') or name (e.g., 'aspirin') for XYZ and PDB data.")
        print("For general queries (e.g., 'hello guys'), get a chat response.")
        print(f"Connecting to server: {self.server_script}")
        print("Type 'exit' to quit.")

        while True:
            composition = input("\nYou: ").strip()
            if composition.lower() == "exit":
                print("Goodbye!")
                break

            result = await self.call_server(composition)
            print("Response from server:")
            if result.get("status") == "success":
                if "chat_message" in result:
                    print(f"Chat Response: {result['chat_message']}")
                else:
                    molecule_data = result.get("molecule_data")
                    xyz_content = result.get("xyz_content")
                    pdb_content = result.get("pdb_content")
                    gemini_composition = result.get("gemini_composition")
                    gemini_message = result.get("gemini_message")
                    xyz_filename = result.get("xyz_filename", f"{gemini_composition}_molecule.xyz")
                    pdb_filename = result.get("pdb_filename", f"{gemini_composition}_molecule.pdb")

                    try:
                        with open(xyz_filename, "w") as f:
                            f.write(xyz_content)
                        with open(pdb_filename, "w") as f:
                            f.write(pdb_content)
                        print(f"XYZ file saved: {xyz_filename}")
                        print(f"PDB file saved: {pdb_filename}")
                        print(f"Composition: {gemini_composition}")
                        print(f"Chemical Formula: {molecule_data.get('chemical_formula')}")
                        print(f"Number of Atoms: {molecule_data.get('num_atoms')}")
                        print(f"Message: {gemini_message}")
                    except Exception as e:
                        print(f"Error saving files: {str(e)}")
            else:
                print(f"Error: {result.get('error')}")
                if result.get("gemini_message"):
                    print(f"Gemini Message: {result.get('gemini_message')}")

if __name__ == "__main__":
    client = MoleculeClient(server_script="server.py", tool_name="search_molecule")
    asyncio.run(client.run())