import logging
import json
import os
from fastmcp import FastMCP
import google.generativeai as genai
from fetch_coordinates import get_2d_coordinates, get_3d_coordinates, coords_to_pdb
from dotenv import load_dotenv
from rdkit import Chem

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()
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
    mcp = FastMCP("Chemical Coordinates Server")
    logger.info("FastMCP server initialized with name: Chemical Coordinates Server")
except Exception as e:
    logger.error(f"Failed to initialize FastMCP: {str(e)}")
    raise

# Define MCP tool
@mcp.tool()
def get_coordinates(query: str) -> dict:
    """
    Fetch 2D and 3D coordinates and PDB content for a chemical compound given its SMILES or name.
    Uses Gemini to validate or convert the query to SMILES if needed.

    Args:
        query (str): SMILES string, chemical name, or description.

    Returns:
        dict: 2D and 3D coordinates, PDB contents, Gemini response, or error message.
    """
    logger.info(f"Received query: {query}")
    try:
        # Validate query length
        if len(query) > 1000:
            logger.warning(f"Query too long: {query[:50]}...")
            return {
                "status": "error",
                "error": "Query exceeds maximum length of 1000 characters",
                "gemini_message": ""
            }

        # Use Gemini to process the query
        prompt = f"""
        You are a chemistry expert. If the input is a SMILES string, confirm its validity.
        If the input is a chemical name or description (including misspellings like 'aspirine' for 'aspirin'), convert it to a SMILES string.
        Input: {query}
        Output format: {{ "smiles": "SMILES string", "is_valid": true/false, "message": "explanation" }}
        """
        logger.debug("Sending prompt to Gemini")
        gemini_response = model.generate_content(prompt)

        # Check Gemini response
        if not hasattr(gemini_response, 'text') or not gemini_response.text:
            logger.error("Invalid Gemini response: No text content")
            return {
                "status": "error",
                "error": "Invalid Gemini response: No text content",
                "gemini_message": ""
            }

        # Extract and clean text
        gemini_data = gemini_response.text.strip()
        if gemini_data.startswith("```json"):
            gemini_data = gemini_data[7:].strip()
        if gemini_data.endswith("```"):
            gemini_data = gemini_data[:-3].strip()
        logger.debug(f"Cleaned Gemini response: {gemini_data}")

        # Parse Gemini response
        try:
            gemini_result = json.loads(gemini_data)
            logger.info("Gemini response parsed successfully")
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse Gemini response: {str(e)}")
            return {
                "status": "error",
                "error": "Failed to parse Gemini response",
                "gemini_message": gemini_data
            }

        # Extract SMILES from Gemini
        smiles = gemini_result.get("smiles")
        if not smiles or not gemini_result.get("is_valid"):
            logger.warning(f"Invalid query from Gemini: {gemini_result.get('message')}")
            return {
                "status": "error",
                "error": gemini_result.get("message", "Invalid query from Gemini"),
                "gemini_message": gemini_result.get("message", "")
            }

        # Fetch 2D coordinates using the SMILES
        logger.info(f"Fetching 2D coordinates for SMILES: {smiles}")
        source, coords_2d, mol_2d = get_2d_coordinates(smiles)
        logger.debug(f"get_2d_coordinates returned: source={source}, coords={coords_2d}, mol={mol_2d}")
        if not coords_2d:
            logger.warning("Failed to retrieve or generate 2D coordinates")
            return {
                "status": "error",
                "error": "Failed to retrieve or generate 2D coordinates",
                "gemini_smiles": smiles,
                "gemini_message": gemini_result.get("message", "")
            }

        # Generate 2D PDB content
        query_clean = "".join(c for c in query if c.isalnum() or c in "_-") or "MOLECULE"
        try:
            if mol_2d is None:
                logger.warning("Invalid molecule for 2D PDB generation")
                return {
                    "status": "error",
                    "error": "Invalid molecule for 2D PDB generation",
                    "gemini_smiles": smiles,
                    "gemini_message": gemini_result.get("message", "")
                }
            pdb_content_2d = coords_to_pdb(smiles, coords_2d, mol_2d, query_clean)
            logger.info("2D PDB content generated successfully")
        except Exception as e:
            logger.error(f"Failed to generate 2D PDB content: {str(e)}")
            return {
                "status": "error",
                "error": f"Failed to generate 2D PDB content: {str(e)}",
                "gemini_smiles": smiles,
                "gemini_message": gemini_result.get("message", "")
            }

        # Generate 3D coordinates
        logger.info(f"Generating 3D coordinates for SMILES: {smiles}")
        try:
            mol_3d, coords_3d = get_3d_coordinates(smiles)
            if not coords_3d or mol_3d is None:
                logger.warning("Failed to generate 3D coordinates")
                return {
                    "status": "error",
                    "error": "Failed to generate 3D coordinates",
                    "gemini_smiles": smiles,
                    "gemini_message": gemini_result.get("message", ""),
                    "coordinates_2d": coords_2d,
                    "pdb_content_2d": pdb_content_2d
                }
            logger.info("3D coordinates generated successfully")
        except Exception as e:
            logger.error(f"Failed to generate 3D coordinates: {str(e)}")
            return {
                "status": "error",
                "error": f"Failed to generate 3D coordinates: {str(e)}",
                "gemini_smiles": smiles,
                "gemini_message": gemini_result.get("message", ""),
                "coordinates_2d": coords_2d,
                "pdb_content_2d": pdb_content_2d
            }

        # Generate 3D PDB content
        try:
            # Verify number of atoms matches coordinates
            if len(coords_3d) != mol_3d.GetNumAtoms():
                logger.error(f"Mismatch: {len(coords_3d)} coordinates but {mol_3d.GetNumAtoms()} atoms")
                raise ValueError("Number of 3D coordinates does not match number of atoms")
            pdb_content_3d = coords_to_pdb(smiles, coords_3d, mol_3d, query_clean)
            logger.info("3D PDB content generated successfully")
        except Exception as e:
            logger.error(f"Failed to generate 3D PDB content: {str(e)}")
            return {
                "status": "error",
                "error": f"Failed to generate 3D PDB content: {str(e)}",
                "gemini_smiles": smiles,
                "gemini_message": gemini_result.get("message", ""),
                "coordinates_2d": coords_2d,
                "pdb_content_2d": pdb_content_2d
            }

        return {
            "status": "success",
            "source": source,
            "coordinates_2d": coords_2d,
            "coordinates_3d": coords_3d,
            "gemini_smiles": smiles,
            "gemini_message": gemini_result.get("message", ""),
            "pdb_content_2d": pdb_content_2d,
            "pdb_content_3d": pdb_content_3d
        }
    except Exception as e:
        logger.error(f"Server error processing query '{query}': {str(e)}")
        return {
            "status": "error",
            "error": f"Server error: {str(e)}",
            "gemini_message": ""
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