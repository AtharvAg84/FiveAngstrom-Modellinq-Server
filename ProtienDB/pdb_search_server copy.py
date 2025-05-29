import json
import os
from fastmcp import FastMCP
import google.generativeai as genai
from search_data_api import PDBSearchAPI
from dotenv import load_dotenv
import base64
from io import BytesIO

path="../.env"
# Load environment variables
load_dotenv(path)

# Configure Gemini API
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if not GEMINI_API_KEY:
    raise ValueError("GEMINI_API_KEY not found in .env file")
genai.configure(api_key=GEMINI_API_KEY)

# Initialize Gemini model
try:
    model = genai.GenerativeModel("gemini-1.5-flash")
except Exception as e:
    raise

# Initialize FastMCP
try:
    mcp = FastMCP("PDB Search Server")
except Exception as e:
    raise

# Initialize PDBSearchAPI
pdb_api = PDBSearchAPI()

# Define MCP tool
@mcp.tool()
def search_pdb_data(query: str, limit: int = 10) -> dict:
    """
    Search for PDB IDs and retrieve experimental data for a given query term.
    Uses Gemini to validate or refine the query (e.g., correct misspellings).
    
    Args:
        query (str): Search term (e.g., 'Hemoglobin', 'Insulin').
        limit (int): Maximum number of PDB IDs to process (default: 10).
    
    Returns:
        dict: PDB IDs, experimental data, CSV content, Gemini response, or error message.
    """
    try:
        # Validate query length and limit
        if len(query) > 1000:
            return {
                "status": "error",
                "error": "Query exceeds maximum length of 1000 characters",
                "gemini_message": ""
            }

        # Use Gemini to process the query
        prompt = f"""
        You are a biochemistry expert. Validate the input query as a valid search term for PDB database searches.
        If the input contains misspellings (e.g., 'Hemoglobine' for 'Hemoglobin') or ambiguous terms, suggest a corrected or clarified term.
        Input: {query}
        Output format: {{ "search_term": "corrected term", "is_valid": true/false, "message": "explanation" }}
        """
        gemini_response = model.generate_content(prompt)

        # Check Gemini response
        if not hasattr(gemini_response, 'text') or not gemini_response.text:
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

        # Parse Gemini response
        try:
            gemini_result = json.loads(gemini_data)
        except json.JSONDecodeError as e:
            return {
                "status": "error",
                "error": "Failed to parse Gemini response",
                "gemini_message": gemini_data
            }

        # Extract search term from Gemini
        search_term = gemini_result.get("search_term")
        if not search_term or not gemini_result.get("is_valid"):
            return {
                "status": "error",
                "error": gemini_result.get("message", "Invalid query from Gemini"),
                "gemini_message": gemini_result.get("message", "")
            }

        # Search for PDB IDs using the class instance
        pdb_ids = pdb_api.search_pdb_ids(search_term)
        pdb_ids= pdb_ids[:limit] if limit else pdb_ids
        if not pdb_ids:
            return {
                "status": "error",
                "error": "No PDB IDs found for the search term",
                "gemini_search_term": search_term,
                "gemini_message": gemini_result.get("message", "")
            }

        # Query experimental data using the class instance
        headers, results = pdb_api.query_experimental_data(pdb_ids, limit)
        if not results:
            return {
                "status": "error",
                "error": "No experimental data retrieved",
                "gemini_search_term": search_term,
                "gemini_message": gemini_result.get("message", ""),
                "pdb_ids": pdb_ids
            }

        # Generate CSV data using the class instance
        try:
            csv_data = pdb_api.generate_csv_data(headers, results)
            if not csv_data:
                return {
                    "status": "error",
                    "error": "Failed to generate CSV data",
                    "gemini_search_term": search_term,
                    "gemini_message": gemini_result.get("message", ""),
                    "pdb_ids": pdb_ids,
                    "experimental_data": results
                }
        except Exception as e:
            return {
                "status": "error",
                "error": f"Failed to generate CSV data: {str(e)}",
                "gemini_search_term": search_term,
                "gemini_message": gemini_result.get("message", ""),
                "pdb_ids": pdb_ids,
                "experimental_data": results
            }

        return {
            "status": "success",
            "gemini_search_term": search_term,
            "gemini_message": gemini_result.get("message", ""),
            "pdb_ids": pdb_ids,
            "experimental_data": results,
            "csv_data": csv_data
        }
    except Exception as e:
        return {
            "status": "error",
            "error": f"Server error: {str(e)}",
            "gemini_message": ""
        }

# Run the server
if __name__ == "__main__":
    try:
        mcp.run(transport="stdio")
    except Exception as e:
        raise