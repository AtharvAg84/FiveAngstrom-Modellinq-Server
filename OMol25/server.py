import logging
import json
import os
import re
from fastmcp import FastMCP
import google.generativeai as genai
from fetch_molecule import MoleculeFetcher
from dotenv import load_dotenv
from typing import Dict, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MoleculeServer:
    """FastMCP server for searching molecules by composition or handling chat queries."""
    
    def __init__(self, db_params: Dict[str, str], gemini_api_key: str):
        """
        Initialize server with database and Gemini API settings.
        
        Args:
            db_params (Dict[str, str]): Database connection parameters (dbname, user, password, host, port).
            gemini_api_key (str): Gemini API key for validation.
        
        Raises:
            ValueError: If Gemini API key is missing.
        """
        self.db_params = db_params
        self.fetcher = MoleculeFetcher(db_params)
        self.mcp = FastMCP("Molecule Search Server")
        logger.info("Initialized FastMCP server")

        # Configure Gemini
        genai.configure(api_key=gemini_api_key)
        self.model = genai.GenerativeModel("gemini-1.5-flash")
        logger.info("Initialized Gemini model")

    def is_chemical_composition(self, composition: str) -> bool:
        """
        Check if the input resembles a chemical composition.
        
        Args:
            composition (str): Input string to validate.
            
        Returns:
            bool: True if the input matches a chemical composition pattern (e.g., 'C9H8O4').
        """
        pattern = r'^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$'
        return bool(re.match(pattern, composition))

    def handle_chat_query(self, query: str) -> Dict:
        """
        Process non-chemical inputs with a conversational Gemini query.
        
        Args:
            query (str): User input query.
            
        Returns:
            Dict: Response containing a chat message or error.
        """
        logger.info(f"Processing chat query: {query}")
        prompt = """
        You are a friendly chemistry assistant. Respond conversationally to the user's query.
        If the query is unrelated to chemistry, provide a helpful and appropriate response.
        Input: {query}
        Output format: {{"message": "response text"}}
        """
        try:
            response = self.model.generate_content(prompt.format(query=query))
            if not hasattr(response, 'text') or not response.text:
                logger.error("No text in Gemini chat response")
                return {"status": "error", "error": "No response from Gemini"}

            text = response.text.strip()
            logger.debug(f"Raw Gemini chat response: {text}")

            if text.startswith("```json"):
                text = text[7:].strip()
            if text.endswith("```"):
                text = text[:-3].strip()

            try:
                result = json.loads(text)
                return {"status": "success", "chat_message": result.get("message", text)}
            except json.JSONDecodeError:
                return {"status": "success", "chat_message": text}
        except Exception as e:
            logger.error(f"Gemini chat error: {str(e)}")
            return {"status": "error", "error": f"Gemini error: {str(e)}"}

    def validate_composition(self, composition: str) -> Dict:
        """
        Validate or correct a chemical composition or name using Gemini API.
        
        Args:
            composition (str): Input composition or chemical name.
            
        Returns:
            Dict: Validation result with composition, is_valid, and message.
        """
        logger.info(f"Validating input: {composition}")
        if not composition or not any(c.isalnum() for c in composition):
            return {
                "is_valid": False,
                "message": "Invalid input: Please provide a chemical composition (e.g., 'C9H8O4') or name (e.g., 'aspirin')."
            }

        if self.is_chemical_composition(composition):
            return {
                "composition": composition,
                "is_valid": True,
                "message": "Valid chemical composition"
            }

        prompt = """
        You are a chemistry expert. Validate the input as a chemical composition or name.
        If it's a valid composition (e.g., 'C9H8O4'), return it as is.
        If it's a chemical name or misspelled (e.g., 'aspirine'), correct it to a composition (e.g., 'C9H8O4').
        If unrecognizable, return is_valid as false.
        Input: {composition}
        Output format: {{"composition": "valid composition", "is_valid": true/false, "message": "explanation"}}
        """
        try:
            response = self.model.generate_content(prompt.format(composition=composition))
            if not hasattr(response, 'text') or not response.text:
                logger.error("No text in Gemini response")
                return {
                    "is_valid": False,
                    "message": "Gemini API returned no content."
                }

            text = response.text.strip()
            logger.debug(f"Raw Gemini response: {text}")

            if not text:
                return {
                    "is_valid": False,
                    "message": "Gemini API returned an empty response."
                }

            if text.startswith("```json"):
                text = text[7:].strip()
            if text.endswith("```"):
                text = text[:-3].strip()

            result = json.loads(text)
            if not isinstance(result, dict) or "is_valid" not in result:
                logger.error("Invalid Gemini response format")
                return {
                    "is_valid": False,
                    "message": "Gemini API returned invalid data."
                }
            return result
        except json.JSONDecodeError as e:
            logger.error(f"Gemini JSON parse error: {str(e)}")
            return {
                "is_valid": False,
                "message": f"Gemini API returned malformed data: {str(e)}."
            }
        except Exception as e:
            logger.error(f"Gemini API error: {str(e)}")
            return {
                "is_valid": False,
                "message": f"Gemini API error: {str(e)}."
            }

    def search_molecule(self, composition: str) -> Dict:
        """
        Search for a molecule by composition or process a chat query.
        
        Args:
            composition (str): Input string (composition, name, or chat query).
            
        Returns:
            Dict: Molecule data with XYZ content or chat response.
        """
        logger.info(f"Processing input: {composition[:50]}{'...' if len(composition) > 50 else ''}")
        if len(composition) > 100:
            return {
                "status": "error",
                "error": "Input exceeds 100 characters",
                "gemini_message": ""
            }

        # Handle non-chemical inputs
        if not self.is_chemical_composition(composition) and not re.search(r'[A-Z][a-z]?\d*', composition):
            return self.handle_chat_query(composition)

        # Handle chemical compositions or names
        gemini_result = self.validate_composition(composition)
        validated_composition = gemini_result.get("composition")
        if not validated_composition or not gemini_result.get("is_valid"):
            return {
                "status": "error",
                "error": gemini_result.get("message", "Invalid composition"),
                "gemini_message": gemini_result.get("message", "")
            }

        molecule_data = self.fetcher.query_molecule(validated_composition)
        if not molecule_data:
            return {
                "status": "error",
                "error": f"No molecule found for composition: {validated_composition}",
                "gemini_composition": validated_composition,
                "gemini_message": gemini_result.get("message", "")
            }

        try:
            xyz_content = self.fetcher.to_xyz(molecule_data)
            if not xyz_content:
                logger.error("Failed to generate XYZ content")
                return {
                    "status": "error",
                    "error": "Failed to generate XYZ content",
                    "gemini_composition": validated_composition,
                    "gemini_message": gemini_result.get("message", "")
                }
        except Exception as e:
            logger.error(f"XYZ generation error: {str(e)}")
            return {
                "status": "error",
                "error": f"Failed to generate XYZ content: {str(e)}",
                "gemini_composition": validated_composition,
                "gemini_message": gemini_result.get("message", "")
            }

        return {
            "status": "success",
            "molecule_data": molecule_data,
            "xyz_content": xyz_content,
            "xyz_filename": f"{validated_composition}_molecule.xyz",  # Use validated composition
            "gemini_composition": validated_composition,
            "gemini_message": gemini_result.get("message", "")
        }

    def run(self) -> None:
        """
        Run the FastMCP server with the search_molecule tool.
        
        Raises:
            RuntimeError: If the server fails to start.
        """
        self.mcp.tool()(self.search_molecule)
        try:
            logger.info("Starting FastMCP server")
            self.mcp.run(transport="stdio")
            logger.info("Server started")
        except Exception as e:
            logger.error(f"Failed to start server: {str(e)}")
            raise RuntimeError(f"Server startup failed: {str(e)}")

if __name__ == "__main__":
    path="../.env"
    load_dotenv(path)
    GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
    if not GEMINI_API_KEY:
        logger.error("GEMINI_API_KEY not found in .env")
        raise ValueError("GEMINI_API_KEY not found")

    db_params = {
        "dbname": "omol4m",
        "user": "atharvag",
        "password": "atharv8484#",
        "host": "localhost",
        "port": "5432"
    }
    server = MoleculeServer(db_params, GEMINI_API_KEY)
    server.run()