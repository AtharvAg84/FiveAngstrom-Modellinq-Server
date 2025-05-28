import logging
import json
import asyncio
import os
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import google.generativeai as genai
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
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if not GEMINI_API_KEY:
    logger.error("GEMINI_API_KEY not found in .env file")
    raise ValueError("GEMINI_API_KEY not found in .env file")
genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel("gemini-1.5-flash")

# FastAPI setup
app = FastAPI(title="PDB Search Bridge")

# Enable CORS for React app
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class QueryRequest(BaseModel):
    query: str
    limit: int = 10

@app.post("/search_pdb")
async def search_pdb(request: QueryRequest):
    query = request.query
    limit = request.limit
    logger.info(f"Received query: {query}, limit: {limit}")
    try:
        # Validate inputs
        if not query or len(query.strip()) == 0:
            logger.warning("Empty or invalid query")
            return {"status": "error", "error": "Query cannot be empty"}
        if len(query) > 500:
            logger.warning(f"Query too long: {query[:50]}...")
            return {"status": "error", "error": "Query exceeds maximum length of 500 characters"}
        if limit < 1 or limit > 100:
            logger.warning(f"Invalid limit: {limit}")
            return {"status": "error", "error": "Limit must be between 1 and 100"}

        # Use Gemini to preprocess query
        logger.debug("Sending query to Gemini for preprocessing")
        prompt = f"""
        You are a chemistry expert. If the input is a chemical name or description (including misspellings like 'aspirine' for 'aspirin'),
        return the corrected or standardized chemical name. If the input is valid, return it as-is.
        Input: {query}
        Output format: {{ "corrected_query": "standardized name", "message": "explanation" }}
        """
        gemini_response = model.generate_content(prompt)
        if not hasattr(gemini_response, 'text') or not gemini_response.text:
            logger.error("Invalid Gemini response: No text content")
            return {"status": "error", "error": "Failed to get response from Gemini"}

        gemini_data = gemini_response.text.strip()
        if gemini_data.startswith("```json"):
            gemini_data = gemini_data[7:].strip()
        if gemini_data.endswith("```"):
            gemini_data = gemini_data[:-3].strip()
        try:
            gemini_result = json.loads(gemini_data)
            corrected_query = gemini_result.get("corrected_query", query)
            gemini_message = gemini_result.get("message", "")
            logger.info(f"Gemini corrected query to: {corrected_query}")
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse Gemini response: {str(e)}")
            return {"status": "error", "error": f"Failed to parse Gemini response: {str(e)}"}

        # Search for PDB IDs with corrected query
        logger.info(f"Searching PDB IDs for query: {corrected_query}")
        pdb_ids = search_pdb_ids(corrected_query)
        logger.debug(f"Found {len(pdb_ids)} PDB IDs")

        if not pdb_ids:
            logger.info("No PDB IDs found for query")
            return {
                "status": "success",
                "message": f"No PDB IDs found for '{corrected_query}'. {gemini_message}",
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
            "message": f"Found {len(pdb_ids)} PDB IDs, processed {len(results)}. {gemini_message}",
            "csv_data": csv_data,
            "corrected_query": corrected_query
        }
    except Exception as e:
        logger.error(f"Error processing query '{query}': {str(e)}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Server error: {str(e)}")

if __name__ == "__main__":
    import uvicorn
    logger.info("Starting FastAPI bridge server")
    uvicorn.run(app, host="0.0.0.0", port=8000)