import logging
import json
import asyncio
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from fastmcp.client import Client
from fetch_pdb import PDBFetcher

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

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

# Initialize PDBFetcher
pdb_fetcher = PDBFetcher()

@app.post("/search")
async def search_pdb(request: QueryRequest):
    """
    Receive HTTP POST request from React app and call FastMCP server via stdio.
    """
    query = request.query
    limit = request.limit
    logger.info(f"Received query: {query}, limit: {limit}")
    try:
        # Initialize FastMCP client to communicate with pdb_search_server.py
        client = Client("pdb_search_server.py")
        logger.info("Connected to FastMCP server via stdio")

        # Call the search_pdb_data tool
        async with client:
            logger.debug(f"Calling tool search_pdb_data with query: {query}, limit: {limit}")
            response = await client.call_tool("search_pdb_data", {"query": query, "limit": limit})
            logger.info("Received response from server")

            # Handle response
            if isinstance(response, list):
                if len(response) == 0:
                    raise ValueError("Empty response list from server")
                response = response[0]
                logger.info("Extracted first element from response list")

            # Handle TextContent-like object
            if hasattr(response, 'text'):
                try:
                    response_text = response.text.strip()
                    if response_text.startswith("```json"):
                        response_text = response_text[7:].strip()
                    if response_text.endswith("```"):
                        response_text = response_text[:-3].strip()
                    response = json.loads(response_text)
                    logger.info("Parsed TextContent text to dictionary")
                except json.JSONDecodeError as e:
                    logger.error(f"Failed to parse TextContent text: {str(e)}")
                    raise HTTPException(status_code=500, detail=f"Failed to parse response: {str(e)}")

            # Ensure response is a dictionary
            if not isinstance(response, dict):
                logger.error(f"Unexpected response type: {type(response)}")
                raise HTTPException(status_code=500, detail=f"Unexpected response type: {type(response)}")

            return response
    except Exception as e:
        logger.error(f"Error processing query '{query}': {str(e)}")
        raise HTTPException(status_code=500, detail=f"Server error: {str(e)}")

@app.get("/fetch-pdb/{pdb_id}")
async def fetch_pdb(pdb_id: str):
    """
    Fetch PDB file content for a given PDB ID using PDBFetcher class.
    
    Args:
        pdb_id (str): The PDB ID (e.g., '1A3N').
    
    Returns:
        dict: Contains the PDB ID and file content, or raises an error if failed.
    """
    logger.info(f"Received request to fetch PDB file for ID: {pdb_id}")
    try:
        pdb_content = pdb_fetcher.fetch_pdb_file(pdb_id.upper())
        if not pdb_content:
            logger.error(f"PDB file not found for ID: {pdb_id}")
            raise HTTPException(status_code=404, detail=f"PDB file not found for ID: {pdb_id}")
        logger.info(f"Successfully fetched PDB file for ID: {pdb_id}")
        return {"pdb_id": pdb_id, "pdb_content": pdb_content}
    except Exception as e:
        logger.error(f"Error fetching PDB file for {pdb_id}: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Error fetching PDB file: {str(e)}")

if __name__ == "__main__":
    import uvicorn
    logger.info("Starting FastAPI bridge server")
    uvicorn.run(app, host="0.0.0.0", port=8000)