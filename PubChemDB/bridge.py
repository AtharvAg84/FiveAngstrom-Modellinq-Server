#bridge.py

import logging
import json
import asyncio
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from fastmcp.client import Client

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# FastAPI setup
app = FastAPI(title="Chemical Coordinates Bridge")

# Enable CORS for React app (e.g., running on localhost:5173)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class QueryRequest(BaseModel):
    query: str

@app.post("/get_coordinates")
async def get_coordinates(request: QueryRequest):
    """
    Receive HTTP POST request from React app and call FastMCP server via stdio.
    """
    query = request.query
    logger.info(f"Received query: {query}")
    try:
        # Initialize FastMCP client to communicate with server.py
        client = Client("server.py")
        logger.info("Connected to FastMCP server via stdio")

        # Call the get_coordinates tool
        async with client:
            logger.debug(f"Calling tool get_coordinates with query: {query}")
            response = await client.call_tool("get_coordinates", {"query": query})
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

if __name__ == "__main__":
    import uvicorn
    logger.info("Starting FastAPI bridge server")
    uvicorn.run(app, host="0.0.0.0", port=8000)
