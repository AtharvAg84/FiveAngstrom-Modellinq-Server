import asyncio
import json
import time
import logging
import os

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Server details
SERVER_SCRIPT = "server.py"
TOOL_NAME = "get_database_info"

async def call_mcp_server(query: str, retries=3) -> dict:
    """
    Send a query to the FastMCP server using the stdio transport.

    Args:
        query (str): Query about the database (e.g., "list tables", "data in table users where name is Alice").
        retries (int): Number of retry attempts for connection failures.

    Returns:
        dict: Server response containing database information or error.
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
                # Call the get_database_info tool with query as a dict
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
    print("Welcome to the FastMCP Database Information Client!")
    print("Enter a query about the PostgreSQL database 'something'.")
    print("Examples: 'list tables', 'schema of table users', 'data in table users where name = ''Alice''', 'stats', 'search user'")
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
            if "tables" in result:
                print("Tables in database:")
                for table in result["tables"]:
                    print(f"- {table}")
            elif "schema" in result:
                print(f"Schema of table:")
                for column in result["schema"]:
                    print(f"Column: {column['column_name']}, Type: {column['data_type']}")
            elif "row_count" in result:
                print(f"Row count: {result['row_count']}")
            elif "data" in result:
                print("Table data:")
                for row in result["data"]:
                    print(row)
            elif "column_info" in result:
                info = result["column_info"]
                print(f"Column Info: {info['column_name']}")
                print(f"Table: {info['table_name']}")
                print(f"Data Type: {info['data_type']}")
                print(f"Nullable: {info['is_nullable']}")
            elif "stats" in result:
                stats = result["stats"]
                print("Database atistics:")
                print(f"Total Tables: {stats['total_tables']}")
                print(f"Total Rows: {stats['total_rows']}")
            elif "search_results" in result:
                print("Search Results:")
                for item in result["search_results"]:
                    if item["type"] == "table":
                        print(f"Table: {item['table_name']}")
                    else:
                        print(f"Column: {item['column_name']} in Table: {item['table_name']}")
            print(f"Gemini Message: {result.get('gemini_message')}")
        else:
            print(f"Error: {result.get('error')}")
            if result.get("gemini_message"):
                print(f"Gemini Message: {result.get('gemini_message')}")

if __name__ == "__main__":
    asyncio.run(main())