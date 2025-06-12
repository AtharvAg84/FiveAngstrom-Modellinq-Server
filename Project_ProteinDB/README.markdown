# Project_ProteinDB

This project leverages the [RCSB Protein Data Bank (PDB)](https://www.rcsb.org/) to query millions of protein records using Python libraries for data fetching and searching. Queries can be performed using either a protein name or a PDB ID. The project integrates licensed DNCS modeling software and uses an MCP server with the Gemini API for enhanced functionality. This README guides you through setting up and running the project.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the Project](#running-the-project)
- [Project Structure](#project-structure)
- [Notes](#notes)
- [Contributing](#contributing)

## Prerequisites
- **Python**: For the server-side environment.
- **Node.js & npm**: For the client-side Vite/React app and TypeScript-based Express server.
- **Git**: To clone the repository.
- **Gemini API Key**: Obtain from [Google AI Studio](https://aistudio.google.com/) and store in a `.env` file.
- **DNCS Software**: Licensed modeling software (ensure proper licensing).
- **RCSB PDB Access**: Ensure access to the RCSB PDB API for querying protein data.

## Installation

1. **Clone the Repository**
   ```bash
   git clone https://github.com/<your-username>/Project_ProteinDB.git
   cd Project_ProteinDB
   ```

2. **Set Up the Gemini API Key**
   - Create a `.env` file in the `Project_ProteinDB` directory.
   - Add your Gemini API key:
     ```
     GEMINI_API_KEY=your_api_key_here
     ```

3. **Server Setup (ProteinDB Backend)**
   - Create and activate a Python virtual environment:
     ```bash
     # Linux/MacOS
     python3 -m venv venv
     source venv/bin/activate

     # Windows
     python -m venv venv
     .\venv\Scripts\activate
     ```
   - Navigate to the `server_backend_ProteinDB` directory and install dependencies:
     ```bash
     cd server_backend_ProteinDB
     pip install -r requirements.txt
     ```

4. **Server Setup (DNCS Backend)**
   - Navigate to the `server_backend_dncs` directory:
     ```bash
     cd server_backend_dncs
     ```
   - Install dependencies:
     ```bash
     npm i
     ```

5. **Client Setup**
   - Open a new terminal and navigate to the client directory:
     ```bash
     cd Project_ProteinDB/client
     ```
   - Install dependencies:
     ```bash
     npm i
     ```

## Running the Project

1. **Start the ProteinDB Backend Server**
   - In the `server_backend_ProteinDB` directory, run:
     ```bash
     python bridge.py
     ```
   - The MCP server is configured in `pdb_search_server.py`.

2. **Start the DNCS Backend Server**
   - In the `server_backend_dncs` directory, run:
     ```bash
     npm run dev
     ```

3. **Start the Client**
   - In the client directory (`Project_ProteinDB/client`), run:
     ```bash
     npm run dev
     ```

4. **Access the Application**
   - Open your browser and navigate to `http://localhost:5173`.
   - Use protein names or PDB IDs to query protein information via the RCSB PDB API.

## Project Structure
```
Project_ProteinDB/
├── client/
│   ├── package.json           # Client dependencies for Vite/React app
│   └── src/                   # Frontend source code
├── server_backend_ProteinDB/
│   ├── bridge.py              # Main server script for API interactions
│   ├── fetch_pdb.py           # Script to fetch PDB data
│   ├── pdb_search_client.py   # Client-side server communication
│   ├── pdb_search_server.py   # MCP server configuration
│   ├── search_data_api.py     # API search functionality
│   └── requirements.txt       # Server dependencies
├── server_backend_dncs/
│   ├── package.json           # Dependencies for TypeScript Express server
│   └── src/                   # TypeScript source code
├── package-lock.json          # Dependency lock file
├── sample_server.py           # Sample server script
├── test.ts                    # Test script
├── .env                       # Environment variables (e.g., Gemini API key)
└── README.markdown            # Project documentation
```

## Notes
- Ensure valid protein names or PDB IDs for successful queries.
- The Gemini API key must be correctly set in the `.env` file.
- The `server_backend_ProteinDB` uses Python for RCSB PDB API interactions, while `server_backend_dncs` is an Express framework in TypeScript for DNCS integration.
- The client is a Vite/React app, ensuring a fast and modern frontend experience.
- Verify that the DNCS software is properly licensed and configured.
- Check terminal logs for debugging issues in the servers or client.

## Contributing
Contributions are welcome! Please fork the repository and submit a pull request with your changes.
