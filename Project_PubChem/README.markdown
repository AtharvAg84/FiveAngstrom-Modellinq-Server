# Project_PubChem

This project leverages the [PubChem database](https://pubchem.ncbi.nlm.nih.gov/) API, providing access to approximately 120 million chemical compound queries. Searches can be performed using either chemical names or SMILES notation. The backend uses an MCP server and integrates the Gemini API for enhanced functionality. This README guides you through setting up and running the project.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the Project](#running-the-project)
- [Project Structure](#project-structure)
- [Notes](#notes)

## Prerequisites
- **Python**: For the server-side environment.
- **Node.js & npm**: For the client-side React Vite app.
- **Git**: To clone the repository.
- **Gemini API Key**: Obtain from [Google AI Studio](https://aistudio.google.com/) and store in a `.env` file.

## Installation

1. **Clone the Repository**
   ```bash
   git clone https://github.com/<your-username>/Project_PubChem.git
   cd Project_PubChem
   ```

2. **Set Up the Gemini API Key**
   - Create a `.env` file in the `Project_PubChem` directory.
   - Add your Gemini API key:
     ```
     GEMINI_API_KEY=your_api_key_here
     ```

3. **Server Setup**
   - Create and activate a Python virtual environment:
     ```bash
     # Linux/MacOS
     python3 -m venv venv
     source venv/bin/activate

     # Windows
     python -m venv venv
     .\venv\Scripts\activate
     ```
   - Navigate to the server directory and install dependencies:
     ```bash
     cd server
     pip install -r requirements.txt
     ```

4. **Client Setup**
   - Open a new terminal and navigate to the client directory:
     ```bash
     cd Project_PubChem/client
     ```
   - Install dependencies:
     ```bash
     npm i
     ```

## Running the Project

1. **Start the Server**
   - In the server directory (`Project_PubChem/server`), run:
     ```bash
     python bridge.py
     ```
   - The MCP server is configured in `server.py` and is necessary for API interactions.

2. **Start the Client**
   - In the client directory (`Project_PubChem/client`), run:
     ```bash
     npm run dev
     ```

3. **Access the Application**
   - Open your browser and navigate to `http://localhost:5173`.
   - Use chemical names or SMILES notation to query molecule information via the PubChem API.

## Project Structure
```
Project_PubChem/
├── server/
│   ├── bridge.py              # Main server script for API interactions
│   ├── client.py              # Client-side server communication
│   ├── fetch_coordinates.py   # Script to fetch chemical coordinates
│   ├── server.py              # MCP server configuration
│   └── requirements.txt       # Server dependencies
├── client/
│   ├── package.json           # Client dependencies for React Vite app
│   └── src/                   # Frontend source code
├── .env                       # Environment variables (e.g., Gemini API key)
└── README.markdown            # Project documentation
```

## Notes
- Ensure valid chemical names or SMILES notation for successful queries.
- The Gemini API key must be correctly set in the `.env` file.
- The MCP server (`server.py`) handles PubChem API requests and Gemini API integration.
- Check terminal logs for debugging issues in the server or client.
- The client is a React Vite app, ensuring a fast and modern frontend experience.
