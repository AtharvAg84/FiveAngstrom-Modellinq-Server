# Project_OMol25

This project leverages a Hugging Face database containing 4 million queries stored in a PostgreSQL database. The data is hashed, and searches are performed using the "chemical coordinates" column. The project uses an MCP server and integrates the Gemini API for enhanced functionality. This README provides instructions to set up and run the project.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the Project](#running-the-project)
- [Project Structure](#project-structure)
- [Notes](#notes)
- [License](#license)

## Prerequisites
- **PostgreSQL Database**: A database with 4 million chemical compound records.
- **Python**: For the server-side environment.
- **Node.js & npm**: For the client-side Vite/React app.
- **Git**: To clone the repository.
- **Gemini API Key**: Obtain from [Google AI Studio](https://aistudio.google.com/) and store in a `.env` file.

## Installation

1. **Clone the Repository**
   ```bash
   git clone https://github.com/<your-username>/Project_OMol25.git
   cd Project_OMol25
   ```

2. **Set Up the Gemini API Key**
   - Create a `.env` file in the `Project_OMol25` directory.
   - Add your Gemini API key:
     ```
     GEMINI_API_KEY=your_api_key_here
     ```

3. **Database Setup**
   - Ensure a PostgreSQL database is set up with 4 million records.
   - The "chemical coordinates" column in the database is hashed and used for search queries.
   - Refer to the database configuration details at the end of `server.py` for setup.

4. **Server Setup**
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

5. **Client Setup**
   - Open a new terminal and navigate to the client directory:
     ```bash
     cd Project_OMol25/client
     ```
   - Install dependencies:
     ```bash
     npm i
     ```

## Running the Project

1. **Start the Server**
   - In the server directory (`Project_OMol25/server`), run:
     ```bash
     python bridge.py
     ```

2. **Start the Client**
   - In the client directory (`Project_OMol25/client`), run:
     ```bash
     npm run dev
     ```

3. **Access the Application**
   - Open your browser and navigate to `http://localhost:5173`.

## Project Structure
```
Project_OMol25/
├── server/
│   ├── bridge.py              # Main server script
│   ├── client.py              # Client-side server communication
│   ├── fetch_molecule.py      # Script to fetch molecule data
│   ├── hashing_dataset.py     # Script for hashing dataset
│   ├── postgres_dataset.py    # Script for PostgreSQL dataset management
│   ├── server.py              # MCP server configuration
│   └── requirements.txt       # Server dependencies
├── client/
│   ├── package.json           # Client dependencies for Vite/React app
│   └── src/                   # Frontend source code
├── .env                       # Environment variables (e.g., Gemini API key)
└── README.markdown            # Project documentation
```

## Notes
- Ensure the PostgreSQL database is properly configured and running before starting the server.
- The search functionality relies on the hashed "chemical coordinates" column.
- The Gemini API key must be correctly set in the `.env` file.
- The client is a Vite/React app, providing a fast and modern frontend experience.
- Check terminal logs for debugging issues in the server or client.

## License
This repository and its projects are licensed under the Five Angstrom LLP License.