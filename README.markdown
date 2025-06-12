# Chemical and Protein Data Analysis Suite

This repository contains three projects focused on querying and analyzing chemical and protein data: **Project_OMol25**, **Project_PubChem**, and **Project_ProteinDB**. Each project leverages specific databases and APIs to enable searches for chemical compounds or proteins, using modern web technologies and Python-based backends. Below is a summary of each project, including their purpose, technologies, and setup instructions.

## Table of Contents
- [Project Overview](#project-overview)
  - [Project_OMol25](#project_omol25)
  - [Project_PubChem](#project_pubchem)
  - [Project_ProteinDB](#project_proteindb)
- [General Prerequisites](#general-prerequisites)
- [Setup Instructions](#setup-instructions)

## Project Overview

### Project_OMol25
**Purpose**: Queries a Hugging Face database with 4 million chemical compound records stored in PostgreSQL. Searches are performed using hashed "chemical coordinates" for efficient retrieval.

**Technologies**:
- **Backend**: Python, MCP server, PostgreSQL, Gemini API
- **Frontend**: Vite/React app
- **Key Features**: Hashed search queries, PostgreSQL integration, Gemini API for enhanced processing

**Directory**: [Project_OMol25](./Project_OMol25)

### Project_PubChem
**Purpose**: Interfaces with the [PubChem database](https://pubchem.ncbi.nlm.nih.gov/) API to query approximately 120 million chemical compounds using chemical names or SMILES notation.

**Technologies**:
- **Backend**: Python, MCP server, Gemini API
- **Frontend**: Vite/React app
- **Key Features**: PubChem API integration, SMILES-based searches, Gemini API support

**Directory**: [Project_PubChem](./Project_PubChem)

### Project_ProteinDB
**Purpose**: Queries the [RCSB Protein Data Bank (PDB)](https://www.rcsb.org/) for millions of protein records using protein names or PDB IDs, with integration of licensed DNCS modeling software.

**Technologies**:
- **Backend**: Python (ProteinDB backend), TypeScript/Express (DNCS backend), Gemini API
- **Frontend**: Vite/React app
- **Key Features**: PDB API integration, DNCS modeling, dual-backend architecture

**Directory**: [Project_ProteinDB](./Project_ProteinDB)

## General Prerequisites
- **Python**: For server-side environments across all projects.
- **Node.js & npm**: For Vite/React frontends and TypeScript/Express backend (Project_ProteinDB).
- **Git**: To clone the repository.
- **Gemini API Key**: Obtain from [Google AI Studio](https://aistudio.google.com/) and store in a `.env` file in each project directory.
- **Database/API Access**:
  - PostgreSQL for Project_OMol25
  - PubChem API for Project_PubChem
  - RCSB PDB API for Project_ProteinDB
- **DNCS Software**: Licensed software for Project_ProteinDB.

## Setup Instructions

1. **Clone the Repository**
   ```bash
   git clone https://github.com/<your-username>/<your-repo-name>.git
   cd <your-repo-name>
   ```

2. **Project-Specific Setup**
   Each project has its own `README.markdown` with detailed setup instructions. Below is a high-level overview:

   **Project_OMol25**
   - Set up a PostgreSQL database with 4M records.
   - Create a `.env` file with the Gemini API key.
   - Activate a Python virtual environment, install server dependencies, and run `bridge.py` in `Project_OMol25/server`.
   - Install and run the Vite/React app in `Project_OMol25/client` using `npm i` and `npm run dev`.
   - Access at `http://localhost:5173`.
   - Details: [Project_OMol25 README](./Project_OMol25/README.markdown)

   **Project_PubChem**
   - Create a `.env` file with the Gemini API key.
   - Activate a Python virtual environment, install server dependencies, and run `bridge.py` in `Project_PubChem/server`.
   - Install and run the Vite/React app in `Project_PubChem/client` using `npm i` and `npm run dev`.
   - Access at `http://localhost:5173`.
   - Details: [Project_PubChem README](./Project_PubChem/README.markdown)

   **Project_ProteinDB**
   - Create a `.env` file with the Gemini API key.
   - Set up two backends:
     - Python backend: Activate a Python virtual environment, install dependencies, and run `bridge.py` in `Project_ProteinDB/server_backend_ProteinDB`.
     - TypeScript/Express backend: Install dependencies and run `npm run dev` in `Project_ProteinDB/server_backend_dncs`.
   - Install and run the Vite/React app in `Project_ProteinDB/client` using `npm i` and `npm run dev`.
   - Access at `http://localhost:5173`.
   - Details: [Project_ProteinDB README](./Project_ProteinDB/README.markdown)
