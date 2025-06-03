import { useState } from 'react';
import './App.css';

function App() {
  const [query, setQuery] = useState('');
  const [response, setResponse] = useState(null);
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (query.length > 100) {
      setError('Query must be 100 characters or less.');
      return;
    }
    setLoading(true);
    setError('');
    setResponse(null);

    try {
      const res = await fetch('http://localhost:8000/search', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query }),
      });

      if (!res.ok) {
        throw new Error(`HTTP error! Status: ${res.status}`);
      }

      const data = await res.json();
      setResponse(data);
    } catch (err) {
      setError(`Failed to fetch response: ${err.message}`);
    } finally {
      setLoading(false);
    }
  };

  const downloadFile = (content, filename) => {
    const blob = new Blob([content], { type: 'text/plain' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    window.URL.revokeObjectURL(url);
  };

  return (
    <div className="app-container">
      <div className="left-column">
        <h1>Molecule Search</h1>
        <p>Enter a chemical composition (e.g., 'C9H8O4'), name (e.g., 'aspirin'), or chat query (e.g., 'hello').</p>
        <form onSubmit={handleSubmit} className="form">
          <input
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="Enter query"
            className="input"
            disabled={loading}
          />
          <button type="submit" className="button" disabled={loading}>
            {loading ? 'Searching...' : 'Search'}
          </button>
        </form>

        {error && <div className="error">{error}</div>}

        {response && (
          <div className="response">
            {response.status === 'success' ? (
              <>
                {response.chat_message ? (
                  <div>
                    <h2>Chat Response</h2>
                    <p>{response.chat_message}</p>
                  </div>
                ) : (
                  <div>
                    <h2>Molecule Data</h2>
                    <p><strong>Composition:</strong> {response.gemini_composition}</p>
                    <p><strong>Chemical Formula:</strong> {response.molecule_data?.chemical_formula}</p>
                    <p><strong>Number of Atoms:</strong> {response.molecule_data?.num_atoms}</p>
                    <p><strong>Message:</strong> {response.gemini_message}</p>
                    <div className="download-buttons">
                      <button
                        onClick={() => downloadFile(response.xyz_content, response.xyz_filename)}
                        className="download-button"
                      >
                        Download XYZ
                      </button>
                      <button
                        onClick={() => downloadFile(response.pdb_content, response.pdb_filename)}
                        className="download-button"
                      >
                        Download PDB
                      </button>
                    </div>
                    <div className="file-content">
                      <h3>XYZ File Content</h3>
                      <pre className="file-content-pre">{response.xyz_content}</pre>
                    </div>
                    <div className="file-content">
                      <h3>PDB File Content</h3>
                      <pre className="file-content-pre">{response.pdb_content}</pre>
                    </div>
                  </div>
                )}
              </>
            ) : (
              <div>
                <h2 className="error-title">Error</h2>
                <p>{response.error}</p>
                {response.gemini_message && (
                  <p><strong>Gemini Message:</strong> {response.gemini_message}</p>
                )}
              </div>
            )}
          </div>
        )}
      </div>
      <div className="right-column">
        <h2>Placeholder</h2>
        <p>This area is reserved for additional content (e.g., 3D visualization).</p>
      </div>
    </div>
  );
}

export default App;
