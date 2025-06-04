import { useState, useEffect, useRef } from 'react';
import * as NGL from 'ngl';
import './App.css';

function NGLViewer({ pdbContent, viewerFeatures, isPdbLoading, setViewerFeatures }) {
  const stageRef = useRef(null);
  const viewerRef = useRef(null);

  useEffect(() => {
    if (!viewerRef.current) return;

    // Initialize NGL Stage
    if (!stageRef.current) {
      stageRef.current = new NGL.Stage(viewerRef.current, {
        backgroundColor: 'black',
        cameraType: 'perspective',
      });

      // Handle window resize
      const handleResize = () => {
        if (stageRef.current) {
          stageRef.current.handleResize();
        }
      };
      window.addEventListener('resize', handleResize);
      return () => window.removeEventListener('resize', handleResize);
    }

    // Load PDB content
    const loadStructure = async () => {
      if (isPdbLoading || !pdbContent) return;

      try {
        // Clear existing components
        stageRef.current.removeAllComponents();

        // Load PDB as a string
        const blob = new Blob([pdbContent], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const structure = await stageRef.current.loadFile(url, { ext: 'pdb' });
        URL.revokeObjectURL(url);

        // Apply representations based on viewerFeatures
        if (viewerFeatures.backbone) {
          structure.addRepresentation('backbone', { sele: ':A', color: 'blue', radius: 0.5 });
        }
        if (viewerFeatures.cartoon) {
          structure.addRepresentation('cartoon', { sele: ':A', color: 'green', radius: 0.3 });
        }
        if (viewerFeatures.line) {
          structure.addRepresentation('line', { sele: 'all', color: 'grey', opacity: 0.7 });
        }
        if (viewerFeatures.ballAndStick) {
          structure.addRepresentation('ball+stick', {
            sele: 'not protein and not water',
            color: 'red',
            radiusScale: 0.5,
          });
        }
        if (viewerFeatures.label) {
          structure.addRepresentation('label', {
            sele: 'sidechain',
            labelType: 'resname',
            color: 'white',
            scale: 0.8,
            showBackground: true,
          });
        }

        // Auto-view the structure
        stageRef.current.autoView(1000);
      } catch (error) {
        console.error('Error loading PDB in NGL:', error);
      }
    };

    loadStructure();

    return () => {
      // Cleanup components (but keep stage for reuse)
      if (stageRef.current) {
        stageRef.current.removeAllComponents();
      }
    };
  }, [pdbContent, isPdbLoading, viewerFeatures]);

  return (
    <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 w-full">
      <h3 className="text-base font-semibold mb-2 text-cyan-300">3D Structure Viewer</h3>
      {isPdbLoading ? (
        <p className="text-gray-200 text-sm">Loading PDB content...</p>
      ) : pdbContent ? (
        <div>
          <div
            ref={viewerRef}
            className="w-full"
            style={{ height: '400px', position: 'relative' }}
          />
          <div className="mt-3 flex flex-wrap gap-2">
            {Object.keys(viewerFeatures).map((feature) => (
              <button
                key={feature}
                onClick={() =>
                  setViewerFeatures((prev) => ({
                    ...prev,
                    [feature]: !prev[feature],
                  }))
                }
                className={`py-1 px-3 text-sm rounded ${
                  viewerFeatures[feature]
                    ? 'bg-cyan-600 text-gray-200'
                    : 'bg-gray-700/50 text-gray-300'
                } hover:bg-cyan-700`}
              >
                {feature.charAt(0).toUpperCase() + feature.slice(1)}
              </button>
            ))}
          </div>
        </div>
      ) : (
        <p className="text-gray-300 text-sm">Enter a query to view 3D structure</p>
      )}
    </div>
  );
}

function App() {
  const [query, setQuery] = useState('');
  const [response, setResponse] = useState(null);
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [viewerFeatures, setViewerFeatures] = useState({
    backbone: false,
    cartoon: true,
    line: false,
    ballAndStick: false,
    label: false,
  });

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (query.length > 100) {
      setError('Query must be 100 characters or less.');
      return;
    }
    setLoading(true);
    setError('');
    setResponse(null);
    setViewerFeatures({
      backbone: false,
      cartoon: true,
      line: false,
      ballAndStick: false,
      label: false,
    });

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
      // Validate pdb_content
      if (data.pdb_content && !data.pdb_content.includes('ATOM') && !data.pdb_content.includes('HETATM')) {
        setError('Invalid PDB data: No ATOM or HETATM records found.');
        setResponse({ ...data, pdb_content: null });
      } else {
        setResponse(data);
      }
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
                        disabled={loading}
                      >
                        Download XYZ
                      </button>
                      <button
                        onClick={() => downloadFile(response.pdb_content, response.pdb_filename)}
                        className="download-button"
                        disabled={loading}
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
        <NGLViewer
          pdbContent={response?.pdb_content}
          viewerFeatures={viewerFeatures}
          isPdbLoading={loading}
          setViewerFeatures={setViewerFeatures}
        />
      </div>
    </div>
  );
}

export default App;