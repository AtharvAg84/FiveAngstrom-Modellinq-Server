import { useState, useEffect, useRef } from 'react';
import * as NGL from 'ngl';

function App() {
  const [query, setQuery] = useState('');
  const [result, setResult] = useState(null);
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [pdbContent2D, setPdbContent2D] = useState('');
  const [pdbContent3D, setPdbContent3D] = useState('');
  const [viewMode, setViewMode] = useState('3D');
  const stageRef = useRef(null);
  const viewerRef = useRef(null);
  const [viewerError, setViewerError] = useState('');

  // Validate PDB content
  const isValidPdb = (content) => {
    if (!content || typeof content !== 'string') {
      console.log('PDB validation failed: Content is empty or not a string');
      return false;
    }
    const lines = content.split('\n').filter(line => line.trim());
    const hasAtoms = lines.some(line => line.startsWith('ATOM') || line.startsWith('HETATM'));
    console.log(`Validating PDB: Has ATOM/HETATM lines: ${hasAtoms}, Total lines: ${lines.length}`);
    return hasAtoms;
  };

  // Check WebGL availability
  useEffect(() => {
    const canvas = document.createElement('canvas');
    const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
    if (!gl) {
      setViewerError('WebGL is not supported in your browser. Please enable WebGL or use a compatible browser.');
      console.error('WebGL not supported');
    }
  }, []);

  // Initialize NGL Stage
  useEffect(() => {
    const initNGL = () => {
      if (viewerRef.current && !stageRef.current) {
        console.log('Initializing NGL Stage');
        try {
          stageRef.current = new NGL.Stage(viewerRef.current, { backgroundColor: 'white' });
          const rect = viewerRef.current.getBoundingClientRect();
          console.log('Viewer container dimensions:', { width: rect.width, height: rect.height });
          if (rect.width === 0 || rect.height === 0) {
            setViewerError('Viewer container has zero dimensions. Check CSS styles or parent container size.');
            console.error('Zero dimensions detected for viewer container');
          }
          stageRef.current.handleResize();
          window.addEventListener('resize', debounce(handleResize, 100));
        } catch (err) {
          console.error('NGL Stage initialization error:', err);
          setViewerError('Failed to initialize NGL Viewer: ' + err.message);
        }
      }
    };

    // Delay initialization to ensure container is rendered
    const timer = setTimeout(initNGL, 100);
    return () => {
      clearTimeout(timer);
      window.removeEventListener('resize', debounce(handleResize, 100));
      if (stageRef.current) {
        stageRef.current.dispose();
        stageRef.current = null;
      }
    };
  }, []);

  // Debounce function for resize
  const debounce = (func, wait) => {
    let timeout;
    return (...args) => {
      clearTimeout(timeout);
      timeout = setTimeout(() => func(...args), wait);
    };
  };

  // Handle window resize
  const handleResize = () => {
    if (stageRef.current) {
      stageRef.current.handleResize();
      console.log('Resized NGL Viewer');
    } else {
      console.log('Resize attempted but stage not initialized');
    }
  };

  // Load PDB content into NGL Viewer
  useEffect(() => {
    if (!stageRef.current) {
      console.log('NGL Stage not initialized');
      return;
    }

    const loadStructure = async () => {
      stageRef.current.removeAllComponents();
      console.log(`Cleared all components for ${viewMode}`);

      const pdbContent = viewMode === '2D' ? pdbContent2D : pdbContent3D;
      console.log(`Attempting to load ${viewMode} PDB content:`, pdbContent ? `${pdbContent.slice(0, 100)}...` : 'Empty');

      if (!pdbContent) {
        setViewerError(`No ${viewMode} PDB content available. Please submit a query.`);
        return;
      }

      if (!isValidPdb(pdbContent)) {
        setViewerError(`Invalid ${viewMode} PDB format. Ensure the server returns valid PDB data.`);
        console.error(`Invalid ${viewMode} PDB content:`, pdbContent);
        return;
      }

      setViewerError('');
      const blob = new Blob([pdbContent], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      try {
        const component = await stageRef.current.loadFile(url, { ext: 'pdb' });
        console.log('Loaded component:', component);
        component.addRepresentation('ball+stick', { sele: 'all' });
        component.addRepresentation('spacefill', { sele: 'all', opacity: 0.6 });
        // Fallback representation
        component.addRepresentation('cartoon', { sele: 'all' });
        stageRef.current.autoView();
        console.log(`Successfully loaded ${viewMode} PDB`);
      } catch (err) {
        console.error(`Error loading ${viewMode} PDB:`, err);
        setViewerError(`Failed to load ${viewMode} PDB: ${err.message}`);
      } finally {
        URL.revokeObjectURL(url);
      }
    };

    loadStructure();
  }, [pdbContent2D, pdbContent3D, viewMode]);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!query.trim()) {
      setError('Please enter a chemical name or SMILES string.');
      return;
    }
    setLoading(true);
    setError('');
    setResult(null);
    setPdbContent2D('');
    setPdbContent3D('');

    try {
      const response = await fetch('http://localhost:8000/get_coordinates', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query: query.trim() }),
      });
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      const data = await response.json();
      console.log('API Response:', data);
      console.log('3D PDB Content:', data.pdb_content_3d);
      console.log('2D PDB Content:', data.pdb_content_2d);
      setResult(data);
      setPdbContent2D(data.pdb_content_2d || '');
      setPdbContent3D(data.pdb_content_3d || '');
    } catch (err) {
      setError('Failed to connect to the server. Please ensure the bridge and server are running.');
      console.error('Fetch Error:', err);
    } finally {
      setLoading(false);
    }
  };

  const downloadPDB = (type) => {
    console.log(`downloadPDB called for ${type}`);
    try {
      if (!result || result.status !== 'success') {
        setError('No valid PDB content available.');
        console.log('Error: Invalid result');
        return;
      }
      const content = type === '2D' ? result.pdb_content_2d : result.pdb_content_3d;
      if (!content) {
        setError(`No valid ${type} PDB content available.`);
        console.log(`Error: No ${type} pdb_content`);
        return;
      }
      const name = query.replace(/[^a-zA-Z0-9_-]/g, '') || 'MOLECULE';
      console.log(`${type} PDB Content:`, content);
      if (type === '2D') {
        setPdbContent2D(content);
      } else {
        setPdbContent3D(content);
      }
      const blob = new Blob([content], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      console.log('Blob URL:', url);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${name}_${type.toLowerCase()}_output.pdb`;
      console.log('Triggering download for:', a.download);
      a.click();
      URL.revokeObjectURL(url);
      console.log('Download triggered');
    } catch (err) {
      setError(`Failed to download ${type} PDB file: ${err.message}`);
      console.error(`Download ${type} PDB Error:`, err);
    }
  };

  return (
    <div className="min-h-screen flex flex-col md:flex-row max-w-full w-full p-0 m-0">
      {/* Left Half: Form and Results */}
      <div className="flex-1 bg-white rounded-lg shadow-lg flex flex-col">
        <h1 className="text-2xl font-bold text-center">Chemical Coordinates</h1>
        <p className="text-gray-600 text-center">
          Enter a SMILES string, chemical name, or description to get 2D and 3D coordinates.
        </p>
        <form onSubmit={handleSubmit} className="space-y-2">
          <input
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="e.g., aspirin, CC(=O)Oc1ccccc1C(=O)O"
            className="w-full p-2 border rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
          <button
            type="submit"
            disabled={loading}
            className="w-full bg-blue-500 text-white p-2 rounded-md hover:bg-blue-600 disabled:bg-blue-300"
          >
            {loading ? 'Processing...' : 'Get Coordinates'}
          </button>
        </form>
        {error && <p className="text-red-500 text-center">{error}</p>}
        {result && (
          <div className="bg-gray-50 rounded-md flex-1">
            {result.status === 'success' ? (
              <div>
                <div>
                  <p><strong>SMILES:</strong> {result.gemini_smiles}</p>
                  <p><strong>Source:</strong> {result.source}</p>
                  <p><strong>Message:</strong> {result.gemini_message}</p>
                </div>
                <div className="flex flex-col md:flex-row">
                  <div className="flex-1">
                    <h2 className="text-lg font-semibold">2D Coordinates</h2>
                    <pre className="bg-gray-200 rounded-md overflow-auto max-h-60">
                      {JSON.stringify(result.coordinates_2d, null, 2)}
                    </pre>
                    <button
                      onClick={() => downloadPDB('2D')}
                      className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full"
                    >
                      Download 2D PDB
                    </button>
                    {pdbContent2D && (
                      <div>
                        <p><strong>2D PDB Content:</strong></p>
                        <pre className="bg-gray-200 rounded-md overflow-auto max-h-60">
                          {pdbContent2D}
                        </pre>
                      </div>
                    )}
                  </div>
                  {result.coordinates_3d && (
                    <div className="hidden md:block w-px bg-gray-300"></div>
                  )}
                  {result.coordinates_3d && (
                    <div className="flex-1">
                      <h2 className="text-lg font-semibold">3D Coordinates</h2>
                      <pre className="bg-gray-200 rounded-md overflow-auto max-h-60">
                        {JSON.stringify(result.coordinates_3d, null, 2)}
                      </pre>
                      <button
                        onClick={() => downloadPDB('3D')}
                        className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full"
                      >
                        Download 3D PDB
                      </button>
                      {pdbContent3D && (
                        <div>
                          <p><strong>3D PDB Content:</strong></p>
                          <pre className="bg-gray-200 rounded-md overflow-auto max-h-60">
                            {pdbContent3D}
                          </pre>
                        </div>
                      )}
                    </div>
                  )}
                </div>
              </div>
            ) : (
              <div>
                <p className="text-red-500">
                  <strong>Error:</strong> {result.error}
                </p>
                {result.gemini_message && (
                  <p><strong>Gemini Message:</strong> {result.gemini_message}</p>
                )}
              </div>
            )}
          </div>
        )}
      </div>
      <div className="hidden md:block w-px bg-gray-300"></div>
      <div className="flex-1 flex flex-col bg-white rounded-lg shadow-lg">
        <div className="flex justify-center space-x-2">
          <button
            onClick={() => setViewMode('2D')}
            className={`p-2 rounded-md ${viewMode === '2D' ? 'bg-blue-500 text-white' : 'bg-gray-200'} disabled:opacity-50`}
            disabled={!pdbContent2D}
          >
            View 2D
          </button>
          <button
            onClick={() => setViewMode('3D')}
            className={`p-2 rounded-md ${viewMode === '3D' ? 'bg-blue-500 text-white' : 'bg-gray-200'} disabled:opacity-50`}
            disabled={!pdbContent3D}
          >
            View 3D
          </button>
        </div>
        <div
          ref={viewerRef}
          className="flex-1 w-full min-h-[400px] bg-gray-100"
          style={{ position: 'relative' }}
        >
          {viewerError && (
            <p className="text-red-500 text-center p-4">{viewerError}</p>
          )}
        </div>
      </div>
    </div>
  );
}

export default App;