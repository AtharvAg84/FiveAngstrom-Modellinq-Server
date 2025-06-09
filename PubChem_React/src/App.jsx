import { useState, useEffect, useRef } from 'react';
import * as NGL from 'ngl';

function App() {
  const [query, setQuery] = useState('');
  const [result, setResult] = useState(null);
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [pdbContent2D, setPdbContent2D] = useState('');
  const [pdbContent3D, setPdbContent3D] = useState('');
  const [image2D, setImage2D] = useState('');
  const [viewMode, setViewMode] = useState('3D');
  const [viewerError, setViewerError] = useState('');
  const [features, setFeatures] = useState({
    backbone: false,
    cartoon: false,
    line: false,
    ballAndStick: false,
    label: false,
  });
  const stageRef = useRef(null);
  const viewerRef = useRef(null);
  const componentRef = useRef(null);

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
        } catch (err) {
          console.error('NGL Stage initialization error:', err);
          setViewerError('Failed to initialize NGL Viewer: ' + err.message);
        }
      }
    };

    const timer = setTimeout(initNGL, 100);
    return () => {
      clearTimeout(timer);
      if (stageRef.current) {
        stageRef.current.removeAllComponents();
        stageRef.current.dispose();
        stageRef.current = null;
      }
      componentRef.current = null;
    };
  }, []);

  // Toggle feature state
  const toggleFeature = (feature) => {
    setFeatures((prev) => ({
      ...prev,
      [feature]: !prev[feature],
    }));
  };

  // Update representations without reloading PDB
  const updateRepresentations = (component, isProtein) => {
    if (!component) return;
    component.removeAllRepresentations();
    console.log('Cleared all representations');

    if (viewMode === '3D') {
      if (features.backbone) {
        component.addRepresentation('backbone', {
          sele: isProtein ? 'protein' : 'all',
          color: 'blue',
          radius: 0.5,
          linewidth: 1.0,
          aspectRatio: 1.0,
        });
        console.log('Applied backbone representation');
      }
      if (features.cartoon) {
        component.addRepresentation('cartoon', {
          sele: isProtein ? 'protein' : 'all',
          color: 'green',
          opacity: 0.8,
        });
        console.log('Applied cartoon representation');
      }
      if (features.line) {
        component.addRepresentation('line', {
          sele: 'all',
          color: 'black',
          linewidth: 3,
        });
        console.log('Applied line representation');
      }
      if (features.ballAndStick) {
        component.addRepresentation('ball+stick', {
          sele: 'all',
          color: 'element',
          radius: 0.2,
          aspectRatio: 2.0,
        });
        console.log('Applied ball+stick representation');
      }
      if (features.label) {
        component.addRepresentation('label', {
          sele: 'all',
          labelType: 'atomname',
          color: 'black',
          fontSize: 0.5,
          showBackground: true,
          backgroundColor: 'white',
          backgroundOpacity: 0.8,
        });
        console.log('Applied label representation');
      }
      if (!Object.values(features).some((v) => v)) {
        component.addRepresentation('spacefill', {
          sele: 'all',
          color: 'element',
          opacity: 0.6,
          radius: 0.9,
        });
        console.log('Applied fallback spacefill representation');
      }
    }
    if (stageRef.current) {
      stageRef.current.autoView();
      console.log('Auto-view triggered');
    }
  };

  // Load PDB content into NGL Viewer (for 3D)
  useEffect(() => {
    if (!stageRef.current) {
      console.log('NGL Stage not initialized');
      return;
    }

    const loadStructure = async () => {
      console.log(`Loading structure for ${viewMode}`);
      stageRef.current.removeAllComponents();
      componentRef.current = null;
      console.log('Cleared all components and componentRef');

      if (viewMode !== '3D') {
        console.log('Not in 3D mode, skipping structure load');
        return;
      }

      const pdbContent = pdbContent3D;
      console.log(`Attempting to load 3D PDB content:`, pdbContent ? `${pdbContent.slice(0, 100)}...` : 'Empty');

      if (!pdbContent) {
        setViewerError('No 3D PDB content available. Please submit a query.');
        return;
      }

      if (!isValidPdb(pdbContent)) {
        setViewerError('Invalid 3D PDB format. Ensure the server returns valid PDB data.');
        console.error('Invalid 3D PDB content:', pdbContent);
        return;
      }

      setViewerError('');
      const blob = new Blob([pdbContent], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      try {
        console.log('Loading PDB file into NGL');
        const component = await stageRef.current.loadFile(url, { ext: 'pdb' });
        console.log('Loaded component:', component);
        componentRef.current = component;

        const isProtein = pdbContent.includes('HELIX') || pdbContent.includes('SHEET') ||
          pdbContent.split('\n').filter(line => line.startsWith('ATOM')).length > 50;
        console.log(`Is protein: ${isProtein}`);

        updateRepresentations(component, isProtein);
      } catch (err) {
        console.error('Error loading 3D PDB:', err);
        setViewerError(`Failed to load 3D PDB: ${err.message}`);
      } finally {
        URL.revokeObjectURL(url);
        console.log('Revoked Blob URL');
      }
    };

    loadStructure();
  }, [pdbContent3D, viewMode]);

  // Update representations when features change
  useEffect(() => {
    if (componentRef.current && stageRef.current && viewMode === '3D') {
      const isProtein = pdbContent3D.includes('HELIX') || pdbContent3D.includes('SHEET') ||
        pdbContent3D.split('\n').filter(line => line.startsWith('ATOM')).length > 50;
      updateRepresentations(componentRef.current, isProtein);
    }
  }, [features, viewMode, pdbContent3D]);

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
    setImage2D('');

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
      console.log('2D Image:', data.image_2d ? 'Present' : 'Absent');
      setResult(data);
      setPdbContent2D(data.pdb_content_2d || '');
      setPdbContent3D(data.pdb_content_3d || '');
      setImage2D(data.image_2d || '');
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

  const downloadImage2D = () => {
    console.log('downloadImage2D called');
    try {
      if (!image2D) {
        setError('No 2D image available.');
        console.log('Error: No 2D image');
        return;
      }
      const name = query.replace(/[^a-zA-Z0-9_-]/g, '') || 'MOLECULE';
      const a = document.createElement('a');
      a.href = `data:image/png;base64,${image2D}`;
      a.download = `${name}_2d_image.png`;
      console.log('Triggering download for:', a.download);
      a.click();
      console.log('2D Image download triggered');
    } catch (err) {
      setError(`Failed to download 2D image: ${err.message}`);
      console.error('Download 2D Image Error:', err);
    }
  };

  return (
    <div className="min-h-screen flex flex-col md:flex-row max-w-full w-full p-0 m-0">
      {/* Left Half: Form and Results */}
      <div className="flex-1 bg-white rounded-lg shadow-lg flex flex-col p-4">
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
            {loading ? 'Processing...' : 'Ready. Steady.. Go...'}
          </button>
        </form>
        {error && <p className="text-red-500 text-center">{error}</p>}
        {result && (
          <div className="bg-gray-50 rounded-md flex-1 p-4">
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
                    <pre className="bg-gray-200 rounded-md overflow-auto max-h-60 p-2">
                      {JSON.stringify(result.coordinates_2d, null, 2)}
                    </pre>
                    <br />
                    <br />
                    {pdbContent2D && (
                      <div>
                        <p><strong>2D PDB Content:</strong></p>
                        <pre className="bg-gray-200 rounded-md overflow-auto max-h-60 p-2">
                          {pdbContent2D}
                        </pre>
                      </div>
                    )}
                    <button
                      onClick={() => downloadPDB('2D')}
                      className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full mt-2"
                    >
                      Download 2D PDB
                    </button>

                    <button
                      onClick={downloadImage2D}
                      className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full mt-2"
                      disabled={!image2D}
                    >
                      Download 2D Image
                    </button>
                  </div>
                  {result.coordinates_3d && (
                    <div className="hidden md:block w-px bg-gray-300"></div>
                  )}
                  {result.coordinates_3d && (
                    <div className="flex-1">
                      <h2 className="text-lg font-semibold">3D Coordinates</h2>
                      <pre className="bg-gray-200 rounded-md overflow-auto max-h-60 p-2">
                        {JSON.stringify(result.coordinates_3d, null, 2)}
                      </pre>
                      <br />
                      <br />
                      {pdbContent3D && (
                        <div>
                          <p><strong>3D PDB Content:</strong></p>
                          <pre className="bg-gray-200 rounded-md overflow-auto max-h-60 p-2">
                            {pdbContent3D}
                          </pre>
                        </div>
                      )}
                      <button
                        onClick={() => downloadPDB('3D')}
                        className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full mt-2"
                      >
                        Download 3D PDB
                      </button>
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
      <div className="flex-1 flex flex-col bg-white rounded-lg shadow-lg p-4">
        <div className="flex flex-wrap justify-center space-x-2 mb-2">
          <button
            onClick={() => setViewMode('2D')}
            className={`p-2 rounded-md ${viewMode === '2D' ? 'bg-blue-500 text-white' : 'bg-gray-200'} disabled:opacity-50`}
            disabled={!image2D && !pdbContent2D}
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
          {viewMode === '3D' && (
            <>
              {/* <button
                onClick={() => toggleFeature('backbone')}
                className={`p-2 rounded-md ${features.backbone ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
              >
                Backbone
              </button>
              <button
                onClick={() => toggleFeature('cartoon')}
                className={`p-2 rounded-md ${features.cartoon ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
              >
                Cartoon
              </button> */}
              <button
                onClick={() => toggleFeature('line')}
                className={`p-2 rounded-md ${features.line ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
              >
                Line
              </button>
              <button
                onClick={() => toggleFeature('ballAndStick')}
                className={`p-2 rounded-md ${features.ballAndStick ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
              >
                Ball+Stick
              </button>
              <button
                onClick={() => toggleFeature('label')}
                className={`p-2 rounded-md ${features.label ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
              >
                Label
              </button>
            </>
          )}
        </div>
        <div
          className="flex-1 w-full bg-gray-100"
          style={{ position: 'relative', width: '600px', height: '400px' }}
        >
          <div
            className="w-full h-full"
            style={{ display: viewMode === '2D' ? 'none' : 'block' }}
            ref={viewerRef}
          >
            {viewerError && viewMode === '3D' && (
              <p className="text-red-500 text-center p-4">{viewerError}</p>
            )}
          </div>
          <div
            className="w-full h-full flex items-center justify-center"
            style={{ display: viewMode === '2D' ? 'block' : 'none' }}
          >
            {image2D ? (
              <img
                src={`data:image/png;base64,${image2D}`}
                alt="2D Molecule"
                className="max-w-full max-h-[400px] m-auto"
              />
            ) : (
              <p className="text-red-500 text-center p-4">No 2D image available. Please submit a query.</p>
            )}
          </div>
        </div>
      </div>
    </div>
  );
}

export default App;