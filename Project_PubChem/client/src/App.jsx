import { useState, useEffect } from 'react';
import NGLViewer from './NGLViewer';

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

  // Toggle feature state
  const toggleFeature = (feature) => {
    setFeatures((prev) => ({
      ...prev,
      [feature]: !prev[feature],
    }));
  };

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
        <h1 className="text-2xl font-bold text-center">Chemical Structure</h1>
        <p className="text-gray-600 text-center">
          Enter a SMILES string, chemical name, or description to get 2D and 3D coordinates.
        </p>
        <form onSubmit={handleSubmit} className="space-y-2">
          <input
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="e.g., aspirin, CC(=O)Oc1ccccc1C(=O)O"
            className="w-full p-2 border rounded-md focus:outline-none focus:ring-2 focus:ring-purple-500"
          />
          <button
            type="submit"
            disabled={loading}
            className="w-full bg-purple-500 text-white p-2 rounded-md hover:bg-purple-600 disabled:bg-purple-300"
          >
            {loading ? 'Processing...' : 'Generate the Structure'}
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
                      className="bg-purple-500 text-white p-2 rounded-md hover:bg-purple-600 w-full mt-2"
                    >
                      Download 2D PDB
                    </button>

                    <button
                      onClick={downloadImage2D}
                      className="bg-purple-500 text-white p-2 rounded-md hover:bg-purple-600 w-full mt-2"
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
                        className="bg-purple-500 text-white p-2 rounded-md hover:bg-purple-600 w-full mt-2"
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
            className={`p-2 rounded-md ${viewMode === '2D' ? 'bg-purple-500 text-white' : 'bg-gray-200'} disabled:opacity-50`}
            disabled={!image2D && !pdbContent2D}
          >
            View 2D
          </button>
          <button
            onClick={() => setViewMode('3D')}
            className={`p-2 rounded-md ${viewMode === '3D' ? 'bg-purple-500 text-white' : 'bg-gray-200'} disabled:opacity-50`}
            disabled={!pdbContent3D}
          >
            View 3D
          </button>
          {viewMode === '3D' && (
            <>
              <button
                onClick={() => toggleFeature('line')}
                className={`p-2 rounded-md ${features.line ? 'bg-purple-500 text-white' : 'bg-gray-200'}`}
              >
                Line
              </button>
              <button
                onClick={() => toggleFeature('ballAndStick')}
                className={`p-2 rounded-md ${features.ballAndStick ? 'bg-purple-500 text-white' : 'bg-gray-200'}`}
              >
                Ball+Stick
              </button>
              <button
                onClick={() => toggleFeature('label')}
                className={`p-2 rounded-md ${features.label ? 'bg-purple-500 text-white' : 'bg-gray-200'}`}
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
          <NGLViewer
            pdbContent3D={pdbContent3D}
            viewMode={viewMode}
            features={features}
            setViewerError={setViewerError}
          />
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