import { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import Papa from 'papaparse';
import * as NGL from 'ngl';
import './App.css'; // Assuming you have a CSS file for styles

function GeminiChat() {
  const [chatQuery, setChatQuery] = useState('');
  const [chatResponse, setChatResponse] = useState(null);
  const [isChatLoading, setIsChatLoading] = useState(false);

  const handleChatSubmit = async (e) => {
    e.preventDefault();
    if (!chatQuery.trim()) {
      alert('Please enter a chat query');
      return;
    }
    setIsChatLoading(true);
    try {
      const response = await axios.post('http://localhost:8000/gemini-chat', { query: chatQuery.trim() });
      setChatResponse(response.data);
    } catch (error) {
      setChatResponse({
        status: 'error',
        error: error.response?.data?.detail || 'Failed to connect to Gemini chat server',
      });
    } finally {
      setIsChatLoading(false);
    }
  };

  return (
    <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 w-full mb-4">
      <h2 className="text-lg font-semibold mb-3 text-cyan-300">Protien-InfoMania</h2>
      <form onSubmit={handleChatSubmit}>
        <div className="mb-3">
          <label className="block text-gray-200 text-sm mb-1" htmlFor="chatQuery">Get Information</label>
          <input
            type="text"
            id="chatQuery"
            value={chatQuery}
            onChange={(e) => setChatQuery(e.target.value)}
            placeholder="e.g., What is hemoglobin?"
            className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
          />
        </div>
        <button
          type="submit"
          disabled={isChatLoading}
          className={`w-full py-2 bg-cyan-600 text-gray-200 rounded text-sm font-medium ${isChatLoading ? 'opacity-50 cursor-not-allowed' : 'hover:bg-cyan-700'}`}
        >
          {isChatLoading ? 'Sending...' : 'Send'}
        </button>
      </form>
      {isChatLoading && (
        <div className="mt-3">
          <p className="text-gray-200 text-sm mb-2">Waiting for response...</p>
          <div className="w-full bg-gray-700/50 rounded-full h-2">
            <div className="bg-cyan-600 h-2 rounded-full animate-pulse" style={{ width: '50%' }}></div>
          </div>
        </div>
      )}
      {chatResponse && (
        <div className="mt-3">
          <h3 className="text-base font-semibold mb-2 text-cyan-400">Response</h3>
          {chatResponse.status === 'error' ? (
            <p className="text-red-400 text-sm">{chatResponse.error}</p>
          ) : (
            <p className="text-gray-200 text-sm whitespace-pre-wrap">{chatResponse.response}</p>
          )}
        </div>
      )}
    </div>
  );
}

function SearchForm({ onSearch }) {
  const [query, setQuery] = useState('');
  const [limit, setLimit] = useState(10);
  const [isLoading, setIsLoading] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!query.trim()) {
      alert('Please enter a search term');
      return;
    }
    setIsLoading(true);
    await onSearch(query.trim(), parseInt(limit));
    setIsLoading(false);
  };

  return (
    <form onSubmit={handleSubmit} className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 w-full">
      <h2 className="text-lg font-semibold mb-3 text-cyan-300">PDB Search</h2>
      <div className="mb-3">
        <label className="block text-gray-200 text-sm mb-1" htmlFor="query">Search Term</label>
        <input
          type="text"
          id="query"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="e.g., Hemoglobin"
          className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
        />
      </div>
      <div className="mb-3">
        <label className="block text-gray-200 text-sm mb-1" htmlFor="limit">Limit (max results)</label>
        <input
          type="number"
          id="limit"
          value={limit}
          onChange={(e) => setLimit(e.target.value)}
          min="1"
          className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
        />
      </div>
      <button
        type="submit"
        disabled={isLoading}
        className={`w-full py-2 bg-cyan-600 text-gray-200 rounded text-sm font-medium ${isLoading ? 'opacity-50 cursor-not-allowed' : 'hover:bg-cyan-700'}`}
      >
        {isLoading ? 'Searching...' : 'Search'}
      </button>
      {isLoading && (
        <div className="mt-3">
          <p className="text-gray-200 text-sm mb-2">Loading content, please wait...</p>
          <div className="w-full bg-gray-700/50 rounded-full h-2">
            <div className="bg-cyan-600 h-2 rounded-full animate-pulse" style={{ width: '50%' }}></div>
          </div>
        </div>
      )}
    </form>
  );
}

function SearchResults({ result }) {
  if (!result) return null;

  // Parse CSV data
  let tableData = { headers: [], rows: [] };
  if (result.csv_data) {
    try {
      const parsed = Papa.parse(result.csv_data, { header: true });
      tableData.headers = parsed.meta.fields || [];
      tableData.rows = parsed.data || [];
    } catch (error) {
      console.error('Error parsing CSV:', error);
    }
  }

  const handleDownloadCSV = () => {
    if (result.csv_data) {
      const blob = new Blob([result.csv_data], { type: 'text/csv' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${result.gemini_search_term || 'search'}_output.csv`;
      a.click();
      URL.revokeObjectURL(url);
    }
  };

  return (
    <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 w-full mt-4">
      <h2 className="text-lg font-semibold mb-3 text-cyan-300">Search Results</h2>
      {result.status === 'error' ? (
        <div className="text-red-400 text-sm">
          <p><strong>Error:</strong> {result.error}</p>
          {result.gemini_message && <p><strong>Message:</strong> {result.gemini_message}</p>}
        </div>
      ) : (
        <div className="text-gray-200 text-sm">
          <p className="mb-2"><strong>Search Term:</strong> {result.gemini_search_term}</p>
          <p className="mb-2">
            <strong>PDB IDs ({result.pdb_ids?.length || 0}):</strong>{' '}
            {result.pdb_ids?.join(', ') || '-'}
          </p>
          <p className="mb-2">
            <strong>Experimental Data ({result.experimental_data?.length || 0} entries):</strong>
          </p>
          <ul className="list-disc pl-5 mb-3 text-gray-300">
            {result.experimental_data?.map((entry, index) => (
              <li key={index} className="mb-1">
                <span className="text-cyan-400">{entry.rcsb_id}</span>: Method=
                {entry.method || '-'}, Details={entry.details || '-'}
              </li>
            )) || <li>No data available</li>}
          </ul>
          <p className="mb-2"><strong>Message:</strong> {result.gemini_message || '-'}</p>
          {tableData.headers.length > 0 && (
            <div className="mb-3">
              <h3 className="text-base font-semibold mb-2 text-cyan-400">Data Table</h3>
              <div className="overflow-x-auto overflow-y-auto" style={{ maxHeight: '12rem' }}>
                <table className="w-full text-left text-sm text-gray-200">
                  <thead className="bg-gray-700/50 text-cyan-400">
                    <tr>
                      <th className="px-3 py-2">S.No</th>
                      {tableData.headers.map((header, index) => (
                        <th key={index} className="px-3 py-2">
                          {header}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {tableData.rows.map((row, rowIndex) => (
                      <tr key={rowIndex} className="border-t border-gray-600/30">
                        <td className="px-3 py-2">{rowIndex + 1}</td>
                        {tableData.headers.map((header, colIndex) => (
                          <td key={colIndex} className="px-3 py-2">
                            {row[header] || '-'}
                          </td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}
          {result.csv_data && (
            <button
              onClick={handleDownloadCSV}
              className="py-2 px-4 bg-cyan-600 hover:bg-cyan-700 text-gray-200 rounded text-sm font-medium"
            >
              Download CSV
            </button>
          )}
        </div>
      )}
    </div>
  );
}

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
        <p className="text-gray-300 text-sm">Select a PDB ID to view 3D structure</p>
      )}
    </div>
  );
}

function App() {
  const [result, setResult] = useState(null);
  const [selectedPdbId, setSelectedPdbId] = useState('');
  const [pdbContent, setPdbContent] = useState('');
  const [pdbFetchError, setPdbFetchError] = useState('');
  const [isPdbLoading, setIsPdbLoading] = useState(false);
  const [viewerFeatures, setViewerFeatures] = useState({
    backbone: false,
    cartoon: true,
    line: false,
    ballAndStick: false,
    label: false,
  });

  const handleSearch = async (query, limit) => {
    try {
      const response = await axios.post('http://localhost:8000/search', { query, limit });
      setResult(response.data);
      setSelectedPdbId('');
      setPdbContent('');
      setPdbFetchError('');
    } catch (error) {
      setResult({
        status: 'error',
        error: error.response?.data?.detail || 'Failed to connect to server',
      });
      setSelectedPdbId('');
      setPdbContent('');
      setPdbFetchError('');
    }
  };

  // Fetch PDB content for the selected PDB ID
  useEffect(() => {
    const fetchPdbContent = async () => {
      if (selectedPdbId) {
        setIsPdbLoading(true);
        try {
          const response = await axios.get(`http://localhost:8000/fetch-pdb/${selectedPdbId}`);
          setPdbContent(response.data.pdb_content);
          setPdbFetchError('');
          console.log('Fetched PDB content for:', selectedPdbId);
        } catch (error) {
          console.error('Error fetching PDB content:', error);
          setPdbContent('');
          setPdbFetchError(error.response?.data?.detail || 'Failed to fetch PDB file');
        } finally {
          setIsPdbLoading(false);
        }
      } else {
        setPdbContent('');
        setPdbFetchError('');
        setIsPdbLoading(false);
      }
    };

    fetchPdbContent();
  }, [selectedPdbId]);

  // Handle PDB file download
  const handleDownloadPdb = () => {
    if (pdbContent && selectedPdbId) {
      const blob = new Blob([pdbContent], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${selectedPdbId}.pdb`;
      a.click();
      URL.revokeObjectURL(url);
    }
  };

  return (
    <div className="flex flex-row w-full min-h-screen">
      <div className="w-2/5 p-4 flex flex-col items-center overflow-y-auto">
        <GeminiChat />
        <SearchForm onSearch={handleSearch} />
        {result && <SearchResults result={result} />}
      </div>
      <div className="w-3/5 p-4 flex flex-col bg-gray-800/90">
        {/* PDB ID Scrollable List */}
        <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 mb-4 max-h-32 overflow-y-auto">
          <h3 className="text-base font-semibold mb-2 text-cyan-300">PDB IDs</h3>
          {result?.status === 'success' && result.pdb_ids?.length > 0 ? (
            <div className="grid grid-cols-4 gap-2">
              {result.pdb_ids.map((pdbId, index) => (
                <button
                  key={pdbId}
                  onClick={() => setSelectedPdbId(pdbId)}
                  className={`p-2 rounded text-sm flex items-center ${selectedPdbId === pdbId
                      ? 'bg-cyan-600 text-gray-200'
                      : 'bg-gray-700/50 text-gray-300 hover:bg-gray-600/50'
                    }`}
                >
                  <span className="mr-2">{index + 1}.</span>
                  <span>{pdbId}</span>
                </button>
              ))}
            </div>
          ) : (
            <p className="text-gray-300 text-sm">No PDB IDs available</p>
          )}
        </div>
        {/* PDB Content Display */}
        <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 mb-4">
          <h3 className="text-base font-semibold mb-2 text-cyan-300">
            PDB Content: {selectedPdbId || 'None Selected'}
          </h3>
          {isPdbLoading ? (
            <p className="text-gray-200 text-sm">Content is being loaded...</p>
          ) : pdbFetchError ? (
            <p className="text-red-400 text-sm">{pdbFetchError}</p>
          ) : pdbContent ? (
            <>
              <pre className="text-gray-200 text-xs overflow-y-auto" style={{ maxHeight: '12rem' }}>
                {pdbContent}
              </pre>
              <button
                onClick={handleDownloadPdb}
                className="mt-2 py-2 px-4 bg-cyan-600 hover:bg-cyan-700 text-gray-200 rounded text-sm font-medium"
              >
                Download PDB
              </button>
            </>
          ) : (
            <p className="text-gray-300 text-sm">Select a PDB ID to view content</p>
          )}
        </div>
        {/* NGL Viewer */}
        <NGLViewer
          pdbContent={pdbContent}
          viewerFeatures={viewerFeatures}
          isPdbLoading={isPdbLoading}
          setViewerFeatures={setViewerFeatures}
        />
      </div>
    </div>
  );
}

export default App;