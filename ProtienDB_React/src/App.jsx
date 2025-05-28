import React, { useState, useEffect, useRef } from 'react';
import * as NGL from 'ngl';
import axios from 'axios';
import './index.css'; // Assuming Tailwind CSS is set up here

const App = () => {
  const [query, setQuery] = useState('');
  const [limit, setLimit] = useState(10);
  const [results, setResults] = useState([]);
  const [message, setMessage] = useState('');
  const [error, setError] = useState('');
  const [selectedPdbId, setSelectedPdbId] = useState('');
  const stageRef = useRef(null);
  const viewerRef = useRef(null);

  // Initialize NGL Viewer
  useEffect(() => {
    if (!viewerRef.current) {
      stageRef.current = new NGL.Stage('ngl-viewer', { backgroundColor: 'white' });
      viewerRef.current = stageRef.current;
      // Handle window resize
      window.addEventListener('resize', () => stageRef.current.handleResize());
    }
    return () => {
      if (viewerRef.current) {
        viewerRef.current.dispose();
        viewerRef.current = null;
      }
    };
  }, []);

  // Load PDB structure when selectedPdbId changes
  useEffect(() => {
    if (selectedPdbId && viewerRef.current) {
      viewerRef.current.removeAllComponents();
      const pdbUrl = `https://files.rcsb.org/download/${selectedPdbId}.pdb`;
      viewerRef.current
        .loadFile(pdbUrl, { ext: 'pdb' })
        .then((component) => {
          component.addRepresentation('cartoon');
          component.addRepresentation('licorice', { sele: 'hetero' });
          component.autoView();
        })
        .catch((err) => {
          setError(`Failed to load PDB ${selectedPdbId}: ${err.message}`);
        });
    }
  }, [selectedPdbId]);

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    setError('');
    setMessage('');
    setResults([]);
    setSelectedPdbId('');

    try {
      const response = await axios.post('http://localhost:8000/search_pdb', {
        query,
        limit,
      });
      const data = response.data;

      if (data.status === 'error') {
        setError(data.error);
        return;
      }

      setMessage(data.message);
      if (data.csv_data) {
        // Parse CSV data
        const rows = data.csv_data.split('\n').map((row) => row.split(','));
        const headers = rows[0];
        const parsedResults = rows.slice(1).map((row) =>
          headers.reduce((obj, header, index) => {
            obj[header] = row[index] || '-';
            return obj;
          }, {})
        );
        setResults(parsedResults);
        if (parsedResults.length > 0) {
          setSelectedPdbId(parsedResults[0].rcsb_id); // Auto-select first PDB ID
        }
      }
    } catch (err) {
      setError(`Failed to fetch data: ${err.message}`);
    }
  };

  return (
    <div className="container mx-auto p-4">
      <h1 className="text-3xl font-bold mb-4">PDB Search and Visualization</h1>

      {/* Form */}
      <form onSubmit={handleSubmit} className="mb-6">
        <div className="flex flex-col md:flex-row gap-4">
          <div className="flex-1">
            <label className="block text-sm font-medium text-gray-700">
              Search Term
            </label>
            <input
              type="text"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              placeholder="e.g., Hemoglobin"
              className="mt-1 block w-full border border-gray-300 rounded-md p-2 focus:outline-none focus:ring-2 focus:ring-blue-500"
              required
            />
          </div>
          <div className="w-32">
            <label className="block text-sm font-medium text-gray-700">
              Limit
            </label>
            <input
              type="number"
              value={limit}
              onChange={(e) => setLimit(Number(e.target.value))}
              min="1"
              max="100"
              className="mt-1 block w-full border border-gray-300 rounded-md p-2 focus:outline-none focus:ring-2 focus:ring-blue-500"
              required
            />
          </div>
        </div>
        <button
          type="submit"
          className="mt-4 bg-blue-500 text-white px-4 py-2 rounded-md hover:bg-blue-600"
        >
          Search
        </button>
      </form>

      {/* Error/Message */}
      {error && <p className="text-red-500 mb-4">{error}</p>}
      {message && <p className="text-green-500 mb-4">{message}</p>}

      {/* Results Table */}
      {results.length > 0 && (
        <div className="mb-6">
          <h2 className="text-xl font-semibold mb-2">Results</h2>
          <div className="overflow-x-auto">
            <table className="min-w-full border-collapse border border-gray-300">
              <thead>
                <tr className="bg-gray-100">
                  {Object.keys(results[0]).map((header) => (
                    <th
                      key={header}
                      className="border border-gray-300 px-4 py-2 text-left"
                    >
                      {header}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {results.map((row, index) => (
                  <tr
                    key={index}
                    className={`${
                      row.rcsb_id === selectedPdbId ? 'bg-blue-100' : ''
                    } hover:bg-gray-50 cursor-pointer`}
                    onClick={() => setSelectedPdbId(row.rcsb_id)}
                  >
                    {Object.values(row).map((value, i) => (
                      <td
                        key={i}
                        className="border border-gray-300 px-4 py-2"
                      >
                        {value}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* NGL Viewer */}
      <div className="mb-6">
        <h2 className="text-xl font-semibold mb-2">PDB Visualization</h2>
        <div
          id="ngl-viewer"
          className="border border-gray-300"
          style={{ width: '100%', height: '400px' }}
        ></div>
        {selectedPdbId && (
          <p className="mt-2 text-sm text-gray-600">
            Viewing PDB ID: {selectedPdbId}
          </p>
        )}
      </div>
    </div>
  );
};

export default App;