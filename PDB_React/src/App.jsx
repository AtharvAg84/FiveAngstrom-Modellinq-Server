import { useState, useEffect } from 'react';
import axios from 'axios';
import NGLViewer from './NGLViewer';

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
        className={`w-full py-2 bg-cyan-600 text-gray-200 rounded text-sm font-medium ${
          isLoading ? 'opacity-50 cursor-not-allowed' : 'hover:bg-cyan-700'
        }`}
      >
        {isLoading ? 'Searching...' : 'Search'}
      </button>
    </form>
  );
}

function SearchResults({ result }) {
  if (!result) return null;

  // Parse CSV data
  let tableData = { headers: [], rows: [] };
  if (result.csv_data) {
    try {
      const parsed = window.Papa.parse(result.csv_data, { header: true });
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
              <div className="overflow-y-auto" style={{ maxHeight: '12rem' }}>
                <table className="w-full text-left text-sm text-gray-200">
                  <thead className="bg-gray-700/50 text-cyan-400">
                    <tr>
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

function App() {
  const [result, setResult] = useState(null);
  const [pdbContent, setPdbContent] = useState('');

  const handleSearch = async (query, limit) => {
    try {
      const response = await axios.post('http://localhost:8000/search', { query, limit });
      setResult(response.data);
    } catch (error) {
      setResult({
        status: 'error',
        error: error.response?.data?.detail || 'Failed to connect to server',
      });
    }
  };

  // Fetch PDB content client-side for the first PDB ID
  useEffect(() => {
    const fetchPdbContent = async () => {
      if (result?.status === 'success' && result.pdb_ids?.length > 0) {
        const pdbId = result.pdb_ids[0];
        try {
          const response = await axios.get(`https://files.rcsb.org/download/${pdbId}.pdb`);
          setPdbContent(response.data);
        } catch (error) {
          console.error('Error fetching PDB content:', error);
          setPdbContent('');
        }
      } else {
        setPdbContent('');
      }
    };

    fetchPdbContent();
  }, [result]);

  return (
    <div className="flex flex-row w-full min-h-screen">
      <div className="w-1/2 p-4 flex flex-col items-center overflow-y-auto">
        <SearchForm onSearch={handleSearch} />
        {result && <SearchResults result={result} />}
      </div>
      <div className="w-1/2 p-4 flex items-center justify-center bg-gray-800/90">
        <NGLViewer
          pdbContent3D={pdbContent}
          viewMode="3D"
          features={{ backbone: true, cartoon: true }}
        />
      </div>
    </div>
  );
}

export default App;