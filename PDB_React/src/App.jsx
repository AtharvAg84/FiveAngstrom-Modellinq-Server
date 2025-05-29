import { useState } from 'react';
import axios from 'axios';

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
    <form onSubmit={handleSubmit} className="bg-blue-900/90 p-4 rounded-md border border-blue-700/50 w-full">
      <h2 className="text-xl font-semibold mb-3 text-cyan-300">PDB Search</h2>
      <div className="mb-3">
        <label className="block text-blue-200 text-xs mb-1" htmlFor="query">Search Term</label>
        <input
          type="text"
          id="query"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="e.g., Hemoglobin"
          className="w-full p-2 bg-blue-950/50 text-blue-100 border border-blue-600/50 rounded text-xs"
        />
      </div>
      <div className="mb-3">
        <label className="block text-blue-200 text-xs mb-1" htmlFor="limit">Limit (max results)</label>
        <input
          type="number"
          id="limit"
          value={limit}
          onChange={(e) => setLimit(e.target.value)}
          min="1"
          className="w-full p-2 bg-blue-950/50 text-blue-100 border border-blue-600/50 rounded text-xs"
        />
      </div>
      <button
        type="submit"
        disabled={isLoading}
        className={`w-full py-2 bg-cyan-600 text-white rounded text-xs font-medium ${
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
    <div className="bg-blue-900/90 p-4 rounded-md border border-blue-700/50 w-full mt-4">
      <h2 className="text-xl font-semibold mb-3 text-cyan-300">Search Results</h2>
      {result.status === 'error' ? (
        <div className="text-red-400 text-sm">
          <p><strong>Error:</strong> {result.error}</p>
          {result.gemini_message && <p><strong>Message:</strong> {result.gemini_message}</p>}
        </div>
      ) : (
        <div className="text-blue-100 text-sm">
          <p className="mb-2"><strong>Search Term:</strong> {result.gemini_search_term}</p>
          <p className="mb-2">
            <strong>PDB IDs ({result.pdb_ids?.length || 0}):</strong>{' '}
            {result.pdb_ids?.join(', ') || '-'}
          </p>
          <p className="mb-2">
            <strong>Experimental Data ({result.experimental_data?.length || 0} entries):</strong>
          </p>
          <ul className="list-disc pl-5 mb-3 text-blue-200">
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
              <div className="max-h-48 overflow-y-auto border border-blue-600/50 rounded">
                <table className="w-full text-left text-xs text-blue-100">
                  <thead className="bg-blue-800/50 text-cyan-400 sticky top-0">
                    <tr>
                      {tableData.headers.map((header, index) => (
                        <th key={index} className="px-2 py-1">
                          {header}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {tableData.rows.map((row, rowIndex) => (
                      <tr key={rowIndex} className="border-t border-blue-600/30 hover:bg-blue-800/30">
                        {tableData.headers.map((header, colIndex) => (
                          <td key={colIndex} className="px-2 py-1">
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
              className="py-2 px-4 bg-cyan-600 hover:bg-cyan-700 text-white rounded text-xs font-medium"
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

  return (
    <div className="flex flex-row w-full min-h-screen">
      <div className="w-1/2 p-4 flex flex-col items-center overflow-y-auto">
        <SearchForm onSearch={handleSearch} />
        {result && <SearchResults result={result} />}
      </div>
      <div className="w-1/2 p-4 flex items-center justify-center bg-blue-900/90">
        <p className="text-blue-200 text-sm">Placeholder for future content</p>
      </div>
    </div>
  );
}

export default App;