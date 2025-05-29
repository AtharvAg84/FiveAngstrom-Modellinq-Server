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
    <form
      onSubmit={handleSubmit}
      className="bg-gradient-to-br from-blue-900/80 to-blue-800/80 p-8 rounded-xl shadow-2xl w-full max-w-md backdrop-blur-md border border-blue-700/50 flex flex-col justify-center h-full"
    >
      <h2 className="text-3xl font-bold mb-6 text-cyan-300 tracking-tight">PDB Search</h2>
      <div className="mb-6">
        <label className="block text-blue-200 mb-2 text-sm font-medium" htmlFor="query">
          Search Term
        </label>
        <input
          type="text"
          id="query"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="e.g., Hemoglobin"
          className="w-full p-3 bg-blue-950/50 text-blue-100 border border-blue-600/50 rounded-lg focus:outline-none focus:ring-2 focus:ring-cyan-400 transition-all duration-300 glow-on-hover"
        />
      </div>
      <div className="mb-6">
        <label className="block text-blue-200 mb-2 text-sm font-medium" htmlFor="limit">
          Limit (max results)
        </label>
        <input
          type="number"
          id="limit"
          value={limit}
          onChange={(e) => setLimit(e.target.value)}
          min="1"
          className="w-full p-3 bg-blue-950/50 text-blue-100 border border-blue-600/50 rounded-lg focus:outline-none focus:ring-2 focus:ring-cyan-400 transition-all duration-300 glow-on-hover"
        />
      </div>
      <button
        type="submit"
        disabled={isLoading}
        className={`w-full py-3 px-4 rounded-lg text-white font-semibold transition-all duration-300 glow-on-hover ${
          isLoading
            ? 'bg-blue-600/50 cursor-not-allowed'
            : 'bg-cyan-600 hover:bg-cyan-700 hover:shadow-lg hover:shadow-cyan-400/30'
        }`}
      >
        {isLoading ? (
          <span className="flex items-center justify-center">
            <svg
              className="animate-spin h-5 w-5 mr-2 text-cyan-300"
              xmlns="http://www.w3.org/2000/svg"
              fill="none"
              viewBox="0 0 24 24"
            >
              <circle
                className="opacity-25"
                cx="12"
                cy="12"
                r="10"
                stroke="currentColor"
                strokeWidth="4"
              ></circle>
              <path
                className="opacity-75"
                fill="currentColor"
                d="M4 12a8 8 0 018-8v8h8a8 8 0 01-16 0z"
              ></path>
            </svg>
            Searching...
          </span>
        ) : (
          'Search'
        )}
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
    <div className="bg-gradient-to-br from-blue-900/80 to-blue-800/80 p-8 rounded-xl shadow-2xl w-full max-w-2xl backdrop-blur-md border border-blue-700/50 overflow-y-auto animate-fade-in h-full flex flex-col">
      <h2 className="text-3xl font-bold mb-6 text-cyan-300 tracking-tight">Search Results</h2>
      {result.status === 'error' ? (
        <div className="text-red-400">
          <p className="text-lg"><strong>Error:</strong> {result.error}</p>
          {result.gemini_message && (
            <p className="text-lg"><strong>Message:</strong> {result.gemini_message}</p>
          )}
        </div>
      ) : (
        <div className="text-blue-100 flex-grow">
          <p className="mb-4 text-lg">
            <strong>Search Term:</strong> {result.gemini_search_term}
          </p>
          <p className="mb-4 text-lg">
            <strong>PDB IDs ({result.pdb_ids?.length || 0}):</strong>{' '}
            {result.pdb_ids?.join(', ') || '-'}
          </p>
          <p className="mb-4 text-lg">
            <strong>Experimental Data ({result.experimental_data?.length || 0} entries):</strong>
          </p>
          <ul className="list-disc pl-6 mb-6 text-blue-200">
            {result.experimental_data?.map((entry, index) => (
              <li key={index} className="mb-2">
                <span className="text-cyan-400">{entry.rcsb_id}</span>: Method=
                {entry.method || '-'}, Details={entry.details || '-'}
              </li>
            )) || <li>No data available</li>}
          </ul>
          <p className="mb-4 text-lg">
            <strong>Message:</strong> {result.gemini_message || '-'}
          </p>
          {tableData.headers.length > 0 && (
            <div className="mb-6">
              <h3 className="text-xl font-semibold mb-4 text-cyan-300">Data Table</h3>
              <div className="max-h-64 overflow-y-auto rounded-lg border border-blue-700/50">
                <table className="w-full text-left text-blue-100">
                  <thead className="bg-blue-950/50 text-cyan-300 sticky top-0">
                    <tr>
                      {tableData.headers.map((header, index) => (
                        <th key={index} className="px-4 py-2 font-medium">
                          {header}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {tableData.rows.map((row, rowIndex) => (
                      <tr
                        key={rowIndex}
                        className="border-t border-blue-700/30 hover:bg-blue-950/30 glow-on-hover transition-all duration-300"
                      >
                        {tableData.headers.map((header, colIndex) => (
                          <td key={colIndex} className="px-4 py-2">
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
              className="mt-4 py-3 px-6 bg-cyan-600 hover:bg-cyan-700 text-white rounded-lg font-semibold transition-all duration-300 glow-on-hover"
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
    <div className="min-h-screen flex flex-row w-full">
      <div className="w-1/2 p-6 flex items-center justify-center">
        <SearchForm onSearch={handleSearch} />
      </div>
      <div className="w-1/2 p-6 flex items-start justify-center">
        <SearchResults result={result} />
      </div>
    </div>
  );
}

export default App;