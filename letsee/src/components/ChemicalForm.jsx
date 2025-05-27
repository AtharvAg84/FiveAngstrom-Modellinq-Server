import { handleSubmit } from '../utils/apiUtils';

function ChemicalForm({
  query,
  setQuery,
  loading,
  setLoading,
  setError,
  setResult,
  setPdbContent2D,
  setPdbContent3D,
  setImage2D,
}) {
  return (
    <form onSubmit={(e) => handleSubmit(e, query, setLoading, setError, setResult, setPdbContent2D, setPdbContent3D, setImage2D)} className="space-y-2">
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
        {loading ? 'Processing...' : 'Submit'}
      </button>
    </form>
  );
}

export default ChemicalForm;