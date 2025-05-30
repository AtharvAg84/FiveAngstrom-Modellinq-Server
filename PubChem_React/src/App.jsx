import { useState } from 'react';
import ChemicalForm from './components/ChemicalForm';
import ResultDisplay from './components/ResultDisplay';
import MoleculeViewer from './components/MoleculeViewer';

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

  return (
    <div className="min-h-screen flex flex-col md:flex-row max-w-full w-full p-0 m-0">
      <div className="flex-1 bg-white rounded-lg shadow-lg flex flex-col p-4">
        <h1 className="text-2xl font-bold text-center">Chemical Coordinates</h1>
        <p className="text-gray-600 text-center">
          Enter a SMILES string, chemical name, or description to get 2D and 3D coordinates.
        </p>
        <ChemicalForm
          query={query}
          setQuery={setQuery}
          loading={loading}
          setLoading={setLoading}
          setError={setError}
          setResult={setResult}
          setPdbContent2D={setPdbContent2D}
          setPdbContent3D={setPdbContent3D}
          setImage2D={setImage2D}
        />
        {error && <p className="text-red-500 text-center">{error}</p>}
        <ResultDisplay
          result={result}
          pdbContent2D={pdbContent2D}
          pdbContent3D={pdbContent3D}
          query={query}
        />
      </div>
      <div className="hidden md:block w-px bg-gray-300"></div>
      <MoleculeViewer
        viewMode={viewMode}
        setViewMode={setViewMode}
        features={features}
        setFeatures={setFeatures}
        image2D={image2D}
        pdbContent3D={pdbContent3D}
        viewerError={viewerError}
        setViewerError={setViewerError}
      />
    </div>
  );
}

export default App;