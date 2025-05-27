import { downloadPDB, downloadImage2D } from '../utils/apiUtils';

function ResultDisplay({ result, pdbContent2D, pdbContent3D, query }) {
  return (
    result && (
      <div className="bg-gray-50 rounded-md flex-1 p-4">
        {result.status === 'success' ? (
          <div>
            <div>
              <p><strong>SMILES:</strong> {result.gemini_smiles}</p>
              <p><strong>Source:</strong> {result.source}</p>
              <p><strong>Information:</strong> {result.gemini_message}</p>
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
                  onClick={() => downloadPDB('2D', result, query)}
                  className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full mt-2"
                >
                  Download 2D PDB
                </button>
                <button
                  onClick={() => downloadImage2D(query, result.image_2d)}
                  className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full mt-2"
                  disabled={!result.image_2d}
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
                    onClick={() => downloadPDB('3D', result, query)}
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
    )
  );
}

export default ResultDisplay;