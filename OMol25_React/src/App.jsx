import { useState } from "react";
import VisualizationControls from "./VisualizationControls";
import NGLViewer from "./NGLViewer";
import "./App.css";

function App() {
  const [query, setQuery] = useState("");
  const [response, setResponse] = useState(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const [viewerFeatures, setViewerFeatures] = useState({
    backbone: false,
    cartoon: true,
    line: false,
    ballAndStick: false,
    label: false,
    surface: false,
    licorice: false,
    spacefill: false,
  });
  const [colorOptions, setColorOptions] = useState({
    backbone: "blue",
    cartoon: "purple",
    line: "grey",
    ballAndStick: "red",
    label: "black",
    surface: "yellow",
    licorice: "magenta",
    spacefill: "teal",
  });
  const [isReloading, setIsReloading] = useState(false);
  const [nglKey, setNglKey] = useState(0);

  const handleFormSubmit = async (formData) => {
    setLoading(true);
    setError("");
    setResponse(null);
    setViewerFeatures({
      backbone: false,
      cartoon: true,
      line: false,
      ballAndStick: false,
      label: false,
      surface: false,
      licorice: false,
      spacefill: false,
    });
    setIsReloading(true);
    // Wait for the NGLViewer to show reload message and cleanup
    await new Promise((res) => setTimeout(res, 1200));
    setIsReloading(false);
    setNglKey((prev) => prev + 1); // increment key to force remount
    try {
      const res = await fetch("http://localhost:8000/search", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(formData),
      });

      if (!res.ok) {
        throw new Error(`HTTP error! Status: ${res.status}`);
      }

      const data = await res.json();
      if (
        data.pdb_content &&
        !data.pdb_content.includes("ATOM") &&
        !data.pdb_content.includes("HETATM")
      ) {
        setError("Invalid PDB data: No ATOM or HETATM records found.");
        setResponse({ ...data, pdb_content: null });
      } else {
        setResponse(data);
      }
    } catch (err) {
      setError(`Failed to fetch response: ${err.message}`);
    } finally {
      setLoading(false);
    }
  };

  const handleQuerySubmit = async (e) => {
    e.preventDefault();
    if (query.length > 100) {
      setError("Query must be 100 characters or less.");
      return;
    }
    handleFormSubmit({ query });
  };

  const downloadFile = (content, filename) => {
    const blob = new Blob([content], { type: "text/plain" });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    a.click();
    window.URL.revokeObjectURL(url);
  };

  return (
    <div className="app-container themed">
      <div className="left-column">
        <div className="molecule-search-header">
          <h1>Molecule Search</h1>
          <p>
            Enter a chemical composition (e.g.,{" "}
            <span style={{ color: "#7c3aed", fontWeight: 600 }}>'C9H8O4'</span>
            ), name (e.g.,{" "}
            <span style={{ color: "#7c3aed", fontWeight: 600 }}>'aspirin'</span>
            ), or chat query (e.g.,{" "}
            <span style={{ color: "#7c3aed", fontWeight: 600 }}>'hello'</span>
            ).
          </p>
        </div>
        <form onSubmit={handleQuerySubmit} className="form">
          <input
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="Enter query"
            className="input"
            disabled={loading}
          />
          <button type="submit" className="button" disabled={loading}>
            {loading ? "Searching..." : "Search"}
          </button>
        </form>
        {error && <div className="error">{error}</div>}
        {response && (
          <div className="response">
            {response.status === "success" ? (
              <>
                {response.chat_message ? (
                  <div>
                    <h2>Chat Response</h2>
                    <p>{response.chat_message}</p>
                  </div>
                ) : (
                  <div>
                    <h2>Molecule Data</h2>
                    <p>
                      <strong>Composition:</strong>{" "}
                      {response.gemini_composition}
                    </p>
                    <p>
                      <strong>Chemical Formula:</strong>{" "}
                      {response.molecule_data?.chemical_formula}
                    </p>
                    <p>
                      <strong>Number of Atoms:</strong>{" "}
                      {response.molecule_data?.num_atoms}
                    </p>
                    <p>
                      <strong>Message:</strong> {response.gemini_message}
                    </p>
                    <div className="download-buttons">
                      <button
                        onClick={() =>
                          downloadFile(
                            response.xyz_content,
                            response.xyz_filename
                          )
                        }
                        className="download-button"
                        disabled={loading}
                      >
                        Download XYZ
                      </button>
                      <button
                        onClick={() =>
                          downloadFile(
                            response.pdb_content,
                            response.pdb_filename
                          )
                        }
                        className="download-button"
                        disabled={loading}
                      >
                        Download PDB
                      </button>
                    </div>
                    <div className="file-content">
                      <h3>XYZ File Content</h3>
                      <pre className="file-content-pre">
                        {response.xyz_content}
                      </pre>
                    </div>
                    <div className="file-content">
                      <h3>PDB File Content</h3>
                      <pre className="file-content-pre">
                        {response.pdb_content}
                      </pre>
                    </div>
                  </div>
                )}
              </>
            ) : (
              <div>
                <h2 className="error-title">Error</h2>
                <p>{response.error}</p>
                {response.gemini_message && (
                  <p>
                    <strong>Gemini Message:</strong> {response.gemini_message}
                  </p>
                )}
              </div>
            )}
          </div>
        )}
      </div>
      <div className="right-column">
        <VisualizationControls
          viewerFeatures={viewerFeatures}
          setViewerFeatures={setViewerFeatures}
          colorOptions={colorOptions}
          setColorOptions={setColorOptions}
        />
        <NGLViewer
          key={nglKey}
          pdbContent={response?.pdb_content}
          viewerFeatures={viewerFeatures}
          colorOptions={colorOptions}
          isPdbLoading={loading}
          isReloading={isReloading}
        />
      </div>
    </div>
  );
}

export default App;
