import { useState, useEffect, useRef } from "react";
import * as NGL from "ngl";

function App() {
  const [query, setQuery] = useState("");
  const [result, setResult] = useState(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const [pdbContent2D, setPdbContent2D] = useState("");
  const [pdbContent3D, setPdbContent3D] = useState("");
  const [image2D, setImage2D] = useState("");
  const [viewMode, setViewMode] = useState("3D");
  const [viewerError, setViewerError] = useState("");
  const [features, setFeatures] = useState({
    backbone: false,
    cartoon: false,
    line: true,
    ballAndStick: false,
    label: false,
  });
  const [colors, setColors] = useState({
    backbone: "blue",
    cartoon: "green",
    line: "black",
    ballAndStick: "element",
    label: "black",
  });
  const [form, setForm] = useState({
    sequence: "",
    samples: 1,
    forceField: "amber03.xml",
    grid: "1",
    minimize: "0",
  });
  const [invalidSequence, setInvalidSequence] = useState(false);
  const stageRef = useRef(null);
  const viewerRef = useRef(null);
  const componentRef = useRef(null);

  // Available color options
  const colorOptions = [
    "red",
    "blue",
    "green",
    "yellow",
    "purple",
    "orange",
    "cyan",
    "magenta",
    "white",
    "rainbow",
    "element",
    "residueindex",
    "chain",
    "sstruc",
  ];

  // Validate PDB content
  const isValidPdb = (content) => {
    if (!content || typeof content !== "string") {
      console.log("PDB validation failed: Content is empty or not a string");
      return false;
    }
    const lines = content.split("\n").filter((line) => line.trim());
    const hasAtoms = lines.some(
      (line) => line.startsWith("ATOM") || line.startsWith("HETATM")
    );
    console.log(
      `Validating PDB: Has ATOM/HETATM lines: ${hasAtoms}, Total lines: ${lines.length}`
    );
    return hasAtoms;
  };

  // Validate sequence
  const validateSequence = (sequence) => {
    const validAminoAcids = /^[A-Z]*$/;
    const isValid = sequence === "" || validAminoAcids.test(sequence);
    setInvalidSequence(!isValid);
    return isValid;
  };

  // Check WebGL availability
  useEffect(() => {
    const canvas = document.createElement("canvas");
    const gl =
      canvas.getContext("webgl") || canvas.getContext("experimental-webgl");
    if (!gl) {
      setViewerError(
        "WebGL is not supported in your browser. Please enable WebGL or use a compatible browser."
      );
      console.error("WebGL not supported");
    }
  }, []);

  // Initialize NGL Stage with black background
  useEffect(() => {
    const initNGL = () => {
      if (viewerRef.current && !stageRef.current) {
        console.log("Initializing NGL Stage");
        try {
          stageRef.current = new NGL.Stage(viewerRef.current, {
            backgroundColor: "black", // Set background to black
          });
          const rect = viewerRef.current.getBoundingClientRect();
          console.log("Viewer container dimensions:", {
            width: rect.width,
            height: rect.height,
          });
          if (rect.width === 0 || rect.height === 0) {
            setViewerError(
              "Viewer container has zero dimensions. Check CSS styles or parent container size."
            );
            console.error("Zero dimensions detected for viewer container");
          }
        } catch (err) {
          console.error("NGL Stage initialization error:", err);
          setViewerError("Failed to initialize NGL Viewer: " + err.message);
        }
      }
    };

    const timer = setTimeout(initNGL, 100);
    return () => {
      clearTimeout(timer);
      if (stageRef.current) {
        stageRef.current.removeAllComponents();
        stageRef.current.dispose();
        stageRef.current = null;
      }
      componentRef.current = null;
    };
  }, []);

  // Toggle feature state
  const toggleFeature = (feature) => {
    setFeatures((prev) => ({
      ...prev,
      [feature]: !prev[feature],
    }));
  };

  // Update color for a specific feature
  const updateColor = (feature, color) => {
    setColors((prev) => ({
      ...prev,
      [feature]: color,
    }));
  };

  // Update representations with enhanced color options
  const updateRepresentations = (component, isProtein) => {
    if (!component) return;
    component.removeAllRepresentations();
    console.log("Cleared all representations");

    if (viewMode === "3D") {
      if (features.backbone) {
        component.addRepresentation("backbone", {
          sele: isProtein ? "protein" : "all",
          color: colors.backbone,
          radius: 0.5,
          linewidth: 1.0,
          aspectRatio: 1.0,
        });
        console.log(
          `Applied backbone representation with color: ${colors.backbone}`
        );
      }
      if (features.cartoon) {
        component.addRepresentation("cartoon", {
          sele: isProtein ? "protein" : "all",
          color: colors.cartoon,
          opacity: 0.8,
        });
        console.log(
          `Applied cartoon representation with color: ${colors.cartoon}`
        );
      }
      if (features.line) {
        component.addRepresentation("line", {
          sele: "all",
          color: colors.line,
          linewidth: 3,
        });
        console.log(`Applied line representation with color: ${colors.line}`);
      }
      if (features.ballAndStick) {
        component.addRepresentation("ball+stick", {
          sele: "all",
          color: colors.ballAndStick,
          radius: 0.2,
          aspectRatio: 2.0,
        });
        console.log(
          `Applied ball+stick representation with color: ${colors.ballAndStick}`
        );
      }
      if (features.label) {
        component.addRepresentation("label", {
          sele: "all",
          labelType: "atomname",
          color: colors.label,
          fontSize: 0.5,
          showBackground: true,
          backgroundColor: "white",
          backgroundOpacity: 0.8,
        });
        console.log(`Applied label representation with color: ${colors.label}`);
      }

      if (!Object.values(features).some((v) => v)) {
        component.addRepresentation("spacefill", {
          sele: "all",
          color: "element",
          opacity: 0.6,
          radius: 0.9,
        });
        console.log("Applied fallback spacefill representation");
      }
    }
    if (stageRef.current) {
      stageRef.current.autoView();
      console.log("Auto-view triggered");
    }
  };

  // Load PDB content into NGL Viewer
  useEffect(() => {
    if (!stageRef.current) {
      console.log("NGL Stage not initialized");
      return;
    }

    const loadStructure = async () => {
      console.log(`Loading structure for ${viewMode}`);
      stageRef.current.removeAllComponents();
      componentRef.current = null;
      console.log("Cleared all components and componentRef");

      if (viewMode !== "3D") {
        console.log("Not in 3D mode, skipping structure load");
        return;
      }

      const pdbContent = pdbContent3D;
      console.log(
        `Attempting to load 3D PDB content:`,
        pdbContent ? `${pdbContent.slice(0, 100)}...` : "Empty"
      );

      if (!pdbContent) {
        setViewerError("No 3D PDB content available. Please submit a query.");
        return;
      }

      if (!isValidPdb(pdbContent)) {
        setViewerError(
          "Invalid 3D PDB format. Ensure the server returns valid PDB data."
        );
        console.error("Invalid 3D PDB content:", pdbContent);
        return;
      }

      setViewerError("");
      const blob = new Blob([pdbContent], { type: "text/plain" });
      const url = URL.createObjectURL(blob);
      try {
        console.log("Loading PDB file into NGL");
        const component = await stageRef.current.loadFile(url, { ext: "pdb" });
        console.log("Loaded component:", component);
        componentRef.current = component;

        const isProtein =
          pdbContent.includes("HELIX") ||
          pdbContent.includes("SHEET") ||
          pdbContent.split("\n").filter((line) => line.startsWith("ATOM"))
            .length > 50;
        console.log(`Is protein: ${isProtein}`);

        updateRepresentations(component, isProtein);
      } catch (err) {
        console.error("Error loading 3D PDB:", err);
        setViewerError(`Failed to load 3D PDB: ${err.message}`);
      } finally {
        URL.revokeObjectURL(url);
        console.log("Revoked Blob URL");
      }
    };

    loadStructure();
  }, [pdbContent3D, viewMode]);

  // Update representations when features or colors change
  useEffect(() => {
    if (componentRef.current && stageRef.current && viewMode === "3D") {
      const isProtein =
        pdbContent3D.includes("HELIX") ||
        pdbContent3D.includes("SHEET") ||
        pdbContent3D.split("\n").filter((line) => line.startsWith("ATOM"))
          .length > 50;
      updateRepresentations(componentRef.current, isProtein);
    }
  }, [features, colors, viewMode, pdbContent3D]);

  // Handle chemical query submission
  const handleQuerySubmit = async (e) => {
    e.preventDefault();
    if (!query.trim()) {
      setError("Please enter a chemical name or SMILES string.");
      return;
    }
    setLoading(true);
    setError("");
    setResult(null);
    setPdbContent2D("");
    setPdbContent3D("");
    setImage2D("");

    try {
      const response = await fetch("http://localhost:8000/get_coordinates", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: query.trim() }),
      });
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      const data = await response.json();
      console.log("API Response:", data);
      setResult(data);
      setPdbContent2D(data.pdb_content_2d || "");
      setPdbContent3D(data.pdb_content_3d || "");
      setImage2D(data.image_2d || "");
    } catch (err) {
      setError(
        "Failed to connect to the server. Please ensure the bridge and server are running."
      );
      console.error("Fetch Error:", err);
    } finally {
      setLoading(false);
    }
  };

  // Handle sidebar form submission
  const handleFormSubmit = async (e) => {
    e.preventDefault();
    if (!validateSequence(form.sequence)) {
      setError("Please enter a valid amino acid sequence.");
      return;
    }
    if (form.samples < 1) {
      setError("Number of samples must be at least 1.");
      return;
    }
    setLoading(true);
    setError("");

    try {
      const response = await fetch("http://localhost:8000/generate_structure", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(form),
      });
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      const data = await response.json();
      console.log("Sidebar Form API Response:", data);
      setPdbContent3D(data.pdb_content_3d || "");
      setResult(data);
    } catch (err) {
      setError("Failed to connect to the server for structure generation.");
      console.error("Sidebar Form Fetch Error:", err);
    } finally {
      setLoading(false);
    }
  };

  // Update form state
  const handleFormChange = (field, value) => {
    setForm((prev) => ({ ...prev, [field]: value }));
    if (field === "sequence") {
      validateSequence(value);
    }
  };

  // Download PDB
  const downloadPDB = (type) => {
    console.log(`downloadPDB called for ${type}`);
    try {
      if (!result || result.status !== "success") {
        setError("No valid PDB content available.");
        console.log("Error: Invalid result");
        return;
      }
      const content =
        type === "2D" ? result.pdb_content_2d : result.pdb_content_3d;
      if (!content) {
        setError(`No valid ${type} PDB content available.`);
        console.log(`Error: No ${type} pdb_content`);
        return;
      }
      const name = query.replace(/[^a-zA-Z0-9_-]/g, "") || "MOLECULE";
      console.log(`${type} PDB Content:`, content);
      const blob = new Blob([content], { type: "text/plain" });
      const url = URL.createObjectURL(blob);
      console.log("Blob URL:", url);
      const a = document.createElement("a");
      a.href = url;
      a.download = `${name}_${type.toLowerCase()}_output.pdb`;
      console.log("Triggering download for:", a.download);
      a.click();
      URL.revokeObjectURL(url);
      console.log("Download triggered");
    } catch (err) {
      setError(`Failed to download ${type} PDB file: ${err.message}`);
      console.error(`Download ${type} PDB Error:`, err);
    }
  };

  // Download 2D Image
  const downloadImage2D = () => {
    console.log("downloadImage2D called");
    try {
      if (!image2D) {
        setError("No 2D image available.");
        console.log("Error: No 2D image");
        return;
      }
      const name = query.replace(/[^a-zA-Z0-9_-]/g, "") || "MOLECULE";
      const a = document.createElement("a");
      a.href = `data:image/png;base64,${image2D}`;
      a.download = `${name}_2d_image.png`;
      console.log("Triggering download for:", a.download);
      a.click();
      console.log("2D Image download triggered");
    } catch (err) {
      setError(`Failed to download 2D image: ${err.message}`);
      console.error("Download 2D Image Error:", err);
    }
  };

  return (
    <div className="min-h-screen flex max-w-full w-full p-0 m-0">
      <div className="flex-1 flex flex-col md:flex-row">
        <div className="flex-1 bg-white rounded-lg shadow-lg flex flex-col p-4">
          <h1 className="text-2xl font-bold text-center">
            Chemical Coordinates
          </h1>
          <p className="text-gray-600 text-center">
            Enter a SMILES string, chemical name, or description to get 2D and
            3D coordinates.
          </p>
          {/* // Inside the App component, add this CSS for the progress bar (you
          can also put it in a separate CSS file) */}
          <style jsx>{`
            .progress-bar {
              width: 100%;
              height: 8px;
              background-color: #e5e7eb;
              border-radius: 4px;
              overflow: hidden;
            }
            .progress-bar-inner {
              height: 100%;
              background-color: #3b82f6;
              animation: progress 2s infinite;
            }
            @keyframes progress {
              0% {
                width: 0%;
              }
              50% {
                width: 80%;
              }
              100% {
                width: 0%;
              }
            }
          `}</style>
          {/* // Update the form section in the return statement */}
          <form onSubmit={handleQuerySubmit} className="space-y-2">
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
              {loading ? "Processing..." : "Ready. Steady.. Go..."}
            </button>
            {loading && (
              <div className="mt-2">
                <div className="progress-bar">
                  <div className="progress-bar-inner"></div>
                </div>
                <p className="text-blue-500 text-center mt-1">
                  Fetching chemical coordinates, please wait...
                </p>
              </div>
            )}
          </form>
          {error && <p className="text-red-500 text-center">{error}</p>}
          {result && (
            <div className="bg-gray-50 rounded-md flex-1 p-4">
              {result.status === "success" ? (
                <div>
                  <div>
                    <p>
                      <strong>SMILES:</strong> {result.gemini_smiles}
                    </p>
                    <p>
                      <strong>Source:</strong> {result.source}
                    </p>
                    <p>
                      <strong>Message:</strong> {result.gemini_message}
                    </p>
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
                          <p>
                            <strong>2D PDB Content:</strong>
                          </p>
                          <pre className="bg-gray-200 rounded-md overflow-auto max-h-60 p-2">
                            {pdbContent2D}
                          </pre>
                        </div>
                      )}
                      <button
                        onClick={() => downloadPDB("2D")}
                        className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full mt-2"
                      >
                        Download 2D PDB
                      </button>
                      <button
                        onClick={downloadImage2D}
                        className="bg-green-500 text-white p-2 rounded-md hover:bg-green-600 w-full mt-2"
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
                        <h2 className="text-lg font-semibold">
                          3D Coordinates
                        </h2>
                        <pre className="bg-gray-200 rounded-md overflow-auto max-h-60 p-2">
                          {JSON.stringify(result.coordinates_3d, null, 2)}
                        </pre>
                        <br />
                        <br />
                        {pdbContent3D && (
                          <div>
                            <p>
                              <strong>3D PDB Content:</strong>
                            </p>
                            <pre className="bg-gray-200 rounded-md overflow-auto max-h-60 p-2">
                              {pdbContent3D}
                            </pre>
                          </div>
                        )}
                        <button
                          onClick={() => downloadPDB("3D")}
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
                    <p>
                      <strong>Gemini Message:</strong> {result.gemini_message}
                    </p>
                  )}
                </div>
              )}
            </div>
          )}
        </div>
        <div className="w-80 bg-gray-100 p-4 flex flex-col space-y-4">
          <h2 className="text-xl font-bold">Structure Generation</h2>
          <form onSubmit={handleFormSubmit} className="space-y-4">
            <div className="form-field">
              <label htmlFor="seq" className="block text-sm font-medium">
                Sequence
              </label>
              <div className="relative">
                <textarea
                  id="seq"
                  value={form.sequence}
                  onChange={(e) => handleFormChange("sequence", e.target.value)}
                  className={`w-full p-2 border rounded-md ${
                    invalidSequence ? "border-red-500" : ""
                  }`}
                  placeholder="Amino Acid Sequence"
                  rows="1"
                  disabled={loading}
                  required
                />
                {invalidSequence && (
                  <div className="text-red-500 text-xs mt-1">
                    Invalid Sequence
                  </div>
                )}
              </div>
            </div>
            <div className="form-field">
              <label htmlFor="samples" className="block text-sm font-medium">
                Samples
              </label>
              <input
                type="number"
                id="samples"
                value={form.samples}
                onChange={(e) =>
                  handleFormChange("samples", parseInt(e.target.value))
                }
                className="w-full p-2 border rounded-md"
                placeholder="No. of Samples"
                min="1"
                disabled={loading}
                required
              />
            </div>
            <div className="form-field">
              <label htmlFor="force" className="block text-sm font-medium">
                Force Field
              </label>
              <select
                id="force"
                value={form.forceField}
                onChange={(e) => handleFormChange("forceField", e.target.value)}
                className="w-full p-2 border rounded-md"
                disabled={loading}
              >
                <option value="amber03.xml">AMBER03</option>
                <option value="amber10.xml">AMBER10</option>
                <option value="amber96.xml">AMBER96</option>
                <option value="amber99sb.xml">AMBER99sb</option>
                <option value="amberfb15.xml">AMBERFB15</option>
              </select>
            </div>
            <div className="form-field">
              <label htmlFor="gridSplit" className="block text-sm font-medium">
                Grid Split
              </label>
              <select
                id="gridSplit"
                value={form.grid}
                onChange={(e) => handleFormChange("grid", e.target.value)}
                className="w-full p-2 border rounded-md"
                disabled={loading}
              >
                <option value="1">1</option>
                <option value="2">2</option>
                <option value="4">4</option>
                <option value="5">Any</option>
              </select>
            </div>
            <div className="form-field">
              <label htmlFor="minimize" className="block text-sm font-medium">
                Minimize
              </label>
              <select
                id="minimize"
                value={form.minimize}
                onChange={(e) => handleFormChange("minimize", e.target.value)}
                className="w-full p-2 border rounded-md"
                disabled={loading}
              >
                <option value="0">None</option>
                <option value="5">Top 5</option>
                <option value="10">Top 10</option>
                <option value="15">Top 15</option>
                <option value="20">Top 20</option>
              </select>
            </div>
            <button
              type="submit"
              disabled={loading}
              className="w-full bg-blue-500 text-white p-2 rounded-md hover:bg-blue-600 disabled:bg-blue-300 flex items-center justify-center"
            >
              <i className="fas fa-play mr-2"></i>
              Generate
            </button>
          </form>
          <h2 className="text-xl font-bold mt-4">3D Visualization Options</h2>
          <div className="flex flex-col space-y-2">
            {[
              // { name: "Backbone", key: "backbone" },
              // { name: "Cartoon", key: "cartoon" },q
              { name: "Line", key: "line" },
              { name: "Ball+Stick", key: "ballAndStick" },
              { name: "Label", key: "label" },
            ].map((feature) => (
              <div key={feature.key} className="space-y-1">
                <button
                  onClick={() => toggleFeature(feature.key)}
                  className={`w-full p-2 rounded-md text-left ${
                    features[feature.key]
                      ? "bg-blue-500 text-white"
                      : "bg-gray-200"
                  }`}
                >
                  {feature.name}
                </button>
                {features[feature.key] && (
                  <select
                    value={colors[feature.key]}
                    onChange={(e) => updateColor(feature.key, e.target.value)}
                    className="w-full p-1 border rounded-md"
                  >
                    {colorOptions.map((color) => (
                      <option key={color} value={color}>
                        {color.charAt(0).toUpperCase() + color.slice(1)}
                      </option>
                    ))}
                  </select>
                )}
              </div>
            ))}
          </div>
        </div>
        <div className="flex-1 flex flex-col bg-white rounded-lg shadow-lg p-4">
          <div className="flex justify-center space-x-2 mb-2">
            <button
              onClick={() => setViewMode("2D")}
              className={`p-2 rounded-md ${
                viewMode === "2D" ? "bg-blue-500 text-white" : "bg-gray-200"
              } disabled:opacity-50`}
              disabled={!image2D && !pdbContent2D}
            >
              View 2D
            </button>
            <button
              onClick={() => setViewMode("3D")}
              className={`p-2 rounded-md ${
                viewMode === "3D" ? "bg-blue-500 text-white" : "bg-gray-200"
              } disabled:opacity-50`}
              disabled={!pdbContent3D}
            >
              View 3D
            </button>
          </div>
          <div
            className="flex-1 w-full bg-gray-100"
            style={{ position: "relative", width: "600px", height: "400px" }}
          >
            <div
              className="w-full h-full"
              style={{ display: viewMode === "2D" ? "none" : "block" }}
              ref={viewerRef}
            >
              {viewerError && viewMode === "3D" && (
                <p className="text-red-500 text-center p-4">{viewerError}</p>
              )}
            </div>
            <div
              className="w-full h-full flex items-center justify-center"
              style={{ display: viewMode === "2D" ? "block" : "none" }}
            >
              {image2D ? (
                <img
                  src={`data:image/png;base64,${image2D}`}
                  alt="2D Molecule"
                  className="max-w-full max-h-[400px] m-auto"
                />
              ) : (
                <p className="text-red-500 text-center p-4">
                  No 2D image available. Please submit a query.
                </p>
              )}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
export default App;
