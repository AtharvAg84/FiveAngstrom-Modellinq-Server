import { useState, useEffect, useRef } from "react";
import * as NGL from "ngl";
import axios from "axios";
import JSZip from "jszip";

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
  const [chemicalLoading, setChemicalLoading] = useState(false);
  const [structureLoading, setStructureLoading] = useState(false);
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
    forceField: "amber03",
    grid: "1",
    minimize: "0",
  });
  const [invalidSequence, setInvalidSequence] = useState(false);
  const stageRef = useRef(null);
  const viewerRef = useRef(null);
  const componentRef = useRef(null);
  const [processId, setProcessId] = useState("");
  const [statusMessage, setStatusMessage] = useState("");
  const [simulationProgress, setSimulationProgress] = useState(0);
  const [resultFolderLocation, setResultFolderLocation] = useState("");
  const [firstPdbContent, setFirstPdbContent] = useState("");

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

  const forceFieldMap = {
    AMBER03: "amber03",
    AMBER10: "amber10",
    AMBER96: "amber96",
    AMBER99sb: "amber99-sb",
    AMBERFB15: "amber-fb15",
  };

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

  const validateSequence = (sequence) => {
    const validAminoAcids = /^[A-Za-z]*$/;
    const isValid = sequence === "" || validAminoAcids.test(sequence);
    setInvalidSequence(!isValid);
    return isValid;
  };

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

  useEffect(() => {
    const initNGL = () => {
      if (viewerRef.current && !stageRef.current) {
        console.log("Initializing NGL Stage");
        try {
          stageRef.current = new NGL.Stage(viewerRef.current, {
            backgroundColor: "black",
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

  const toggleFeature = (feature) => {
    setFeatures((prev) => ({
      ...prev,
      [feature]: !prev[feature],
    }));
  };

  const updateColor = (feature, color) => {
    setColors((prev) => ({
      ...prev,
      [feature]: color,
    }));
  };

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

  const handleQuerySubmit = async (e) => {
    e.preventDefault();
    if (!query.trim()) {
      setError("Please enter a chemical name or SMILES string.");
      return;
    }
    setChemicalLoading(true);
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
      setResult(data);
      setPdbContent2D(data.pdb_content_2d || "");
      setPdbContent3D(data.pdb_content_3d || "");
      setImage2D(data.image_2d || "");
    } catch (err) {
      setError(
        "Failed to connect to the server. Please ensure the bridge and server are running."
      );
    } finally {
      setChemicalLoading(false);
    }
  };

  const validateInputs = () => {
    if (!validateSequence(form.sequence)) {
      setError("Please enter a valid amino acid sequence.");
      return false;
    }
    if (form.samples < 1 || form.samples > 1000) {
      setError("Number of samples must be between 1 and 1000.");
      return false;
    }
    return true;
  };

  const sendSimulationRequest = async () => {
    const payload = {
      sequence: form.sequence.toUpperCase(),
      sample: parseInt(form.samples),
      forcefield: form.forceField,
      grid: form.grid,
    };

    try {
      const response = await axios.post("http://localhost:3000/api/sample", payload);
      console.log("Simulation Request HTTP Status Code:", response.status);
      console.log("Simulation Request Raw Response:", response.data);
      const processId = response.data.toString().trim();
      if (!processId) {
        setError("Empty process ID received from server.");
        console.log("Warning: Empty process ID");
        return null;
      }
      console.log("Process ID received:", processId);
      setProcessId(processId);
      return processId;
    } catch (err) {
      setError(
        err.response?.data?.detail ||
          "Failed to connect to the server for structure generation."
      );
      console.error("Error sending simulation request:", err);
      return null;
    }
  };

  const pollStatus = async (processId, maxDuration = 300, pollInterval = 3, maxInvalidResponses = 10) => {
    const startTime = Date.now();
    let invalidResponseCount = 0;

    console.log("Starting status polling for process ID:", processId);
    setStatusMessage("Waiting for server to register the process...");
    setSimulationProgress(10);

    await new Promise((resolve) => setTimeout(resolve, 3000));

    while ((Date.now() - startTime) / 1000 < maxDuration) {
      try {
        const response = await axios.get(`http://localhost:3000/api/status?id=${processId}`, {
          headers: {
            Accept: "text/plain, application/json, */*",
            "Cache-Control": "no-cache",
            Pragma: "no-cache",
          },
          timeout: 10000,
        });

        console.log("Poll HTTP Status:", response.status);
        console.log("Poll Response:", response.data);

        if (response.status === 200) {
          const newStatus = response.data.toString().trim();
          if (newStatus && !["", "null", "undefined", "not found html"].includes(newStatus.toLowerCase())) {
            invalidResponseCount = 0;
            if (newStatus !== statusMessage) {
              setStatusMessage(newStatus);
              const statusLower = newStatus.toLowerCase();
              if (statusLower.includes("waiting") || statusLower.includes("queue")) {
                setSimulationProgress(10);
              } else if (statusLower.includes("generating") || statusLower.includes("sample")) {
                setSimulationProgress(33);
              } else if (statusLower.includes("minimizing") || statusLower.includes("minimize")) {
                setSimulationProgress(66);
              } else if (
                statusLower.includes("completed") ||
                statusLower.includes("complete") ||
                statusLower.includes("done")
              ) {
                setSimulationProgress(100);
                setStatusMessage("Simulation completed successfully!");
                return true;
              } else if (
                statusLower.includes("error") ||
                statusLower.includes("failed") ||
                statusLower.includes("fail")
              ) {
                setError(`Simulation failed: ${newStatus}`);
                return false;
              }
            }
          } else {
            invalidResponseCount++;
            console.log(`Invalid response #${invalidResponseCount}:`, newStatus);
          }
        } else if (response.status === 404) {
          invalidResponseCount++;
          console.log("Process ID not found, retrying...");
        } else {
          invalidResponseCount++;
          console.log("Unexpected HTTP status:", response.status);
        }
      } catch (err) {
        invalidResponseCount++;
        console.error("Poll error:", err);
      }

      if (invalidResponseCount >= maxInvalidResponses) {
        setError(`Too many invalid responses (${invalidResponseCount}). Aborting polling.`);
        return false;
      }

      await new Promise((resolve) => setTimeout(resolve, pollInterval * 1000));
    }

    setError(`Polling timed out after ${maxDuration} seconds. Last status: ${statusMessage}`);
    return false;
  };

  const loadFirstSample = async (processId) => {
    try {
      const response = await axios.get(
        `http://localhost:3000/api/pdb?path=Result/${processId}/sample/sample_0000.pdb`,
        {
          headers: { "Cache-Control": "no-cache", Pragma: "no-cache" },
        }
      );
      console.log("First Sample HTTP Status:", response.status);
      console.log("First Sample Response:", response.data.slice(0, 100), "...");
      setFirstPdbContent(response.data);
      return response.data;
    } catch (err) {
      setError("Failed to load first sample PDB.");
      console.error("Error fetching first sample:", err);
      return null;
    }
  };

  const downloadData = async (processId, samples) => {
    const zip = new JSZip();
    const folder = zip.folder(processId).folder("sample");

    for (let i = 0; i < samples; i++) {
      const paddedNumber = i.toString().padStart(4, "0");
      const sampleUrl = `http://localhost:3000/api/pdb?path=Result/${processId}/sample/sample_${paddedNumber}.pdb`;
      try {
        const response = await axios.get(sampleUrl);
        console.log(`Sample ${paddedNumber} HTTP Status:`, response.status);
        if (response.data) {
          folder.file(`sample_${paddedNumber}.pdb`, response.data);
        }
      } catch (err) {
        console.error(`Error fetching sample ${paddedNumber}:`, err);
      }
    }

    try {
      const sampleOutUrl = `http://localhost:3000/api/pdb?path=Result/${processId}/sample/sampled.out`;
      const response = await axios.get(sampleOutUrl);
      console.log("Sampled.out HTTP Status:", response.status);
      if (response.data) {
        folder.file("sampled.out", response.data);
      }
    } catch (err) {
      console.error("Error fetching sampled.out:", err);
    }

    try {
      const content = await zip.generateAsync({ type: "blob" });
      let savePath = "";

      if (window.showDirectoryPicker) {
        try {
          const dirHandle = await window.showDirectoryPicker();
          const fileHandle = await dirHandle.getFileHandle(`${processId}.zip`, { create: true });
          const writable = await fileHandle.createWritable();
          await writable.write(content);
          await writable.close();
          savePath = dirHandle.name;
        } catch (err) {
          console.error("Error saving with File System Access API:", err);
          setError("Failed to save results folder. Falling back to download.");
        }
      }

      if (!savePath) {
        const url = URL.createObjectURL(content);
        const a = document.createElement("a");
        a.href = url;
        a.download = `${processId}.zip`;
        a.click();
        URL.revokeObjectURL(url);
        savePath = "Browser Downloads";
      }

      setResultFolderLocation(savePath);
      return true;
    } catch (err) {
      setError("Failed to generate or save ZIP file.");
      console.error("Error generating ZIP:", err);
      return false;
    }
  };

  const handleFormSubmit = async (e) => {
    e.preventDefault();
    if (!validateInputs()) {
      return;
    }

    setStructureLoading(true);
    setError("");
    setProcessId("");
    setStatusMessage("");
    setSimulationProgress(0);
    setResultFolderLocation("");
    setFirstPdbContent("");

    try {
      const processId = await sendSimulationRequest();
      if (!processId) {
        setStructureLoading(false);
        return;
      }

      const success = await pollStatus(processId);
      if (!success) {
        setStructureLoading(false);
        return;
      }

      await loadFirstSample(processId);
      await downloadData(processId, form.samples);
    } catch (err) {
      setError("Error in simulation workflow: " + err.message);
      console.error("Simulation workflow error:", err);
    } finally {
      setStructureLoading(false);
    }
  };

  const handleFormChange = (field, value) => {
    setForm((prev) => ({ ...prev, [field]: value }));
    if (field === "sequence") {
      validateSequence(value);
    }
  };

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
              width: ${simulationProgress}%;
              transition: width 0.3s ease-in-out;
            }
          `}</style>
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
              disabled={chemicalLoading}
              className="w-full bg-blue-500 text-white p-2 rounded-md hover:bg-blue-600 disabled:bg-blue-300"
            >
              {chemicalLoading ? "Processing..." : "Ready. Steady.. Go..."}
            </button>
            {chemicalLoading && (
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
                Sequence *(Capital)
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
                  disabled={structureLoading}
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
                max="1000"
                disabled={structureLoading}
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
                disabled={structureLoading}
              >
                <option value="amber03">AMBER03</option>
                <option value="amber10">AMBER10</option>
                <option value="amber96">AMBER96</option>
                <option value="amber99-sb">AMBER99sb</option>
                <option value="amber-fb15">AMBERFB15</option>
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
                disabled={structureLoading}
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
                disabled={structureLoading}
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
              disabled={structureLoading}
              className="w-full bg-blue-500 text-white p-2 rounded-md hover:bg-blue-600 disabled:bg-blue-300 flex items-center justify-center"
            >
              {structureLoading ? (
                <>
                  <i className="fas fa-spinner fa-spin mr-2"></i>
                  Generating...
                </>
              ) : (
                <>
                  <i className="fas fa-play mr-2"></i>
                  Generate Multiple Structure
                </>
              )}
            </button>
            {structureLoading && (
              <div className="mt-2">
                <div className="progress-bar">
                  <div className="progress-bar-inner"></div>
                </div>
                <p className="text-blue-500 text-center mt-1">
                  {statusMessage || "Generating protein structure, please wait..."}
                </p>
              </div>
            )}
            {processId && (
              <p className="text-gray-700 text-center mt-1">
                <strong>Process ID:</strong> {processId}
              </p>
            )}
            {statusMessage && !structureLoading && (
              <p className="text-gray-700 text-center mt-1">
                <strong>Status:</strong> {statusMessage}
              </p>
            )}
            {resultFolderLocation && (
              <p className="text-gray-700 text-center mt-1">
                <strong>Results Saved At:</strong> {resultFolderLocation}
              </p>
            )}
          </form>
          {firstPdbContent && (
            <div className="mt-4">
              <h2 className="text-xl font-bold">First Sample PDB Content</h2>
              <pre className="bg-gray-200 rounded-md overflow-auto max-h-60 p-2">
                {firstPdbContent}
              </pre>
            </div>
          )}
          <h2 className="text-xl font-bold mt-4">3D Visualization Options</h2>
          <div className="flex flex-col space-y-2">
            {[
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