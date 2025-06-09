import { useState, useEffect, useRef } from "react";
import axios from "axios";
import Papa from "papaparse";
import * as NGL from "ngl";
import JSZip from "jszip";
import "./App.css";

function GeminiChat() {
  const [chatQuery, setChatQuery] = useState("");
  const [chatResponse, setChatResponse] = useState(null);
  const [isChatLoading, setIsChatLoading] = useState(false);

  const handleChatSubmit = async (e) => {
    e.preventDefault();
    if (!chatQuery.trim()) {
      alert("Please enter a chat query");
      return;
    }
    setIsChatLoading(true);
    try {
      const response = await axios.post("http://localhost:8000/gemini-chat", {
        query: chatQuery.trim(),
      });
      setChatResponse(response.data);
    } catch (error) {
      setChatResponse({
        status: "error",
        error:
          error.response?.data?.detail ||
          "Failed to connect to chat server",
      });
    } finally {
      setIsChatLoading(false);
    }
  };

  return (
    <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 w-full mb-4">
      <h2 className="text-lg font-semibold mb-3 text-cyan-300">
        Protein-InfoMania
      </h2>
      <form onSubmit={handleChatSubmit}>
        <div className="mb-3">
          <label
            className="block text-gray-200 text-sm mb-1"
            htmlFor="chatQuery"
          >
            Get Information
          </label>
          <input
            type="text"
            id="chatQuery"
            value={chatQuery}
            onChange={(e) => setChatQuery(e.target.value)}
            placeholder="e.g., What is hemoglobin?"
            className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
          />
        </div>
        <button
          type="submit"
          disabled={isChatLoading}
          className={`w-full py-2 bg-cyan-600 text-gray-200 rounded text-sm font-medium ${
            isChatLoading
              ? "opacity-50 cursor-not-allowed"
              : "hover:bg-cyan-700"
          }`}
        >
          {isChatLoading ? "Sending..." : "Send"}
        </button>
      </form>
      {isChatLoading && (
        <div className="mt-3">
          <p className="text-gray-200 text-sm mb-2">Waiting for response...</p>
          <div className="w-full bg-gray-700/50 rounded-full h-2">
            <div
              className="bg-cyan-600 h-2 rounded-full animate-pulse"
              style={{ width: "50%" }}
            ></div>
          </div>
        </div>
      )}
      {chatResponse && (
        <div className="mt-3">
          <h3 className="text-base font-semibold mb-2 text-cyan-400">
            Response
          </h3>
          {chatResponse.status === "error" ? (
            <p className="text-red-400 text-sm">{chatResponse.error}</p>
          ) : (
            <p className="text-gray-200 text-sm whitespace-pre-wrap">
              {chatResponse.response}
            </p>
          )}
        </div>
      )}
    </div>
  );
}

function SearchForm({ onSearch }) {
  const [query, setQuery] = useState("");
  const [limit, setLimit] = useState(10);
  const [isLoading, setIsLoading] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!query.trim()) {
      alert("Please enter a search term");
      return;
    }
    setIsLoading(true);
    await onSearch(query.trim(), parseInt(limit));
    setIsLoading(false);
  };

  return (
    <form
      onSubmit={handleSubmit}
      className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 w-full"
    >
      <h2 className="text-lg font-semibold mb-3 text-cyan-300">PDB Search</h2>
      <div className="mb-3">
        <label className="block text-gray-200 text-sm mb-1" htmlFor="query">
          Search Term
        </label>
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
        <label className="block text-gray-200 text-sm mb-1" htmlFor="limit">
          Limit (max results)
        </label>
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
          isLoading ? "opacity-50 cursor-not-allowed" : "hover:bg-cyan-700"
        }`}
      >
        {isLoading ? "Searching..." : "Search"}
      </button>
      {isLoading && (
        <div className="mt-3">
          <p className="text-gray-200 text-sm mb-2">
            Loading content, please wait...
          </p>
          <div className="w-full bg-gray-700/50 rounded-full h-2">
            <div
              className="bg-cyan-600 h-2 rounded-full animate-pulse"
              style={{ width: "50%" }}
            ></div>
          </div>
        </div>
      )}
    </form>
  );
}

function SearchResults({ result }) {
  if (!result) return null;

  let tableData = { headers: [], rows: [] };
  if (result.csv_data) {
    try {
      const parsed = Papa.parse(result.csv_data, { header: true });
      tableData.headers = parsed.meta.fields || [];
      tableData.rows = parsed.data || [];
    } catch (error) {
      console.error("Error parsing CSV:", error);
    }
  }

  const handleDownloadCSV = () => {
    if (result.csv_data) {
      const blob = new Blob([result.csv_data], { type: "text/csv" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `${result.gemini_search_term || "search"}_output.csv`;
      a.click();
      URL.revokeObjectURL(url);
    }
  };

  return (
    <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 w-full mt-4">
      <h2 className="text-lg font-semibold mb-3 text-cyan-300">
        Search Results
      </h2>
      {result.status === "error" ? (
        <div className="text-red-400 text-sm">
          <p>
            <strong>Error:</strong> {result.error}
          </p>
          {result.gemini_message && (
            <p>
              <strong>Message:</strong> {result.gemini_message}
            </p>
          )}
        </div>
      ) : (
        <div className="text-gray-200 text-sm">
          <p className="mb-2">
            <strong>Search Term:</strong> {result.gemini_search_term}
          </p>
          <p className="mb-2">
            <strong>PDB IDs ({result.pdb_ids?.length || 0}):</strong>{" "}
            {result.pdb_ids?.join(", ") || "-"}
          </p>
          <p className="mb-2">
            <strong>
              Experimental Data ({result.experimental_data?.length || 0}{" "}
              entries):
            </strong>
          </p>
          <ul className="list-disc pl-5 mb-3 text-gray-300">
            {result.experimental_data?.map((entry, index) => (
              <li key={index} className="mb-1">
                <span className="text-cyan-400">{entry.rcsb_id}</span>: Method=
                {entry.method || "-"}, Details={entry.details || "-"}
              </li>
            )) || <li>No data available</li>}
          </ul>
          <p className="mb-2">
            <strong>Message:</strong> {result.gemini_message || "-"}
          </p>
          {tableData.headers.length > 0 && (
            <div className="mb-3">
              <h3 className="text-base font-semibold mb-2 text-cyan-400">
                Data Table
              </h3>
              <div
                className="overflow-x-auto overflow-y-auto"
                style={{ maxHeight: "12rem" }}
              >
                <table className="w-full text-left text-sm text-gray-200">
                  <thead className="bg-gray-700/50 text-cyan-400">
                    <tr>
                      <th className="px-3 py-2">S.No</th>
                      {tableData.headers.map((header, index) => (
                        <th key={index} className="px-3 py-2">
                          {header}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {tableData.rows.map((row, rowIndex) => (
                      <tr
                        key={rowIndex}
                        className="border-t border-gray-600/30"
                      >
                        <td className="px-3 py-2">{rowIndex + 1}</td>
                        {tableData.headers.map((header, colIndex) => (
                          <td key={colIndex} className="px-3 py-2">
                            {row[header] || "-"}
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

function NGLViewer({
  pdbContent,
  viewerFeatures,
  isPdbLoading,
  setViewerFeatures,
  rotationActive,
  setRotationActive,
  colorScheme,
}) {
  const stageRef = useRef(null);
  const viewerRef = useRef(null);

  useEffect(() => {
    if (!viewerRef.current) return;

    if (!stageRef.current) {
      stageRef.current = new NGL.Stage(viewerRef.current, {
        backgroundColor: "black",
        cameraType: "perspective",
      });

      const handleResize = () => {
        if (stageRef.current) {
          stageRef.current.handleResize();
        }
      };
      window.addEventListener("resize", handleResize);
      return () => window.removeEventListener("resize", handleResize);
    }

    const loadStructure = async () => {
      if (isPdbLoading || !pdbContent) return;

      try {
        stageRef.current.removeAllComponents();
        const blob = new Blob([pdbContent], { type: "text/plain" });
        const url = URL.createObjectURL(blob);
        const structure = await stageRef.current.loadFile(url, { ext: "pdb" });
        URL.revokeObjectURL(url);

        const color = colorScheme === "default" ? "element" : colorScheme;

        if (viewerFeatures.backbone) {
          structure.addRepresentation("backbone", {
            sele: ":A",
            color,
            radius: 0.5,
          });
        }
        if (viewerFeatures.cartoon) {
          structure.addRepresentation("cartoon", {
            sele: ":A",
            color,
            radius: 0.3,
          });
        }
        if (viewerFeatures.line) {
          structure.addRepresentation("line", {
            sele: "all",
            color,
            opacity: 0.7,
          });
        }
        if (viewerFeatures.ballAndStick) {
          structure.addRepresentation("ball+stick", {
            sele: "not protein and not water",
            color,
            radiusScale: 0.5,
          });
        }
        if (viewerFeatures.label) {
          structure.addRepresentation("label", {
            sele: "sidechain",
            labelType: "resname",
            color: "white",
            scale: 0.8,
            showBackground: true,
          });
        }
        if (viewerFeatures.surface) {
          structure.addRepresentation("surface", {
            sele: "all",
            color,
            opacity: 0.6,
          });
        }
        if (viewerFeatures.spacefill) {
          structure.addRepresentation("spacefill", {
            sele: "all",
            color,
            radiusScale: 0.8,
          });
        }
        if (viewerFeatures.licorice) {
          structure.addRepresentation("licorice", {
            sele: "not protein",
            color,
            radiusScale: 0.5,
          });
        }

        stageRef.current.autoView(1000);
      } catch (error) {
        console.error("Error loading PDB in NGL:", error);
      }
    };

    loadStructure();

    return () => {
      if (stageRef.current) {
        stageRef.current.removeAllComponents();
      }
    };
  }, [pdbContent, isPdbLoading, viewerFeatures, colorScheme]);

  useEffect(() => {
    let animationId;
    const animate = () => {
      if (rotationActive && stageRef.current) {
        stageRef.current.viewerControls.rotate({ x: 0.01, y: 0.01, z: 0 });
        animationId = requestAnimationFrame(animate);
      }
    };

    if (rotationActive) {
      animationId = requestAnimationFrame(animate);
    }

    return () => cancelAnimationFrame(animationId);
  }, [rotationActive]);

  return (
    <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 w-full">
      <h3 className="text-base font-semibold mb-2 text-cyan-300">
        3D Structure Viewer
      </h3>
      {isPdbLoading ? (
        <p className="text-gray-200 text-sm">Loading PDB content...</p>
      ) : pdbContent ? (
        <div>
          <div
            ref={viewerRef}
            className="w-full"
            style={{ height: "400px", position: "relative" }}
          />
        </div>
      ) : (
        <p className="text-gray-300 text-sm">
          Select a PDB ID to view 3D structure
        </p>
      )}
    </div>
  );
}

function Sidebar({
  setViewerFeatures,
  viewerFeatures,
  rotationActive,
  setRotationActive,
  colorScheme,
  setColorScheme,
}) {
  const [form, setForm] = useState({
    sequence: "",
    samples: 10,
    forceField: "amber03",
    grid: "1",
    minimize: "0",
  });
  const [invalidSequence, setInvalidSequence] = useState(false);
  const [loading, setLoading] = useState(false);
  const [processId, setProcessId] = useState("");
  const [statusMessage, setStatusMessage] = useState("");
  const [simulationProgress, setSimulationProgress] = useState(0);
  const [resultFolderLocation, setResultFolderLocation] = useState("");
  const [firstPdbContent, setFirstPdbContent] = useState("");
  const [error, setError] = useState("");

  const forceFieldMap = {
    amber03: "amber03",
    amber10: "amber10",
    amber96: "amber96",
    "amber99sb": "amber99-sb",
    amberfb15: "amber-fb15",
  };

  const validateSequence = (value) => {
    const validAminoAcids = /^[A-Za-z]*$/;
    const isValid = value === "" || validAminoAcids.test(value);
    setInvalidSequence(!isValid);
    return isValid;
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
      forcefield: forceFieldMap[form.forceField.replace(".xml", "")],
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

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!validateInputs()) {
      return;
    }

    setLoading(true);
    setError("");
    setProcessId("");
    setStatusMessage("");
    setSimulationProgress(0);
    setResultFolderLocation("");
    setFirstPdbContent("");

    try {
      const processId = await sendSimulationRequest();
      if (!processId) {
        setLoading(false);
        return;
      }

      const success = await pollStatus(processId);
      if (!success) {
        setLoading(false);
        return;
      }

      await loadFirstSample(processId);
      await downloadData(processId, form.samples);
    } catch (err) {
      setError("Error in simulation workflow: " + err.message);
      console.error("Simulation workflow error:", err);
    } finally {
      setLoading(false);
    }
  };

  const handleFormChange = (e) => {
    const { name, value } = e.target;
    setForm((prev) => ({ ...prev, [name]: value }));
    if (name === "sequence") {
      validateSequence(value);
    }
  };

  return (
    <div className="w-64 bg-gray-900/90 p-4 flex flex-col h-full border-r border-gray-700/50">
      <div className="flex justify-between items-center mb-4">
        <h2 className="text-lg font-semibold text-cyan-300">Configuration</h2>
      </div>
      <style jsx>{`
        .progress-bar {
          width: 100%;
          height: 8px;
          background-color: #4a5568;
          border-radius: 4px;
          overflow: hidden;
        }
        .progress-bar-inner {
          height: 100%;
          background-color: #00b7eb;
          width: ${simulationProgress}%;
          transition: width 0.3s ease-in-out;
        }
      `}</style>
      <form onSubmit={handleSubmit} className="flex-1 flex flex-col space-y-3">
        <div>
          <label className="block text-gray-200 text-sm mb-1" htmlFor="seq">
            Sequence
          </label>
          <textarea
            id="seq"
            name="sequence"
            value={form.sequence}
            onChange={handleFormChange}
            placeholder="Amino Acid Sequence"
            rows="2"
            className={`w-full p-2 bg-gray-900/50 text-gray-200 border ${
              invalidSequence ? "border-red-500" : "border-gray-600/50"
            } rounded text-sm`}
            disabled={loading}
            required
          />
          {invalidSequence && (
            <p className="text-red-400 text-xs mt-1">Invalid Sequence</p>
          )}
        </div>
        <div>
          <label className="block text-gray-200 text-sm mb-1" htmlFor="samples">
            Samples
          </label>
          <input
            type="number"
            id="samples"
            name="samples"
            value={form.samples}
            onChange={handleFormChange}
            placeholder="No. of Samples"
            min="1"
            max="1000"
            className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
            disabled={loading}
            required
          />
        </div>
        <div>
          <label className="block text-gray-200 text-sm mb-1" htmlFor="force">
            Force Field
          </label>
          <select
            id="force"
            name="forceField"
            value={form.forceField}
            onChange={handleFormChange}
            className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
            disabled={loading}
          >
            <option value="amber03">AMBER03</option>
            <option value="amber10">AMBER10</option>
            <option value="amber96">AMBER96</option>
            <option value="amber99sb">AMBER99sb</option>
            <option value="amberfb15">AMBERFB15</option>
          </select>
        </div>
        <div>
          <label
            className="block text-gray-200 text-sm mb-1"
            htmlFor="gridSplit"
          >
            Grid Split
          </label>
          <select
            id="gridSplit"
            name="grid"
            value={form.grid}
            onChange={handleFormChange}
            className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
            disabled={loading}
          >
            <option value="1">1</option>
            <option value="2">2</option>
            <option value="4">4</option>
            <option value="5">Any</option>
          </select>
        </div>
        <div>
          <label
            className="block text-gray-200 text-sm mb-1"
            htmlFor="minimize"
          >
            Minimize
          </label>
          <select
            id="minimize"
            name="minimize"
            value={form.minimize}
            onChange={handleFormChange}
            className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
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
          className={`w-full py-2 bg-cyan-600 text-gray-200 rounded text-sm font-medium ${
            loading ? "opacity-50 cursor-not-allowed" : "hover:bg-cyan-700"
          }`}
        >
          {loading ? (
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
        {error && (
          <p className="text-red-400 text-sm text-center mt-2">{error}</p>
        )}
        {loading && (
          <div className="mt-2">
            <div className="progress-bar">
              <div className="progress-bar-inner"></div>
            </div>
            <p className="text-cyan-400 text-sm text-center mt-1">
              {statusMessage || "Generating protein structure, please wait..."}
            </p>
          </div>
        )}
        {processId && (
          <p className="text-gray-200 text-sm text-center mt-1">
            <strong>Process ID:</strong> {processId}
          </p>
        )}
        {statusMessage && !loading && (
          <p className="text-gray-200 text-sm text-center mt-1">
            <strong>Status:</strong> {statusMessage}
          </p>
        )}
        {resultFolderLocation && (
          <p className="text-gray-200 text-sm text-center mt-1">
            <strong>Results Saved At:</strong> {resultFolderLocation}
          </p>
        )}
      </form>
      {firstPdbContent && (
        <div className="mt-4">
          <h3 className="text-base font-semibold mb-2 text-cyan-300">
            First Sample PDB Content
          </h3>
          <pre
            className="bg-gray-900/50 text-gray-200 text-xs overflow-y-auto rounded p-2"
            style={{ maxHeight: "12rem" }}
          >
            {firstPdbContent}
          </pre>
        </div>
      )}
      <div className="mt-4">
        <h3 className="text-base font-semibold mb-2 text-cyan-300">
          3D Viewer Controls
        </h3>
        <div className="flex flex-wrap gap-2 mb-3">
          {Object.keys(viewerFeatures).map((feature) => (
            <button
              key={feature}
              type="button"
              onClick={() =>
                setViewerFeatures((prev) => ({
                  ...prev,
                  [feature]: !prev[feature],
                }))
              }
              className={`py-1 px-3 text-sm rounded ${
                viewerFeatures[feature]
                  ? "bg-cyan-600 text-gray-200"
                  : "bg-gray-700/50 text-gray-300"
              } hover:bg-cyan-700`}
            >
              {feature.charAt(0).toUpperCase() + feature.slice(1)}
            </button>
          ))}
        </div>
        <div className="mb-3">
          <label
            className="block text-gray-200 text-sm mb-1"
            htmlFor="colorScheme"
          >
            Color Scheme
          </label>
          <select
            id="colorScheme"
            value={colorScheme}
            onChange={(e) => setColorScheme(e.target.value)}
            className="w-full p-2 bg-gray-900/50 text-gray-200 border border-gray-600/50 rounded text-sm"
          >
            <option value="default">Default (Element)</option>
            <option value="red">Red</option>
            <option value="blue">Blue</option>
            <option value="green">Green</option>
            <option value="yellow">Yellow</option>
          </select>
        </div>
      </div>
    </div>
  );
}

function App() {
  const [result, setResult] = useState(null);
  const [selectedPdbId, setSelectedPdbId] = useState("");
  const [pdbContent, setPdbContent] = useState("");
  const [pdbFetchError, setPdbFetchError] = useState("");
  const [isPdbLoading, setIsPdbLoading] = useState(false);
  const [viewerFeatures, setViewerFeatures] = useState({
    backbone: false,
    cartoon: false,
    line: true,
    ballAndStick: false,
    label: false,
    surface: false,
    spacefill: false,
    licorice: false,
  });
  const [rotationActive, setRotationActive] = useState(false);
  const [colorScheme, setColorScheme] = useState("default");
  const [isSidebarOpen, setIsSidebarOpen] = useState(true);

  const handleSearch = async (query, limit) => {
    try {
      const response = await axios.post("http://localhost:8000/search", {
        query,
        limit,
      });
      setResult(response.data);
      setSelectedPdbId("");
      setPdbContent("");
      setPdbFetchError("");
    } catch (error) {
      setResult({
        status: "error",
        error: error.response?.data?.detail || "Failed to connect to server",
      });
      setSelectedPdbId("");
      setPdbContent("");
      setPdbFetchError("");
    }
  };

  useEffect(() => {
    const fetchPdbContent = async () => {
      if (selectedPdbId) {
        setIsPdbLoading(true);
        try {
          const response = await axios.get(
            `http://localhost:8000/fetch-pdb/${selectedPdbId}`
          );
          setPdbContent(response.data.pdb_content);
          setPdbFetchError("");
          console.log("Fetched PDB content for:", selectedPdbId);
        } catch (error) {
          console.error("Error fetching PDB content:", error);
          setPdbContent("");
          setPdbFetchError(
            error.response?.data?.detail || "Failed to fetch PDB file"
          );
        } finally {
          setIsPdbLoading(false);
        }
      } else {
        setPdbContent("");
        setPdbFetchError("");
        setIsPdbLoading(false);
      }
    };

    fetchPdbContent();
  }, [selectedPdbId]);

  const handleDownloadPdb = () => {
    if (pdbContent && selectedPdbId) {
      const blob = new Blob([pdbContent], { type: "text/plain" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `${selectedPdbId}.pdb`;
      a.click();
      URL.revokeObjectURL(url);
    }
  };

  return (
    <div className="flex w-full min-h-screen">
      <div
        className="w-2/5 p-4 flex flex-col items-center overflow-y-auto"
        style={{ minWidth: "350px", maxWidth: "500px", zIndex: 2 }}
      >
        <GeminiChat />
        <SearchForm onSearch={handleSearch} />
        {result && <SearchResults result={result} />}
      </div>
      <div
        className={`${
          isSidebarOpen ? "w-64" : "w-0"
        } transition-all duration-300 overflow-hidden`}
        style={{ minWidth: isSidebarOpen ? "256px" : "0", zIndex: 1 }}
      >
        <Sidebar
          setViewerFeatures={setViewerFeatures}
          viewerFeatures={viewerFeatures}
          rotationActive={rotationActive}
          setRotationActive={setRotationActive}
          colorScheme={colorScheme}
          setColorScheme={setColorScheme}
        />
      </div>
      <div
        className="w-3/5 p-4 flex flex-col bg-gray-800/90"
        style={{ minWidth: "400px", zIndex: 0 }}
      >
        <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 mb-4 max-h-32 overflow-y-auto">
          <h3 className="text-base font-semibold mb-2 text-cyan-300">
            PDB IDs
          </h3>
          {result?.status === "success" && result.pdb_ids?.length > 0 ? (
            <div className="grid grid-cols-4 gap-2">
              {result.pdb_ids.map((pdbId, index) => (
                <button
                  key={pdbId}
                  onClick={() => setSelectedPdbId(pdbId)}
                  className={`p-2 rounded text-sm flex items-center ${
                    selectedPdbId === pdbId
                      ? "bg-cyan-600 text-gray-200"
                      : "bg-gray-700/50 text-gray-300 hover:bg-gray-600/50"
                  }`}
                >
                  <span className="mr-2">{index + 1}.</span>
                  <span>{pdbId}</span>
                </button>
              ))}
            </div>
          ) : (
            <p className="text-gray-300 text-sm">No PDB IDs available</p>
          )}
        </div>
        <div className="bg-gray-800/90 p-4 rounded-md border border-gray-700/50 mb-4">
          <h3 className="text-base font-semibold mb-2 text-cyan-300">
            PDB Content: {selectedPdbId || "None Selected"}
          </h3>
          {isPdbLoading ? (
            <p className="text-gray-200 text-sm">Content is being loaded...</p>
          ) : pdbFetchError ? (
            <p className="text-red-400 text-sm">{pdbFetchError}</p>
          ) : pdbContent ? (
            <>
              <pre
                className="text-gray-200 text-xs overflow-y-auto"
                style={{ maxHeight: "12rem" }}
              >
                {pdbContent}
              </pre>
              <button
                onClick={handleDownloadPdb}
                className="mt-2 py-2 px-4 bg-cyan-600 hover:bg-cyan-700 text-gray-200 rounded text-sm font-medium"
              >
                Download PDB
              </button>
            </>
          ) : (
            <p className="text-gray-300 text-sm">
              Select a PDB ID to view content
            </p>
          )}
        </div>
        <NGLViewer
          pdbContent={pdbContent}
          viewerFeatures={viewerFeatures}
          isPdbLoading={isPdbLoading}
          setViewerFeatures={setViewerFeatures}
          rotationActive={rotationActive}
          setRotationActive={setRotationActive}
          colorScheme={colorScheme}
        />
      </div>
    </div>
  );
}

export default App;