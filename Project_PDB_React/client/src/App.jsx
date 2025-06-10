    import { useState, useEffect, useRef } from "react";
    import axios from "axios";
    import Papa from "papaparse";
    import * as NGL from "ngl";
    import JSZip from "jszip";

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
              error.response?.data?.detail || "Failed to connect to chat server",
          });
        } finally {
          setIsChatLoading(false);
        }
      };

      return (
        <div className="bg-white p-4 rounded-md border border-gray-200 w-full mb-4 shadow-sm">
          <h2 className="text-lg font-semibold mb-3 text-purple-700">
            Search for Protein Information
          </h2>
          <form onSubmit={handleChatSubmit}>
            <div className="mb-3">
              <label
                className="block text-gray-800 text-sm mb-1"
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
                className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
              />
            </div>
            <button
              type="submit"
              disabled={isChatLoading}
              className={`w-full py-2 bg-purple-600 text-white rounded text-sm font-medium ${
                isChatLoading
                  ? "opacity-50 cursor-not-allowed"
                  : "hover:bg-purple-700"
              }`}
            >
              {isChatLoading ? "Sending..." : "Send"}
            </button>
          </form>
          {isChatLoading && (
            <div className="mt-3">
              <p className="text-gray-800 text-sm mb-2">Waiting for response...</p>
              <div className="progress-bar">
                <div className="progress-bar-inner" style={{ width: "50%" }}></div>
              </div>
            </div>
          )}
          {chatResponse && (
            <div className="mt-3">
              <h3 className="text-base font-semibold mb-2 text-purple-600">
                Response
              </h3>
              {chatResponse.status === "error" ? (
                <p className="text-red-500 text-sm">{chatResponse.error}</p>
              ) : (
                <p className="text-gray-800 text-sm whitespace-pre-wrap">
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
          className="bg-white p-4 rounded-md border border-gray-200 w-full shadow-sm"
        >
          <h2 className="text-lg font-semibold mb-3 text-purple-700">
            Search for PDB Search
          </h2>
          <div className="mb-3">
            <label className="block text-gray-800 text-sm mb-1" htmlFor="query">
              Search Term
            </label>
            <input
              type="text"
              id="query"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              placeholder="e.g., Hemoglobin"
              className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
            />
          </div>
          <div className="mb-3">
            <label className="block text-gray-800 text-sm mb-1" htmlFor="limit">
              Limit (max results)
            </label>
            <input
              type="number"
              id="limit"
              value={limit}
              onChange={(e) => setLimit(e.target.value)}
              min="1"
              className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
            />
          </div>
          <button
            type="submit"
            disabled={isLoading}
            className={`w-full py-2 bg-purple-600 text-white rounded text-sm font-medium ${
              isLoading ? "opacity-50 cursor-not-allowed" : "hover:bg-purple-700"
            }`}
          >
            {isLoading ? "Searching..." : "Search"}
          </button>
          {isLoading && (
            <div className="mt-3">
              <p className="text-gray-800 text-sm mb-2">
                Loading content, please wait...
              </p>
              <div className="progress-bar">
                <div className="progress-bar-inner" style={{ width: "50%" }}></div>
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
        <div className="bg-white p-4 rounded-md border border-gray-200 w-full mt-4 shadow-sm">
          <h2 className="text-lg font-semibold mb-3 text-purple-700">
            Search Results
          </h2>
          {result.status === "error" ? (
            <div className="text-red-500 text-sm">
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
            <div className="text-gray-800 text-sm">
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
              <ul className="list-disc pl-5 mb-3 text-gray-600">
                {result.experimental_data?.map((entry, index) => (
                  <li key={index} className="mb-1">
                    <span className="text-purple-600">{entry.rcsb_id}</span>:
                    Method=
                    {entry.method || "-"}, Details={entry.details || "-"}
                  </li>
                )) || <li>No data available</li>}
              </ul>
              <p className="mb-2">
                <strong>Message:</strong> {result.gemini_message || "-"}
              </p>
              {tableData.headers.length > 0 && (
                <div className="mb-3">
                  <h3 className="text-base font-semibold mb-2 text-purple-600">
                    Data Table
                  </h3>
                  <div
                    className="overflow-x-auto overflow-y-auto"
                    style={{ maxHeight: "12rem" }}
                  >
                    <table className="w-full text-left text-sm text-gray-800">
                      <thead className="bg-gray-100 text-purple-600">
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
                          <tr key={rowIndex} className="border-t border-gray-200">
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
                  className="py-2 px-4 bg-purple-600 hover:bg-purple-700 text-white rounded text-sm font-medium"
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
      setColorScheme,
      resetViewer,
    }) {
      const stageRef = useRef(null);
      const viewerRef = useRef(null);

      const handleReset = () => {
        if (stageRef.current) {
          stageRef.current.removeAllComponents();
          stageRef.current.viewerControls.orient([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
          stageRef.current.autoView(0);
        }
        setViewerFeatures({
          backbone: false,
          cartoon: false,
          line: true,
          ballAndStick: false,
          label: false,
          surface: false,
          spacefill: false,
        });
        setColorScheme("default");
        setRotationActive(false);
        resetViewer(); // Clear pdbContent and related states
      };

      useEffect(() => {
        if (!viewerRef.current) return;

        if (!stageRef.current) {
          stageRef.current = new NGL.Stage(viewerRef.current, {
            backgroundColor: "white",
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
              const ballStickRep = structure.addRepresentation("ball+stick", {
                sele: "hetero and not water",
                color,
                radiusScale: 1.0,
                aspectRatio: 2.0,
              });
              console.log("Ball+Stick representation added:", ballStickRep);
              console.log(
                "Selected atoms for ball+stick:",
                ballStickRep.repr.selection
              );
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
        <div className="bg-white p-4 rounded-md border border-gray-200 w-full shadow-sm">
          <div className="flex justify-between items-center mb-2">
            <h3 className="text-base font-semibold text-purple-700">
              3D Structure Viewer
            </h3>
            <button
              onClick={handleReset}
              className="py-1 px-3 bg-red-500 hover:bg-red-600 text-white rounded text-sm font-medium"
            >
              Reset Viewer
            </button>
          </div>
          {isPdbLoading ? (
            <p className="text-gray-800 text-sm">Loading PDB content...</p>
          ) : pdbContent ? (
            <div>
              <div
                ref={viewerRef}
                className="w-full"
                style={{ height: "400px", position: "relative" }}
              />
            </div>
          ) : (
            <p className="text-gray-600 text-sm">
              Select a PDB ID or simulation sample to view 3D structure
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
      pdbIds,
      setSimulationPdbs,
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
      const [mode, setMode] = useState("manual");
      const [selectedFastaPdbId, setSelectedFastaPdbId] = useState("");
      const [fastaLoading, setFastaLoading] = useState(false);

      const forceFieldMap = {
        amber03: "amber03",
        amber10: "amber10",
        amber96: "amber96",
        amber99sb: "amber99-sb",
        amberfb15: "amber-fb15",
      };

      const validateSequence = (value) => {
        const validAminoAcids = /^[A-Za-z]*$/;
        const isValid = value === "" || validAminoAcids.test(value);
        setInvalidSequence(!isValid);
        return isValid;
      };

      const validateInputs = () => {
        if (mode === "manual" && !validateSequence(form.sequence)) {
          setError("Please enter a valid amino acid sequence.");
          return false;
        }
        if (form.samples < 1 || form.samples > 1000) {
          setError("Number of samples must be between 1 and 1000.");
          return false;
        }
        return true;
      };

      const fetchFastaSequence = async (pdbId) => {
        if (!pdbId) return;
        setFastaLoading(true);
        setError("");
        try {
          const response = await axios.get(
            `http://localhost:8000/fetch-fasta/${pdbId}`
          );
          setForm((prev) => ({ ...prev, sequence: response.data.fasta_sequence }));
          setInvalidSequence(false);
        } catch (error) {
          setError(
            error.response?.data?.detail || "Failed to fetch FASTA sequence"
          );
          setForm((prev) => ({ ...prev, sequence: "" }));
        } finally {
          setFastaLoading(false);
        }
      };

      useEffect(() => {
        if (mode === "pdb" && selectedFastaPdbId) {
          fetchFastaSequence(selectedFastaPdbId);
        } else if (mode === "manual") {
          setForm((prev) => ({ ...prev, sequence: "" }));
          setSelectedFastaPdbId("");
        }
      }, [mode, selectedFastaPdbId]);

      const sendSimulationRequest = async () => {
        const payload = {
          sequence: form.sequence.toUpperCase(),
          sample: parseInt(form.samples),
          forcefield: forceFieldMap[form.forceField.replace(".xml", "")],
          grid: form.grid,
        };

        try {
          const response = await axios.post(
            "http://localhost:3000/api/sample",
            payload
          );
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

      const pollStatus = async (
        processId,
        maxDuration = 300,
        pollInterval = 3,
        maxInvalidResponses = 10
      ) => {
        const startTime = Date.now();
        let invalidResponseCount = 0;

        console.log("Starting status polling for process ID:", processId);
        setStatusMessage("Waiting for server to register the process...");
        setSimulationProgress(10);

        await new Promise((resolve) => setTimeout(resolve, 3000));

        while ((Date.now() - startTime) / 1000 < maxDuration) {
          try {
            const response = await axios.get(
              `http://localhost:3000/api/status?id=${processId}`,
              {
                headers: {
                  Accept: "text/plain, application/json, */*",
                  "Cache-Control": "no-cache",
                  Pragma: "no-cache",
                },
                timeout: 10000,
              }
            );

            console.log("Poll HTTP Status:", response.status);
            console.log("Poll Response:", response.data);

            if (response.status === 200) {
              const newStatus = response.data.toString().trim();
              if (
                newStatus &&
                !["", "null", "undefined", "not found html"].includes(
                  newStatus.toLowerCase()
                )
              ) {
                invalidResponseCount = 0;
                if (newStatus !== statusMessage) {
                  setStatusMessage(newStatus);
                  const statusLower = newStatus.toLowerCase();
                  if (
                    statusLower.includes("waiting") ||
                    statusLower.includes("queue")
                  ) {
                    setSimulationProgress(10);
                  } else if (
                    statusLower.includes("generating") ||
                    statusLower.includes("sample")
                  ) {
                    setSimulationProgress(33);
                  } else if (
                    statusLower.includes("minimizing") ||
                    statusLower.includes("minimize")
                  ) {
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
                console.log(
                  `Invalid response #${invalidResponseCount}:`,
                  newStatus
                );
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
            setError(
              `Too many invalid responses (${invalidResponseCount}). Aborting polling.`
            );
            return false;
          }

          await new Promise((resolve) => setTimeout(resolve, pollInterval * 1000));
        }

        setError(
          `Polling timed out after ${maxDuration} seconds. Last status: ${statusMessage}`
        );
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
        const pdbFiles = [];

        for (let i = 0; i < samples; i++) {
          const paddedNumber = i.toString().padStart(4, "0");
          const sampleUrl = `http://localhost:3000/api/pdb?path=Result/${processId}/sample/sample_${paddedNumber}.pdb`;
          try {
            const response = await axios.get(sampleUrl);
            console.log(`Sample ${paddedNumber} HTTP Status:`, response.status);
            if (response.data) {
              folder.file(`sample_${paddedNumber}.pdb`, response.data);
              pdbFiles.push({
                name: `sample_${paddedNumber}.pdb`,
                content: response.data,
              });
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
          setSimulationPdbs(pdbFiles);
          const content = await zip.generateAsync({ type: "blob" });
          let savePath = "";

          if (window.showDirectoryPicker) {
            try {
              const dirHandle = await window.showDirectoryPicker();
              const fileHandle = await dirHandle.getFileHandle(`${processId}.zip`, {
                create: true,
              });
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
        setSimulationPdbs([]);

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
        if (name === "sequence" && mode === "manual") {
          validateSequence(value);
        }
      };

      return (
        <div className="w-64 bg-white p-4 flex flex-col h-full border-r border-gray-200 shadow-sm">
          <div className="flex justify-between items-center mb-4">
            <h2 className="text-lg font-semibold text-purple-700">
              Model Configuration
            </h2>
          </div>
          <form onSubmit={handleSubmit} className="flex-1 flex flex-col space-y-3">
            <div className="flex gap-2 mb-3">
              <button
                type="button"
                onClick={() => setMode("manual")}
                className={`flex-1 py-2 text-sm rounded ${
                  mode === "manual"
                    ? "bg-purple-600 text-white"
                    : "bg-gray-200 text-gray-800 hover:bg-gray-300"
                }`}
              >
                User defined Sequence
              </button>
              <button
                type="button"
                onClick={() => setMode("pdb")}
                className={`flex-1 py-2 text-sm rounded ${
                  mode === "pdb"
                    ? "bg-purple-600 text-white"
                    : "bg-gray-200 text-gray-800 hover:bg-gray-300"
                }`}
              >
                Sequence fetched from PDB
              </button>
            </div>
            {mode === "pdb" && (
              <div>
                <label className="block text-gray-800 text-sm mb-1" htmlFor="pdbId">
                  Select PDB ID
                </label>
                <select
                  id="pdbId"
                  value={selectedFastaPdbId}
                  onChange={(e) => setSelectedFastaPdbId(e.target.value)}
                  className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
                  disabled={fastaLoading || !pdbIds?.length}
                >
                  <option value="">Select a PDB ID</option>
                  {pdbIds?.map((pdbId) => (
                    <option key={pdbId} value={pdbId}>
                      {pdbId}
                    </option>
                  ))}
                </select>
                {!pdbIds?.length && (
                  <p className="text-red-500 text-xs mt-1">No PDB IDs available</p>
                )}
              </div>
            )}
            <div>
              <label className="block text-gray-800 text-sm mb-1" htmlFor="seq">
                Sequence
              </label>
              <textarea
                id="seq"
                name="sequence"
                value={form.sequence}
                onChange={handleFormChange}
                placeholder={
                  mode === "pdb"
                    ? "FASTA sequence will appear here"
                    : "Amino Acid Sequence"
                }
                rows="2"
                className={`w-full p-2 bg-gray-50 text-gray-800 border ${
                  invalidSequence ? "border-red-500" : "border-gray-300"
                } rounded text-sm focus:ring-purple-500 focus:border-purple-500`}
                disabled={loading || mode === "pdb" || fastaLoading}
                required
              />
              {invalidSequence && mode === "manual" && (
                <p className="text-red-500 text-xs mt-1">Invalid Sequence</p>
              )}
              {fastaLoading && (
                <p className="text-purple-600 text-xs mt-1">
                  Fetching FASTA sequence...
                </p>
              )}
            </div>
            <div>
              <label className="block text-gray-800 text-sm mb-1" htmlFor="samples">
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
                className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
                disabled={loading}
                required
              />
            </div>
            <div>
              <label className="block text-gray-800 text-sm mb-1" htmlFor="force">
                Force Field
              </label>
              <select
                id="force"
                name="forceField"
                value={form.forceField}
                onChange={handleFormChange}
                className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
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
                className="block text-gray-800 text-sm mb-1"
                htmlFor="gridSplit"
              >
                Grid Split
              </label>
              <select
                id="gridSplit"
                name="grid"
                value={form.grid}
                onChange={handleFormChange}
                className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
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
                className="block text-gray-800 text-sm mb-1"
                htmlFor="minimize"
              >
                Minimize
              </label>
              <select
                id="minimize"
                name="minimize"
                value={form.minimize}
                onChange={handleFormChange}
                className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
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
              className={`w-full py-2 bg-purple-600 text-white rounded text-sm font-medium ${
                loading ? "opacity-50 cursor-not-allowed" : "hover:bg-purple-700"
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
              <p className="text-red-500 text-sm text-center mt-2">{error}</p>
            )}
            {loading && (
              <div className="mt-2">
                <div className="progress-bar">
                  <div className="progress-bar-inner"></div>
                </div>
                <p className="text-purple-600 text-sm text-center mt-1">
                  {statusMessage || "Generating protein structure, please wait..."}
                </p>
              </div>
            )}
            {processId && (
              <p className="text-gray-800 text-sm text-center mt-1">
                <strong>Process ID:</strong> {processId}
              </p>
            )}
            {statusMessage && !loading && (
              <p className="text-gray-800 text-sm text-center mt-1">
                <strong>Status:</strong> {statusMessage}
              </p>
            )}
            {resultFolderLocation && (
              <p className="text-gray-800 text-sm text-center mt-1">
                <strong>Results Saved At:</strong> {resultFolderLocation}
              </p>
            )}
            {firstPdbContent && (
              <div className="mt-4">
                <h3 className="text-base font-semibold mb-2 text-purple-700">
                  First Sample PDB Content
                </h3>
                <pre
                  className="bg-gray-50 text-gray-800 text-xs overflow-y-auto rounded p-2 border border-gray-200"
                  style={{ maxHeight: "12rem" }}
                >
                  {firstPdbContent}
                </pre>
              </div>
            )}
            <div className="mt-4">
              <h3 className="text-base font-semibold mb-2 text-purple-700">
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
                        ? "bg-purple-600 text-white"
                        : "bg-gray-200 text-gray-800"
                    } hover:bg-purple-700 hover:text-white`}
                  >
                    {feature.charAt(0).toUpperCase() + feature.slice(1)}
                  </button>
                ))}
              </div>
              <div className="mb-3">
                <label
                  className="block text-gray-800 text-sm mb-1"
                  htmlFor="colorScheme"
                >
                  Color Scheme
                </label>
                <select
                  id="colorScheme"
                  value={colorScheme}
                  onChange={(e) => setColorScheme(e.target.value)}
                  className="w-full p-2 bg-gray-50 text-gray-800 border border-gray-300 rounded text-sm focus:ring-purple-500 focus:border-purple-500"
                >
                  <option value="default">Default (Element)</option>
                  <option value="red">Red</option>
                  <option value="blue">Blue</option>
                  <option value="green">Green</option>
                  <option value="yellow">Yellow</option>
                </select>
              </div>
            </div>
          </form>
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
      });
      const [rotationActive, setRotationActive] = useState(false);
      const [colorScheme, setColorScheme] = useState("default");
      const [isSidebarOpen, setIsSidebarOpen] = useState(true);
      const [simulationPdbs, setSimulationPdbs] = useState([]);
      const [currentSimulationIndex, setCurrentSimulationIndex] = useState(-1);

      const resetViewer = () => {
        setPdbContent("");
        setSelectedPdbId("");
        setCurrentSimulationIndex(-1);
        setPdbFetchError("");
      };

      // Reset viewer controls when pdbContent is cleared
      useEffect(() => {
        if (!pdbContent) {
          setViewerFeatures({
            backbone: false,
            cartoon: false,
            line: true,
            ballAndStick: false,
            label: false,
            surface: false,
            spacefill: false,
          });
          setColorScheme("default");
          setRotationActive(false);
        }
      }, [pdbContent]);

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
          setSimulationPdbs([]);
          setCurrentSimulationIndex(-1);
        } catch (error) {
          setResult({
            status: "error",
            error: error.response?.data?.detail || "Failed to connect to server",
          });
          setSelectedPdbId("");
          setPdbContent("");
          setPdbFetchError("");
          setSimulationPdbs([]);
          setCurrentSimulationIndex(-1);
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
          } else if (
            currentSimulationIndex >= 0 &&
            simulationPdbs[currentSimulationIndex]
          ) {
            setPdbContent(simulationPdbs[currentSimulationIndex].content);
            setPdbFetchError("");
            setIsPdbLoading(false);
          } else {
            setPdbContent("");
            setPdbFetchError("");
            setIsPdbLoading(false);
          }
        };

        fetchPdbContent();
      }, [selectedPdbId, currentSimulationIndex, simulationPdbs]);

      const handleDownloadPdb = () => {
        if (pdbContent && (selectedPdbId || currentSimulationIndex >= 0)) {
          const blob = new Blob([pdbContent], { type: "text/plain" });
          const url = URL.createObjectURL(blob);
          const a = document.createElement("a");
          a.href = url;
          a.download = selectedPdbId
            ? `${selectedPdbId}.pdb`
            : simulationPdbs[currentSimulationIndex]?.name || "sample.pdb";
          a.click();
          URL.revokeObjectURL(url);
        }
      };

      const handlePreviousPdb = () => {
        if (currentSimulationIndex > 0) {
          setCurrentSimulationIndex(currentSimulationIndex - 1);
          setSelectedPdbId("");
        }
      };

      const handleNextPdb = () => {
        if (currentSimulationIndex < simulationPdbs.length - 1) {
          setCurrentSimulationIndex(currentSimulationIndex + 1);
          setSelectedPdbId("");
        }
      };

      return (
        <div className="flex w-full min-h-screen bg-gray-50">
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
              pdbIds={result?.pdb_ids || []}
              setSimulationPdbs={setSimulationPdbs}
            />
          </div>
          <div
            className="w-3/5 p-4 flex flex-col bg-white shadow-sm"
            style={{ minWidth: "400px", zIndex: 0 }}
          >
            <div className="bg-white p-4 rounded-md border border-gray-200 mb-4 max-h-32 overflow-y-auto shadow-sm">
              <h3 className="text-base font-semibold mb-2 text-purple-700">
                Click to display PDB IDs
              </h3>
              {result?.status === "success" && result.pdb_ids?.length > 0 ? (
                <div className="grid grid-cols-4 gap-2">
                  {result.pdb_ids.map((pdbId, index) => (
                    <button
                      key={pdbId}
                      onClick={() => {
                        setSelectedPdbId(pdbId);
                        setCurrentSimulationIndex(-1);
                      }}
                      className={`p-2 rounded text-sm flex items-center ${
                        selectedPdbId === pdbId
                          ? "bg-purple-600 text-white"
                          : "bg-gray-200 text-gray-800 hover:bg-gray-300"
                      }`}
                    >
                      <span className="mr-2">{index + 1}.</span>
                      <span>{pdbId}</span>
                    </button>
                  ))}
                </div>
              ) : (
                <p className="text-gray-600 text-sm">No PDB IDs available</p>
              )}
            </div>
            <div className="bg-white p-4 rounded-md border border-gray-200 mb-4 shadow-sm">
              <h3 className="text-base font-semibold mb-2 text-purple-700">
                Simulation Samples
              </h3>
              {simulationPdbs.length > 0 ? (
                <div className="flex items-center justify-between mb-3">
                  <button
                    onClick={handlePreviousPdb}
                    disabled={currentSimulationIndex <= 0}
                    className={`py-2 px-4 rounded text-sm font-medium flex items-center ${
                      currentSimulationIndex <= 0
                        ? "bg-gray-300 text-gray-500 cursor-not-allowed"
                        : "bg-purple-600 text-white hover:bg-purple-700"
                    }`}
                  >
                    <i className="fas fa-arrow-left mr-2"></i>
                    Previous
                  </button>
                  <span className="text-gray-800 text-sm">
                    {currentSimulationIndex >= 0
                      ? simulationPdbs[currentSimulationIndex].name
                      : "Select a sample"}
                    {currentSimulationIndex >= 0 &&
                      ` (${currentSimulationIndex + 1}/${simulationPdbs.length})`}
                  </span>
                  <button
                    onClick={handleNextPdb}
                    disabled={currentSimulationIndex >= simulationPdbs.length - 1}
                    className={`py-2 px-4 rounded text-sm font-medium flex items-center ${
                      currentSimulationIndex >= simulationPdbs.length - 1
                        ? "bg-gray-300 text-gray-500 cursor-not-allowed"
                        : "bg-purple-600 text-white hover:bg-purple-700"
                    }`}
                  >
                    Next
                    <i className="fas fa-arrow-right ml-2"></i>
                  </button>
                </div>
              ) : (
                <p className="text-gray-600 text-sm">
                  Run a simulation to view samples
                </p>
              )}
            </div>
            <div className="bg-white p-4 rounded-md border border-gray-200 mb-4 shadow-sm">
              <h3 className="text-base font-semibold mb-2 text-purple-700">
                PDB Content:{" "}
                {selectedPdbId ||
                  (currentSimulationIndex >= 0
                    ? simulationPdbs[currentSimulationIndex]?.name
                    : "None Selected")}
              </h3>
              {isPdbLoading ? (
                <p className="text-gray-800 text-sm">Content is being loaded...</p>
              ) : pdbFetchError ? (
                <p className="text-red-500 text-sm">{pdbFetchError}</p>
              ) : pdbContent ? (
                <>
                  <pre
                    className="text-gray-800 text-xs overflow-y-auto bg-gray-50 border border-gray-200 rounded p-2"
                    style={{ maxHeight: "12rem" }}
                  >
                    {pdbContent}
                  </pre>
                  <button
                    onClick={handleDownloadPdb}
                    className="mt-2 py-2 px-4 bg-purple-600 hover:bg-purple-700 text-white rounded text-sm font-medium"
                  >
                    Download PDB
                  </button>
                </>
              ) : (
                <p className="text-gray-600 text-sm">
                  Select a PDB ID or simulation sample to view content
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
              setColorScheme={setColorScheme}
              resetViewer={resetViewer}
            />
          </div>
        </div>
      );
    }

    export default App;