import React from "react";
import "./App.css";

export default function VisualizationControls({ viewerFeatures, setViewerFeatures, colorOptions, setColorOptions }) {
  const representationList = [
    // { key: "backbone", label: "Backbone" }, // COMMENTED OUT
    // { key: "cartoon", label: "Cartoon" }, // COMMENTED OUT
    { key: "line", label: "Line" },
    { key: "ballAndStick", label: "Ball+Stick" },
    { key: "label", label: "Label" },
    // { key: "surface", label: "Surface" }, // COMMENTED OUT
    { key: "licorice", label: "Licorice" },
    { key: "spacefill", label: "Spacefill" },
  ];
  // Only a few color options, plus 'element' as default
  const colorChoices = [
    { value: "element", label: "Element (default)" },
    { value: "purple", label: "Purple" },
    { value: "grey", label: "Grey" },
    { value: "red", label: "Red" },
    { value: "blue", label: "Blue" },
    { value: "green", label: "Green" },
  ];
  return (
    <div className="visualization-controls-bar">
      {representationList.map((rep) => (
        <div key={rep.key} className="vis-control-group">
          <button
            onClick={() => setViewerFeatures((prev) => ({ ...prev, [rep.key]: !prev[rep.key] }))}
            className={`visualization-button ${viewerFeatures[rep.key] ? "active" : ""}`}
          >
            {rep.label}
          </button>
          <select
            className="color-dropdown"
            value={colorOptions[rep.key] || "element"}
            onChange={e => setColorOptions((prev) => ({ ...prev, [rep.key]: e.target.value }))}
            disabled={!viewerFeatures[rep.key]}
          >
            {colorChoices.map(opt => (
              <option key={opt.value} value={opt.value}>{opt.label}</option>
            ))}
          </select>
        </div>
      ))}
    </div>
  );
}
