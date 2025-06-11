import { useEffect, useRef, useState } from "react";
import * as NGL from "ngl";
import "./App.css";

export default function NGLViewer({ pdbContent, viewerFeatures, colorOptions, isPdbLoading, isReloading }) {
  const stageRef = useRef(null);
  const viewerRef = useRef(null);
  const [showReloadMsg, setShowReloadMsg] = useState(false);

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
    let cleanupTimeout;
    const loadStructure = async () => {
      if (isPdbLoading || !pdbContent) return;
      if (isReloading) {
        setShowReloadMsg(true);
        await new Promise(res => setTimeout(res, 1200));
        setShowReloadMsg(false);
      }
      try {
        stageRef.current.removeAllComponents();
        const blob = new Blob([pdbContent], { type: "text/plain" });
        const url = URL.createObjectURL(blob);
        const structure = await stageRef.current.loadFile(url, { ext: "pdb" });
        URL.revokeObjectURL(url);
        // Helper to resolve color
        const resolveColor = (repKey, fallback) => {
          const c = colorOptions[repKey] || "element";
          return c === "element" ? undefined : c;
        };
        if (viewerFeatures.line) {
          structure.addRepresentation("line", {
            sele: "all",
            color: resolveColor("line", "grey"),
            opacity: 0.7,
          });
        }
        if (viewerFeatures.ballAndStick) {
          structure.addRepresentation("ball+stick", {
            sele: "not protein and not water",
            color: resolveColor("ballAndStick", "red"),
            radiusScale: 0.5,
          });
        }
        if (viewerFeatures.licorice) {
          structure.addRepresentation("licorice", {
            sele: "all",
            color: resolveColor("licorice", "purple"),
          });
        }
        if (viewerFeatures.spacefill) {
          structure.addRepresentation("spacefill", {
            sele: "all",
            color: resolveColor("spacefill", "green"),
          });
        }
        if (viewerFeatures.label) {
          structure.addRepresentation("label", {
            sele: "all",
            labelType: "atom",
            color: resolveColor("label", "black"),
            scale: 1.0,
            showBackground: true,
            backgroundColor: "white",
            borderColor: "#7c3aed",
            borderWidth: 0.5,
            zOffset: 1,
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
      if (cleanupTimeout) clearTimeout(cleanupTimeout);
    };
  }, [pdbContent, isPdbLoading, viewerFeatures, colorOptions, isReloading]);

  return (
    <div className="ngl-viewer-container">
      <h3 className="ngl-title">3D Structure Viewer</h3>
      {showReloadMsg && (
        <div className="ngl-reload-msg">Reloading 3D viewer, please wait...</div>
      )}
      {isPdbLoading ? (
        <p className="ngl-loading">Loading PDB content...</p>
      ) : pdbContent ? (
        <div
          ref={viewerRef}
          className="ngl-viewer"
          style={{ minHeight: "400px", position: "relative" }}
        />
      ) : (
        <p className="ngl-placeholder">Enter a sequence to view 3D structure</p>
      )}
    </div>
  );
}
