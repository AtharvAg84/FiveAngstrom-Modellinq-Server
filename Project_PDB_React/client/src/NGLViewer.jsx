import { useRef, useEffect, useCallback } from "react";
import * as NGL from "ngl";

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

  const cleanupStage = useCallback(() => {
    if (stageRef.current) {
      stageRef.current.removeAllComponents();
      stageRef.current.dispose();
      stageRef.current = null;
    }
  }, []);

  const initStage = useCallback(() => {
    if (!viewerRef.current) return;

    // Clean up existing stage before creating a new one
    cleanupStage();

    // Create new stage
    stageRef.current = new NGL.Stage(viewerRef.current, {
      backgroundColor: "white",
      cameraType: "perspective",
      lightIntensity: 1.0,
    });

    console.log("New NGL Stage Initialized");
  }, [cleanupStage]);

  const loadStructure = useCallback(async () => {
    if (!pdbContent || isPdbLoading || !stageRef.current) return;

    try {
      // Initialize new stage for each structure load
      if (!stageRef.current) {
        initStage();
      }

      stageRef.current.removeAllComponents();

      const blob = new Blob([pdbContent], { type: "text/plain" });
      const url = URL.createObjectURL(blob);
      const structure = await stageRef.current.loadFile(url, { ext: "pdb" });
      URL.revokeObjectURL(url);

      const color = colorScheme === "default" ? "element" : colorScheme;

      // Add representations based on viewerFeatures
      if (viewerFeatures.backbone) {
        structure.addRepresentation("backbone", { color, radius: 0.5 });
      }
      if (viewerFeatures.cartoon) {
        structure.addRepresentation("cartoon", { color, radius: 0.3 });
      }
      if (viewerFeatures.line) {
        structure.addRepresentation("line", { color });
      }
      if (viewerFeatures.ballAndStick) {
        structure.addRepresentation("ball+stick", {
          color,
          radiusScale: 1.0,
          aspectRatio: 2.0,
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
        structure.addRepresentation("surface", { color, opacity: 0.8 });
      }
      if (viewerFeatures.spacefill) {
        structure.addRepresentation("spacefill", { color, radiusScale: 0.8 });
      }

      stageRef.current.autoView(2000);
    } catch (error) {
      console.error("Error loading PDB:", error);
    }
  }, [pdbContent, isPdbLoading, viewerFeatures, colorScheme, initStage]);

  const handleReset = useCallback(() => {
    cleanupStage();
    initStage();
    
    setViewerFeatures({
      backbone: false,
      cartoon: true,
      line: false,
      ballAndStick: false,
      label: false,
      surface: false,
      spacefill: false,
    });
    setColorScheme("default");
    setRotationActive(false);
    resetViewer();
  }, [cleanupStage, initStage, setViewerFeatures, setColorScheme, setRotationActive, resetViewer]);

  // Initialize stage when component mounts
  useEffect(() => {
    initStage();
    
    // Handle window resize
    const handleResize = () => {
      if (stageRef.current) {
        stageRef.current.handleResize();
      }
    };
    window.addEventListener("resize", handleResize);

    return () => {
      window.removeEventListener("resize", handleResize);
      cleanupStage();
    };
  }, [initStage, cleanupStage]);

  // Load structure when pdbContent changes
  useEffect(() => {
    if (!stageRef.current) { // Check if stage is initialized
      initStage(); // If not, initialize it
    }
    loadStructure();
  }, [pdbContent, loadStructure, initStage]);

  // Handle rotation animation
  useEffect(() => {
    let animationId;
    
    const animate = () => {
      if (rotationActive && stageRef.current) {
        stageRef.current.viewerControls.rotate({ x: 0.005, y: 0.005, z: 0 });
        animationId = requestAnimationFrame(animate);
      }
    };

    if (rotationActive) {
      animationId = requestAnimationFrame(animate);
    }

    return () => {
      if (animationId) {
        cancelAnimationFrame(animationId);
      }
    };
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
        <div
          ref={viewerRef}
          className="w-full border border-gray-200 rounded"
          style={{ height: "400px" }}
        />
      ) : (
        <p className="text-gray-600 text-sm">
          Select a PDB ID or simulation sample to view 3D structure
        </p>
      )}
    </div>  
  );
}

export default NGLViewer;