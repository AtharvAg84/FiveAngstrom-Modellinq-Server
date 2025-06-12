import { useEffect, useRef } from 'react';
import * as NGL from 'ngl';

function NGLViewer({ pdbContent3D, viewMode, features, setViewerError, viewerError }) {
  const stageRef = useRef(null);
  const viewerRef = useRef(null);
  const componentRef = useRef(null);

  // Check WebGL availability
  useEffect(() => {
    const canvas = document.createElement('canvas');
    const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
    if (!gl) {
      setViewerError('WebGL is not supported in your browser. Please enable WebGL or use a compatible browser.');
      console.error('WebGL not supported');
    }
  }, [setViewerError]);

  // Initialize NGL Stage
  useEffect(() => {
    const initNGL = () => {
      if (viewerRef.current && !stageRef.current) {
        console.log('Initializing NGL Stage');
        try {
          stageRef.current = new NGL.Stage(viewerRef.current, { backgroundColor: 'white' });
          const rect = viewerRef.current.getBoundingClientRect();
          console.log('Viewer container dimensions:', { width: rect.width, height: rect.height });
          if (rect.width === 0 || rect.height === 0) {
            setViewerError('Viewer container has zero dimensions. Check CSS styles or parent container size.');
            console.error('Zero dimensions detected for viewer container');
          }
        } catch (err) {
          console.error('NGL Stage initialization error:', err);
          setViewerError('Failed to initialize NGL Viewer: ' + err.message);
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
  }, [setViewerError]);

  // Validate PDB content
  const isValidPdb = (content) => {
    if (!content || typeof content !== 'string') {
      console.log('PDB validation failed: Content is empty or not a string');
      return false;
    }
    const lines = content.split('\n').filter(line => line.trim());
    const hasAtoms = lines.some(line => line.startsWith('ATOM') || line.startsWith('HETATM'));
    console.log(`Validating PDB: Has ATOM/HETATM lines: ${hasAtoms}, Total lines: ${lines.length}`);
    return hasAtoms;
  };

  // Update representations
  const updateRepresentations = (component, isProtein) => {
    if (!component) return;
    component.removeAllRepresentations();
    console.log('Cleared all representations');

    if (viewMode === '3D') {
      if (features.backbone) {
        component.addRepresentation('backbone', {
          sele: isProtein ? 'protein' : 'all',
          color: 'blue',
          radius: 0.5,
          linewidth: 1.0,
          aspectRatio: 1.0,
        });
        console.log('Applied backbone representation');
      }
      if (features.cartoon) {
        component.addRepresentation('cartoon', {
          sele: isProtein ? 'protein' : 'all',
          color: 'green',
          opacity: 0.8,
        });
        console.log('Applied cartoon representation');
      }
      if (features.line) {
        component.addRepresentation('line', {
          sele: 'all',
          color: 'black',
          linewidth: 3,
        });
        console.log('Applied line representation');
      }
      if (features.ballAndStick) {
        component.addRepresentation('ball+stick', {
          sele: 'all',
          color: 'element',
          radius: 0.2,
          aspectRatio: 2.0,
        });
        console.log('Applied ball+stick representation');
      }
      if (features.label) {
        component.addRepresentation('label', {
          sele: 'all',
          labelType: 'atomname',
          color: 'black',
          fontSize: 0.5,
          showBackground: true,
          backgroundColor: 'white',
          backgroundOpacity: 0.8,
        });
        console.log('Applied label representation');
      }
      if (!Object.values(features).some((v) => v)) {
        component.addRepresentation('spacefill', {
          sele: 'all',
          color: 'element',
          opacity: 0.6,
          radius: 0.9,
        });
        console.log('Applied fallback spacefill representation');
      }
    }
    if (stageRef.current) {
      stageRef.current.autoView();
      console.log('Auto-view triggered');
    }
  };

  // Load PDB content into NGL Viewer (for 3D)
  useEffect(() => {
    if (!stageRef.current) {
      console.log('NGL Stage not initialized');
      return;
    }

    const loadStructure = async () => {
      console.log(`Loading structure for ${viewMode}`);
      stageRef.current.removeAllComponents();
      componentRef.current = null;
      console.log('Cleared all components and componentRef');

      if (viewMode !== '3D') {
        console.log('Not in 3D mode, skipping structure load');
        return;
      }

      const pdbContent = pdbContent3D;
      console.log(`Attempting to load 3D PDB content:`, pdbContent ? `${pdbContent.slice(0, 100)}...` : 'Empty');

      if (!pdbContent) {
        setViewerError('No 3D PDB content available. Please submit a query.');
        return;
      }

      if (!isValidPdb(pdbContent)) {
        setViewerError('Invalid 3D PDB format. Ensure the server returns valid PDB data.');
        console.error('Invalid 3D PDB content:', pdbContent);
        return;
      }

      setViewerError('');
      const blob = new Blob([pdbContent], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      try {
        console.log('Loading PDB file into NGL');
        const component = await stageRef.current.loadFile(url, { ext: 'pdb' });
        console.log('Loaded component:', component);
        componentRef.current = component;

        const isProtein = pdbContent.includes('HELIX') || pdbContent.includes('SHEET') ||
          pdbContent.split('\n').filter(line => line.startsWith('ATOM')).length > 50;
        console.log(`Is protein: ${isProtein}`);

        updateRepresentations(component, isProtein);
      } catch (err) {
        console.error('Error loading 3D PDB:', err);
        setViewerError(`Failed to load 3D PDB: ${err.message}`);
      } finally {
        URL.revokeObjectURL(url);
        console.log('Revoked Blob URL');
      }
    };

    loadStructure();
  }, [pdbContent3D, viewMode, setViewerError]);

  // Update representations when features change
  useEffect(() => {
    if (componentRef.current && stageRef.current && viewMode === '3D') {
      const isProtein = pdbContent3D.includes('HELIX') || pdbContent3D.includes('SHEET') ||
        pdbContent3D.split('\n').filter(line => line.startsWith('ATOM')).length > 50;
      updateRepresentations(componentRef.current, isProtein);
    }
  }, [features, viewMode, pdbContent3D]);

  return (
    <div
      className="w-full h-full"
      style={{ display: viewMode === '2D' ? 'none' : 'block' , backgroundColor: 'white' }}
      ref={viewerRef}
    >
      {viewerError && viewMode === '3D' && (
        <p className="text-red-500 text-center p-4">{viewerError}</p>
      )}
    </div>
  );
}

export default NGLViewer;