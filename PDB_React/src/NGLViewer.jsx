import { useEffect, useRef, useState } from 'react';
import * as NGL from 'ngl';

function NGLViewer({ pdbContent3D, viewMode, features }) {
  const stageRef = useRef(null);
  const viewerRef = useRef(null);
  const componentRef = useRef(null);
  const [viewerError, setViewerError] = useState('');

  // Validate PDB content
  const isValidPdb = (content) => {
    if (!content || typeof content !== 'string') {
      return false;
    }
    const lines = content.split('\n').filter(line => line.trim());
    return lines.some(line => line.startsWith('ATOM') || line.startsWith('HETATM'));
  };

  // Check WebGL availability
  useEffect(() => {
    const canvas = document.createElement('canvas');
    const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
    if (!gl) {
      setViewerError('WebGL is not supported in your browser. Please enable WebGL or use a compatible browser.');
      console.error('WebGL not supported');
    }
  }, []);

  // Initialize NGL Stage
  useEffect(() => {
    if (viewerRef.current && !stageRef.current) {
      try {
        stageRef.current = new NGL.Stage(viewerRef.current, { backgroundColor: '#2a2a2a' });
        const rect = viewerRef.current.getBoundingClientRect();
        if (rect.width === 0 || rect.height === 0) {
          setViewerError('Viewer container has zero dimensions. Check CSS styles or parent container size.');
          console.error('Zero dimensions detected for viewer container');
        }
      } catch (err) {
        setViewerError('Failed to initialize NGL Viewer: ' + err.message);
        console.error('NGL Stage initialization error:', err);
      }
    }
    return () => {
      if (stageRef.current) {
        stageRef.current.removeAllComponents();
        stageRef.current.dispose();
        stageRef.current = null;
      }
      componentRef.current = null;
    };
  }, []);

  // Update representations
  const updateRepresentations = (component, isProtein) => {
    if (!component) return;
    component.removeAllRepresentations();

    if (viewMode === '3D') {
      if (features.backbone) {
        component.addRepresentation('backbone', {
          sele: isProtein ? 'protein' : 'all',
          color: 'blue',
          radius: 0.5,
          linewidth: 1.0,
          aspectRatio: 1.0,
        });
      }
      if (features.cartoon) {
        component.addRepresentation('cartoon', {
          sele: isProtein ? 'protein' : 'all',
          color: 'green',
          opacity: 0.8,
        });
      }
      if (features.line) {
        component.addRepresentation('line', {
          sele: 'all',
          color: 'black',
          linewidth: 3,
        });
      }
      if (features.ballAndStick) {
        component.addRepresentation('ball+stick', {
          sele: 'all',
          color: 'element',
          radius: 0.2,
          aspectRatio: 2.0,
        });
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
      }
      if (!Object.values(features).some((v) => v)) {
        component.addRepresentation('spacefill', {
          sele: 'all',
          color: 'element',
          opacity: 0.6,
          radius: 0.9,
        });
      }
    }
    if (stageRef.current) {
      stageRef.current.autoView();
    }
  };

  // Load PDB content into NGL Viewer (for 3D)
  useEffect(() => {
    if (!stageRef.current || viewMode !== '3D') return;

    const loadStructure = async () => {
      stageRef.current.removeAllComponents();
      componentRef.current = null;

      if (!pdbContent3D) {
        setViewerError('No 3D PDB content available. Please submit a query.');
        return;
      }

      if (!isValidPdb(pdbContent3D)) {
        setViewerError('Invalid 3D PDB format. Ensure valid PDB data is provided.');
        console.error('Invalid 3D PDB content:', pdbContent3D);
        return;
      }

      setViewerError('');
      const blob = new Blob([pdbContent3D], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      try {
        const component = await stageRef.current.loadFile(url, { ext: 'pdb' });
        componentRef.current = component;

        const isProtein = pdbContent3D.includes('HELIX') || pdbContent3D.includes('SHEET') ||
          pdbContent3D.split('\n').filter(line => line.startsWith('ATOM')).length > 50;

        updateRepresentations(component, isProtein);
      } catch (err) {
        setViewerError(`Failed to load 3D PDB: ${err.message}`);
        console.error('Error loading 3D PDB:', err);
      } finally {
        URL.revokeObjectURL(url);
      }
    };

    loadStructure();
  }, [pdbContent3D, viewMode, features]);

  return (
    <div style={{ width: '100%', height: '100%' }} ref={viewerRef}>
      {viewerError && <p style={{ color: 'red', textAlign: 'center', padding: '1em' }}>{viewerError}</p>}
    </div>
  );
}

export default NGLViewer;