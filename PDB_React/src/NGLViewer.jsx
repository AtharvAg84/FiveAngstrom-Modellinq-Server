import { useEffect, useRef, useState } from 'react';
import * as NGL from 'ngl';

function NGLViewer({ pdbContent3D, viewMode, features, onFeatureToggle }) {
  const stageRef = useRef(null);
  const viewerRef = useRef(null);
  const componentRef = useRef(null);
  const [viewerError, setViewerError] = useState('');

  // Validate PDB content
  const isValidPdb = (content) => {
    if (!content || typeof content !== 'string') {
      console.error('Invalid PDB content: empty or not a string');
      return false;
    }
    const lines = content.split('\n').filter(line => line.trim());
    const isValid = lines.some(line => line.startsWith('ATOM') || line.startsWith('HETATM'));
    if (!isValid) {
      console.error('No ATOM or HETATM lines found in PDB content');
    }
    return isValid;
  };

  // Check WebGL availability
  useEffect(() => {
    const canvas = document.createElement('canvas');
    const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
    if (!gl) {
      setViewerError('WebGL is not supported. Please enable WebGL or use a compatible browser.');
      console.error('WebGL not supported');
    } else {
      console.log('WebGL is supported');
    }
  }, []);

  // Initialize NGL Stage
  useEffect(() => {
    if (viewerRef.current && !stageRef.current) {
      try {
        stageRef.current = new NGL.Stage(viewerRef.current, { backgroundColor: '#2a2a2a' });
        stageRef.current.setSize('100%', '100%'); // Explicitly set size
        const rect = viewerRef.current.getBoundingClientRect();
        if (rect.width === 0 || rect.height === 0) {
          setViewerError('Viewer container has zero dimensions. Check CSS or parent size.');
          console.error('Zero dimensions detected', rect);
        } else {
          console.log('NGL Stage initialized', rect);
        }
      } catch (err) {
        setViewerError('Failed to initialize NGL Viewer: ' + err.message);
        console.error('NGL Stage init error:', err);
      }
    }
    return () => {
      if (stageRef.current) {
        stageRef.current.removeAllComponents();
        stageRef.current.dispose();
        stageRef.current = null;
      }
      componentRef.current = null;
      console.log('NGL Stage cleaned up');
    };
  }, []);

  // Update representations
  const updateRepresentations = (component, isProtein) => {
    if (!component) {
      console.error('No component for representations');
      return;
    }
    component.removeAllRepresentations();
    console.log('Updating representations:', features);

    if (viewMode === '3D') {
      if (features.backbone) {
        component.addRepresentation('backbone', {
          sele: isProtein ? 'protein' : 'all',
          color: 'blue',
          radius: 0.5,
          linewidth: 1.0,
          aspectRatio: 1.0,
        });
        console.log('Added backbone');
      }
      if (features.cartoon) {
        component.addRepresentation('cartoon', {
          sele: isProtein ? 'protein' : 'all',
          color: 'green',
          opacity: 0.8,
        });
        console.log('Added cartoon');
      }
      if (features.line) {
        component.addRepresentation('line', {
          sele: 'all',
          color: 'black',
          linewidth: 3,
        });
        console.log('Added line');
      }
      if (features.ballAndStick) {
        component.addRepresentation('ball+stick', {
          sele: 'all',
          color: 'element',
          radius: 0.2,
          aspectRatio: 2.0,
        });
        console.log('Added ball+stick');
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
        console.log('Added label');
      }
      if (!Object.values(features).some((v) => v)) {
        component.addRepresentation('ball+stick', {
          sele: 'all',
          color: 'element',
          radius: 0.2,
          aspectRatio: 2.0,
        });
        console.log('Added fallback ball+stick');
      }
    }
    if (stageRef.current) {
      stageRef.current.autoView();
      console.log('Auto-view applied');
    }
  };

  // Load PDB content
  useEffect(() => {
    if (!stageRef.current || viewMode !== '3D') {
      console.log('Skipping load: no stage or not 3D');
      return;
    }

    const loadStructure = async () => {
      stageRef.current.removeAllComponents();
      componentRef.current = null;
      console.log('Cleared components');

      if (!pdbContent3D) {
        setViewerError('No PDB content available. Please select a PDB ID.');
        console.log('No PDB content');
        return;
      }

      if (!isValidPdb(pdbContent3D)) {
        setViewerError('Invalid PDB format. Ensure valid PDB data.');
        console.log('Invalid PDB content:', pdbContent3D.substring(0, 100) + '...');
        return;
      }

      setViewerError('');
      console.log('Loading PDB:', pdbContent3D.substring(0, 100) + '...');
      const blob = new Blob([pdbContent3D], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      try {
        const component = await stageRef.current.loadFile(url, { ext: 'pdb' });
        componentRef.current = component;
        console.log('PDB loaded successfully');

        const isProtein = pdbContent3D.includes('HELIX') || pdbContent3D.includes('SHEET') ||
          pdbContent3D.split('\n').filter(line => line.startsWith('ATOM')).length > 50;

        updateRepresentations(component, isProtein);
      } catch (err) {
        setViewerError(`Failed to load PDB: ${err.message}`);
        console.error('Load error:', err);
      } finally {
        URL.revokeObjectURL(url);
      }
    };

    loadStructure();
  }, [pdbContent3D, viewMode, features]);

  return (
    <div className="w-full h-full flex flex-col">
      <div className="bg-gray-800/90 p-3 rounded-md border border-gray-700/50 mb-3">
        <h3 className="text-base font-semibold mb-2 text-cyan-300">Viewer Options</h3>
        <div className="grid grid-cols-3 gap-2 text-sm text-gray-200">
          {['backbone', 'cartoon', 'line', 'ballAndStick', 'label'].map((feature) => (
            <label key={feature} className="flex items-center">
              <input
                type="checkbox"
                checked={features[feature]}
                onChange={() => onFeatureToggle(feature)}
                className="mr-2 accent-cyan-600"
              />
              {feature.charAt(0).toUpperCase() + feature.slice(1).replace(/([A-Z])/g, ' $1')}
            </label>
          ))}
        </div>
      </div>
      <div className="flex-grow rounded-md border border-gray-700/50 overflow-hidden min-h-[500px]">
        <div className="w-full h-full" ref={viewerRef}>
          {viewerError && <p className="text-red-400 text-sm text-center p-4">{viewerError}</p>}
        </div>
      </div>
    </div>
  );
}

export default NGLViewer;