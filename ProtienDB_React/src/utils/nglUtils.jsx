import * as NGL from 'ngl';

export const isValidPdb = (content) => {
  if (!content || typeof content !== 'string') {
    console.log('PDB validation failed: Content is empty or not a string');
    return false;
  }
  const lines = content.split('\n').filter(line => line.trim());
  const hasAtoms = lines.some(line => line.startsWith('ATOM') || line.startsWith('HETATM'));
  console.log(`Validating PDB: Has ATOM/HETATM lines: ${hasAtoms}, Total lines: ${lines.length}`);
  return hasAtoms;
};

export const initializeNGL = (viewerRef, stageRef, setViewerError, componentRef) => {
  const canvas = document.createElement('canvas');
  const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
  if (!gl) {
    setViewerError('WebGL is not supported in your browser. Please enable WebGL or use a compatible browser.');
    console.error('WebGL not supported');
    return;
  }

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
  return () => clearTimeout(timer);
};

export const updateRepresentations = (component, isProtein, features, stageRef, viewMode) => {
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