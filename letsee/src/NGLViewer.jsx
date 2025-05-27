import { useEffect, useRef } from 'react';
import * as NGL from 'ngl';

export default function NGLViewer({ pdbData }) {
  const viewerRef = useRef();

  useEffect(() => {
    if (!pdbData) return;

    const stage = new NGL.Stage(viewerRef.current, { backgroundColor: "white" });

    const blob = new Blob([pdbData], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);

    stage.loadFile(url, { ext: 'pdb' }).then((component) => {
      component.addRepresentation('cartoon', { colorScheme: 'chainid' });
      component.autoView();
    });

    return () => {
      stage.dispose();
      URL.revokeObjectURL(url);
    };
  }, [pdbData]);

  return <div ref={viewerRef} style={{ width: '100%', height: '400px' }} />;
}
