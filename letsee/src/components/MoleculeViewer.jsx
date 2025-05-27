import { useRef, useEffect } from 'react';
import { initializeNGL, isValidPdb, updateRepresentations } from '../utils/nglUtils';

function MoleculeViewer({
    viewMode,
    setViewMode,
    features,
    setFeatures,
    image2D,
    pdbContent3D,
    viewerError,
    setViewerError,
}) {
    const stageRef = useRef(null);
    const viewerRef = useRef(null);
    const componentRef = useRef(null);

    const toggleFeature = (feature) => {
        setFeatures((prev) => ({
            ...prev,
            [feature]: !prev[feature],
        }));
    };

    useEffect(() => {
        initializeNGL(viewerRef, stageRef, setViewerError, componentRef);
        return () => {
            if (stageRef.current) {
                stageRef.current.removeAllComponents();
                stageRef.current.dispose();
                stageRef.current = null;
            }
            componentRef.current = null;
        };
    }, []);

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
                setViewerError('Invalid 3D PDB format. Ensure the server returns valid PDB data.');
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
                updateRepresentations(component, isProtein, features, stageRef, viewMode);
            } catch (err) {
                setViewerError(`Failed to load 3D PDB: ${err.message}`);
            } finally {
                URL.revokeObjectURL(url);
            }
        };

        loadStructure();
    }, [pdbContent3D, viewMode]);

    useEffect(() => {
        if (componentRef.current && stageRef.current && viewMode === '3D') {
            const isProtein = pdbContent3D.includes('HELIX') || pdbContent3D.includes('SHEET') ||
                pdbContent3D.split('\n').filter(line => line.startsWith('ATOM')).length > 50;
            updateRepresentations(componentRef.current, isProtein, features, stageRef, viewMode);
        }
    }, [features, viewMode, pdbContent3D]);

    return (
        <div className="flex-1 flex flex-col bg-white rounded-lg shadow-lg p-4">
            <div className="flex flex-wrap justify-center space-x-2 mb-2">
                <button
                    onClick={() => setViewMode('2D')}
                    className={`p-2 rounded-md ${viewMode === '2D' ? 'bg-blue-500 text-white' : 'bg-gray-200'} disabled:opacity-50`}
                    disabled={!image2D && !pdbContent3D}
                >
                    View 2D
                </button>
                <button
                    onClick={() => setViewMode('3D')}
                    className={`p-2 rounded-md ${viewMode === '3D' ? 'bg-blue-500 text-white' : 'bg-gray-200'} disabled:opacity-50`}
                    disabled={!pdbContent3D}
                >
                    View 3D
                </button>
                {viewMode === '3D' && (
                    <div>
                        <button
                            onClick={() => toggleFeature('backbone')}
                            className={`p-2 rounded-md ${features.backbone ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
                        >
                            Backbone
                        </button>
                        <button
                            onClick={() => toggleFeature('cartoon')}
                            className={`p-2 rounded-md ${features.cartoon ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
                        >
                            Cartoon
                        </button>
                        <button
                            onClick={() => toggleFeature('line')}
                            className={`p-2 rounded-md ${features.line ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
                        >
                            Line
                        </button>
                        <button
                            onClick={() => toggleFeature('ballAndStick')}
                            className={`p-2 rounded-md ${features.ballAndStick ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
                        >
                            Ball+Stick
                        </button>
                        <button
                            onClick={() => toggleFeature('label')}
                            className={`p-2 rounded-md ${features.label ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
                        >
                            Label
                        </button>
                    </div>
                )}
            </div>
            <div
                className="flex-1 w-full bg-gray-100"
                style={{ position: 'relative', width: '600px', height: '400px' }}
            >
                <div
                    className="w-full h-full"
                    style={{ display: viewMode === '2D' ? 'none' : 'block' }}
                    ref={viewerRef}
                >
                    {viewerError && viewMode === '3D' && (
                        <p className="text-red-500 text-center p-4">{viewerError}</p>
                    )}
                </div>
                <div
                    className="w-full h-full flex items-center justify-center"
                    style={{ display: viewMode === '2D' ? 'block' : 'none' }}
                >
                    {image2D ? (
                        <img
                            src={`data:image/png;base64,${image2D}`}
                            alt="2D Molecule"
                            className="max-w-full max-h-[400px] m-auto"
                        />
                    ) : (
                        <p className="text-red-500 text-center p-4">No 2D image available. Please submit a query.</p>
                    )}
                </div>
            </div>
        </div>
    );
}

export default MoleculeViewer;