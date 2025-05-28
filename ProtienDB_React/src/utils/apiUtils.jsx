export const handleSubmit = async (
    e,
    query,
    setLoading,
    setError,
    setResult,
    setPdbContent2D,
    setPdbContent3D,
    setImage2D
  ) => {
    e.preventDefault();
    if (!query.trim()) {
      setError('Please enter a chemical name or SMILES string.');
      return;
    }
    setLoading(true);
    setError('');
    setResult(null);
    setPdbContent2D('');
    setPdbContent3D('');
    setImage2D('');
  
    try {
      const response = await fetch('http://localhost:8000/get_coordinates', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query: query.trim() }),
      });
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      const data = await response.json();
      console.log('API Response:', data);
      console.log('3D PDB Content:', data.pdb_content_3d);
      console.log('2D PDB Content:', data.pdb_content_2d);
      console.log('2D Image:', data.image_2d ? 'Present' : 'Absent');
      setResult(data);
      setPdbContent2D(data.pdb_content_2d || '');
      setPdbContent3D(data.pdb_content_3d || '');
      setImage2D(data.image_2d || '');
    } catch (err) {
      setError('Failed to connect to the server. Please ensure the bridge and server are running.');
      console.error('Fetch Error:', err);
    } finally {
      setLoading(false);
    }
  };
  
  export const downloadPDB = (type, result, query) => {
    console.log(`downloadPDB called for ${type}`);
    try {
      if (!result || result.status !== 'success') {
        throw new Error('No valid PDB content available.');
      }
      const content = type === '2D' ? result.pdb_content_2d : result.pdb_content_3d;
      if (!content) {
        throw new Error(`No valid ${type} PDB content available.`);
      }
      const name = query.replace(/[^a-zA-Z0-9_-]/g, '') || 'MOLECULE';
      console.log(`${type} PDB Content:`, content);
      const blob = new Blob([content], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      console.log('Blob URL:', url);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${name}_${type.toLowerCase()}_output.pdb`;
      console.log('Triggering download for:', a.download);
      a.click();
      URL.revokeObjectURL(url);
      console.log('Download triggered');
    } catch (err) {
      throw new Error(`Failed to download ${type} PDB file: ${err.message}`);
    }
  };
  
  export const downloadImage2D = (query, image2D) => {
    console.log('downloadImage2D called');
    try {
      if (!image2D) {
        throw new Error('No 2D image available.');
      }
      const name = query.replace(/[^a-zA-Z0-9_-]/g, '') || 'MOLECULE';
      const a = document.createElement('a');
      a.href = `data:image/png;base64,${image2D}`;
      a.download = `${name}_2d_image.png`;
      console.log('Triggering download for:', a.download);
      a.click();
      console.log('2D Image download triggered');
    } catch (err) {
      throw new Error(`Failed to download 2D image: ${err.message}`);
    }
  };