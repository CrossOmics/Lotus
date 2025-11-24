// Lotus Embedding Projector - Main Application Logic with Backend Integration

const API_BASE = '/api';
let vectors = null;
let metadata = null;
let currentPlot = null;
let colorMapping = {};
let sessionId = 'default';
let clusterInfo = {}; // Store cluster information

// Function to reset all state when new data is loaded
function resetAllState() {
    console.log('[RESET] Resetting all state...');
    
    // Clear cluster info
    clusterInfo = {};
    
    // Disable all pipeline buttons except preprocess
    document.getElementById('btn-preprocess').disabled = false;
    document.getElementById('btn-cluster').disabled = true;
    document.getElementById('btn-core-selection').disabled = true;
    document.getElementById('btn-visualize').disabled = true;
    document.getElementById('btn-marker-genes').disabled = true;
    
    // Clear visualization results
    const plotDiv = document.getElementById('plot');
    if (plotDiv) {
        plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #64748b;">No visualization yet. Run preprocessing and clustering first.</div>';
    }
    
    // Clear DEG results
    const markerGenesResults = document.getElementById('marker-genes-results');
    if (markerGenesResults) {
        markerGenesResults.style.display = 'none';
        const tableBody = document.getElementById('marker-genes-table-body');
        if (tableBody) {
            tableBody.innerHTML = '<tr><td colspan="4" style="padding: 16px; text-align: center; color: #64748b;">Run Marker Genes to see results</td></tr>';
        }
    }
    
    // Reset cluster keys dropdown
    const degClusterKeySelect = document.getElementById('deg-cluster-key');
    if (degClusterKeySelect) {
        degClusterKeySelect.innerHTML = '<option value="">No cluster keys available - Run clustering first</option>';
        degClusterKeySelect.disabled = true;
    }
    
    // Hide DEG step initially
    const degStep = document.getElementById('deg-step');
    if (degStep) {
        degStep.style.display = 'none';
    }
    
    // Clear coremap visualization if exists
    const coremapViz = document.getElementById('coremap-visualization');
    if (coremapViz) {
        coremapViz.innerHTML = '';
        coremapViz.style.display = 'none';
    }
    
    console.log('[RESET] All state reset complete');
}

// Function to update available cluster keys dropdown
async function updateAvailableClusterKeys() {
    try {
        const response = await fetch(`${API_BASE}/available-cluster-keys`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                session_id: sessionId
            })
        });

        const data = await response.json();

        if (!response.ok) {
            console.warn('Could not fetch available cluster keys:', data.error);
            return;
        }

        // Update the cluster key selector
        const select = document.getElementById('deg-cluster-key');
        if (!select) return;

        // Clear existing options
        select.innerHTML = '';

        // Add available cluster keys (only show existing ones)
        if (data.cluster_keys && data.cluster_keys.length > 0) {
            data.cluster_keys.forEach(clusterKeyInfo => {
                const option = document.createElement('option');
                option.value = clusterKeyInfo.key;
                const label = `${clusterKeyInfo.key} (${clusterKeyInfo.method}, ${clusterKeyInfo.n_clusters} clusters)`;
                option.textContent = label;
                option.selected = clusterKeyInfo.key === clusterInfo.cluster_key;
                select.appendChild(option);
            });
            // Enable the select and marker genes button when cluster keys are available
            select.disabled = false;
            const btnMarkerGenes = document.getElementById('btn-marker-genes');
            if (btnMarkerGenes) {
                btnMarkerGenes.disabled = false;
            }
        } else {
            // Show message if no cluster keys found
            const option = document.createElement('option');
            option.value = '';
            option.textContent = 'No cluster keys available - Run clustering first';
            option.disabled = true;
            select.appendChild(option);
            // Disable the select and marker genes button when no cluster keys
            select.disabled = true;
            const btnMarkerGenes = document.getElementById('btn-marker-genes');
            if (btnMarkerGenes) {
                btnMarkerGenes.disabled = true;
            }
        }
    } catch (error) {
        console.warn('Error updating cluster keys:', error);
    }
}

// Color palette for categorical data
const colorPalette = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d3', '#c7c7c7', '#dbdb8d', '#9edae5'
];

// Initialize
document.addEventListener('DOMContentLoaded', () => {
    checkHealth();
    setupFileInputs();
    setupControls();
    setupPipelineButtons();
    // Initialize cluster keys dropdown (will be empty until clustering is run)
    updateAvailableClusterKeys();
});

async function checkHealth() {
    try {
        const response = await fetch(`${API_BASE}/health`);
        const data = await response.json();
        if (!data.lotus_available) {
            showStatus('Warning: Lotus not available. Some features may be disabled.', 'info');
        }
    } catch (error) {
        console.error('Health check failed:', error);
    }
}

function setupFileInputs() {
    const fileInput = document.getElementById('data-file');
    const uploadArea = document.getElementById('file-upload-area');
    
    fileInput.addEventListener('change', handleDataFile);
    
    // Drag and drop
    uploadArea.addEventListener('dragover', (e) => {
        e.preventDefault();
        uploadArea.classList.add('dragover');
    });
    
    uploadArea.addEventListener('dragleave', () => {
        uploadArea.classList.remove('dragover');
    });
    
    uploadArea.addEventListener('drop', (e) => {
        e.preventDefault();
        uploadArea.classList.remove('dragover');
        if (e.dataTransfer.files.length > 0) {
            fileInput.files = e.dataTransfer.files;
            handleDataFile({ target: fileInput });
        }
    });
    
    uploadArea.addEventListener('click', () => {
        fileInput.click();
    });
}

function setupControls() {
    document.getElementById('color-by').addEventListener('change', updatePlot);
    document.getElementById('dimension').addEventListener('change', updatePlot);
    
    const pointSizeInput = document.getElementById('point-size');
    const pointSizeValue = document.getElementById('point-size-value');
    pointSizeInput.addEventListener('input', (e) => {
        pointSizeValue.textContent = e.target.value;
        updatePlot();
    });
    
    const opacityInput = document.getElementById('opacity');
    const opacityValue = document.getElementById('opacity-value');
    opacityInput.addEventListener('input', (e) => {
        opacityValue.textContent = e.target.value;
        updatePlot();
    });
    
    document.getElementById('cluster-method').addEventListener('change', () => {
        // Update UI based on method
    });
}

function setupPipelineButtons() {
    document.getElementById('btn-preprocess').addEventListener('click', runPreprocess);
    document.getElementById('btn-cluster').addEventListener('click', runClustering);
    document.getElementById('btn-visualize').addEventListener('click', runVisualization);
    document.getElementById('btn-core-selection').addEventListener('click', runCoreSelection);
    document.getElementById('btn-marker-genes').addEventListener('click', runMarkerGenes);
}

async function handleDataFile(event) {
    const file = event.target.files[0];
    if (!file) {
        console.log('[UPLOAD] No file selected');
        return;
    }

    console.log('[UPLOAD] File selected:', file.name, 'Size:', file.size, 'bytes');
    showStatus('Uploading data...', 'info');
    
    // Determine file type from extension
    const fileName = file.name.toLowerCase();
    let fileType = 'h5ad';
    if (fileName.endsWith('.h5')) {
        fileType = 'h5';
    } else if (fileName.endsWith('.csv')) {
        fileType = 'csv';
    } else if (fileName.endsWith('.tsv')) {
        fileType = 'tsv';
    }
    
    console.log('[UPLOAD] Detected file type:', fileType);
    
    try {
        const formData = new FormData();
        formData.append('file', file);
        formData.append('type', fileType);
        formData.append('session_id', sessionId);

        console.log('[UPLOAD] Sending request to:', `${API_BASE}/upload`);
        const response = await fetch(`${API_BASE}/upload`, {
            method: 'POST',
            body: formData
        });

        console.log('[UPLOAD] Response status:', response.status, response.statusText);

        let data;
        try {
            data = await response.json();
            console.log('[UPLOAD] Response data:', data);
        } catch (e) {
            const text = await response.text();
            console.error('[UPLOAD] Failed to parse JSON:', text);
            throw new Error(`Server error: ${response.status} ${response.statusText}. Response: ${text.substring(0, 200)}`);
        }

        if (!response.ok) {
            console.error('[UPLOAD] Error response:', data);
            throw new Error(data.error || `Upload failed: ${response.status}`);
        }

        console.log('[UPLOAD] Success:', data);
        showStatus(data.message || `✓ Loaded ${data.shape[0]} cells, ${data.shape[1]} genes`, 'success');
        
        // Reset all state when new data is loaded
        resetAllState();
        
        // Enable pipeline buttons
        document.getElementById('btn-preprocess').disabled = false;
        document.getElementById('pipeline-section').style.display = 'block';
        
        // Update metadata select
        if (data.obs_columns) {
            updateMetadataSelect(data.obs_columns);
        }
        
        // Update available cluster keys (in case uploaded data already has clustering results)
        await updateAvailableClusterKeys();
        
    } catch (error) {
        console.error('[UPLOAD] Exception:', error);
        showStatus(`Error: ${error.message}`, 'error');
    }
}

async function runPreprocess() {
    showStatus('Running preprocessing... This may take a while.', 'info');
    const btn = document.getElementById('btn-preprocess');
    const originalText = btn.textContent;
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner"></span> Processing...';

    try {
        const n_pcs = parseInt(document.getElementById('n-pcs').value) || 20;
        const n_neighbors = parseInt(document.getElementById('n-neighbors').value) || 15;

        const response = await fetch(`${API_BASE}/preprocess`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                session_id: sessionId,
                n_pcs: n_pcs,
                n_neighbors: n_neighbors
            })
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Preprocessing failed');
        }

        showStatus('✓ Preprocessing complete!', 'success');
        // Enable core selection and cluster buttons (core selection is independent)
        document.getElementById('btn-core-selection').disabled = false;
        document.getElementById('btn-cluster').disabled = false;
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Preprocess error:', error);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

async function runClustering() {
    showStatus('Running clustering...', 'info');
    const btn = document.getElementById('btn-cluster');
    const originalText = btn.textContent;
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner"></span> Clustering...';

    try {
        const method = document.getElementById('cluster-method').value;
        const resolution = parseFloat(document.getElementById('cluster-resolution').value) || 1.2;
        const use_rep = document.getElementById('use-rep').value || null;

        const response = await fetch(`${API_BASE}/cluster`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                session_id: sessionId,
                method: method,
                resolution: resolution,
                use_rep: use_rep || null
            })
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Clustering failed');
        }

        showStatus(`✓ Clustering complete! Found ${data.n_clusters} clusters using ${data.cluster_key}.`, 'success');
        document.getElementById('btn-visualize').disabled = false;
        document.getElementById('btn-marker-genes').disabled = false;
        document.getElementById('deg-step').style.display = 'block';
        
        // Store cluster info
        clusterInfo = {
            cluster_key: data.cluster_key,
            clusters: data.clusters,
            cluster_counts: data.cluster_counts,
            has_model: data.has_model || false
        };
        
        // Update available cluster keys for DEG analysis (only existing ones)
        await updateAvailableClusterKeys();
        
        // Update color-by select to include cluster key
        updateMetadataSelect([data.cluster_key]);
        document.getElementById('color-by').value = data.cluster_key;
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Clustering error:', error);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

async function runVisualization() {
    showStatus('Computing visualization...', 'info');
    const btn = document.getElementById('btn-visualize');
    const originalText = btn.textContent;
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner"></span> Computing UMAP...';

    try {
        const cluster_key = clusterInfo.cluster_key || 'cplearn';
        const n_components = document.getElementById('dimension').value === '3d' ? 3 : 2;
        const min_dist = parseFloat(document.getElementById('umap-min-dist').value) || 0.5;
        const spread = parseFloat(document.getElementById('umap-spread').value) || 1.0;

        const response = await fetch(`${API_BASE}/visualize`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                session_id: sessionId,
                cluster_key: cluster_key,
                n_components: n_components,
                min_dist: min_dist,
                spread: spread
            })
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Visualization failed');
        }

        // Store data
        vectors = data.coordinates;
        metadata = data.metadata;
        
        showStatus('✓ Visualization ready!', 'success');
        updatePlot();
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Visualization error:', error);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

function updateMetadataSelect(columns = []) {
    const select = document.getElementById('color-by');
    const currentValue = select.value;
    
    // Keep existing options, add new ones
    const existingOptions = Array.from(select.options).map(opt => opt.value);
    columns.forEach(col => {
        if (!existingOptions.includes(col)) {
            const option = document.createElement('option');
            option.value = col;
            option.textContent = col;
            select.appendChild(option);
        }
    });
}

function updatePlot() {
    if (!vectors || vectors.length === 0) {
        return;
    }

    const dimension = document.getElementById('dimension').value;
    const colorBy = document.getElementById('color-by').value;
    const pointSize = parseFloat(document.getElementById('point-size').value);
    const opacity = parseFloat(document.getElementById('opacity').value);

    // Show visualization controls
    document.getElementById('visualization-section').style.display = 'block';
    document.getElementById('search-section').style.display = 'block';
    document.getElementById('stats-section').style.display = 'block';

    // Prepare data
    let x = vectors.map(v => v[0]);
    let y = vectors.map(v => v[1]);
    let z = dimension === '3d' && vectors[0].length >= 3 ? vectors.map(v => v[2]) : null;

    // Prepare colors
    let colors = null;
    if (colorBy !== 'none' && metadata && metadata[colorBy]) {
        colors = metadata[colorBy];
        colorMapping = createColorMapping(colors);
    }

    // Create trace
    const trace = {
        x: x,
        y: y,
        mode: 'markers',
        type: dimension === '3d' ? 'scatter3d' : 'scatter',
        marker: {
            size: pointSize,
            opacity: opacity,
            color: colors ? colors.map(c => getColorForValue(c)) : '#1f77b4',
            colorscale: colorBy === 'none' ? 'Viridis' : undefined,
            showscale: colorBy === 'none',
            line: {
                width: 0.5,
                color: 'rgba(0,0,0,0.1)'
            }
        },
        text: metadata ? vectors.map((v, i) => {
            const info = Object.entries(metadata).map(([k, arr]) => `${k}: ${arr[i]}`).join('<br>');
            return `Point ${i}<br>${info}`;
        }) : undefined,
        hovertemplate: '<b>%{text}</b><br>X: %{x:.2f}<br>Y: %{y:.2f}' + (z ? '<br>Z: %{z:.2f}' : '') + '<extra></extra>',
        customdata: vectors.map((v, i) => i) // Store index for click events
    };

    if (z) {
        trace.z = z;
    }

    const layout = {
        title: {
            text: `Embedding Visualization (${vectors.length} points)`,
            font: { size: 16 }
        },
        margin: { l: 50, r: 50, t: 50, b: 50 },
        hovermode: 'closest',
        scene: dimension === '3d' ? {
            xaxis: { title: 'Dimension 1' },
            yaxis: { title: 'Dimension 2' },
            zaxis: { title: 'Dimension 3' },
            camera: {
                eye: { x: 1.5, y: 1.5, z: 1.5 }
            }
        } : {
            xaxis: { title: 'Dimension 1' },
            yaxis: { title: 'Dimension 2' }
        }
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['lasso2d', 'select2d']
    };

    Plotly.newPlot('plot', [trace], layout, config).then(() => {
        // Add click event listener
        const plotDiv = document.getElementById('plot');
        plotDiv.on('plotly_click', handlePlotClick);
    });

    // Update stats
    updateStats();
}

function handlePlotClick(data) {
    if (!data || !data.points || data.points.length === 0) return;
    
    const pointIndex = data.points[0].customdata;
    if (pointIndex === undefined) return;
    
    // Get cluster ID from metadata
    if (!metadata || !clusterInfo.cluster_key) return;
    
    const clusterId = metadata[clusterInfo.cluster_key][pointIndex];
    
    // Fetch cluster information
    fetchClusterInfo(clusterId);
}

async function fetchClusterInfo(clusterId) {
    try {
        const response = await fetch(`${API_BASE}/cluster-info`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                session_id: sessionId,
                cluster_id: clusterId,
                cluster_key: clusterInfo.cluster_key
            })
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Failed to fetch cluster info');
        }

        // Display cluster information
        displayClusterInfo(data);
        
    } catch (error) {
        console.error('Error fetching cluster info:', error);
        showStatus(`Error: ${error.message}`, 'error');
    }
}

function displayClusterInfo(info) {
    const infoDiv = document.getElementById('cluster-info');
    if (!infoDiv) {
        // Create info div if it doesn't exist
        const newDiv = document.createElement('div');
        newDiv.id = 'cluster-info';
        newDiv.className = 'cluster-info-panel';
        document.querySelector('.sidebar').appendChild(newDiv);
    }
    
    let html = `<div class="section">
        <div class="section-title">Cluster ${info.cluster_id}</div>
        <div class="info">
            <strong>Cells:</strong> ${info.stats.n_cells} (${info.stats.percentage.toFixed(1)}%)<br>`;
    
    if (info.cell_types && Object.keys(info.cell_types).length > 0) {
        html += `<br><strong>Cell Types:</strong><br>`;
        Object.entries(info.cell_types).forEach(([col, types]) => {
            html += `<br><em>${col}:</em><br>`;
            Object.entries(types).forEach(([type, count]) => {
                html += `  ${type}: ${count}<br>`;
            });
        });
    }
    
    if (info.top_genes && info.top_genes.length > 0) {
        html += `<br><strong>Top Expressed Genes:</strong><br>`;
        info.top_genes.forEach(gene => {
            html += `  ${gene}<br>`;
        });
    }
    
    html += `</div></div>`;
    
    document.getElementById('cluster-info').innerHTML = html;
    document.getElementById('cluster-info').style.display = 'block';
}

function createColorMapping(values) {
    const uniqueValues = [...new Set(values.filter(v => v !== '' && v !== null && v !== undefined))];
    const mapping = {};
    
    uniqueValues.forEach((value, index) => {
        if (typeof value === 'number' || !isNaN(value)) {
            mapping[value] = value;
        } else {
            mapping[value] = colorPalette[index % colorPalette.length];
        }
    });
    
    return mapping;
}

function getColorForValue(value) {
    if (value === '' || value === null || value === undefined) {
        return '#cccccc';
    }
    
    if (typeof value === 'number' || !isNaN(value)) {
        return value;
    }
    
    return colorMapping[value] || '#cccccc';
}

function updateStats() {
    const statsDiv = document.getElementById('stats');
    if (!vectors) return;

    const stats = {
        'Total Points': vectors.length,
        'Dimensions': vectors[0].length,
        'X Range': `[${Math.min(...vectors.map(v => v[0])).toFixed(2)}, ${Math.max(...vectors.map(v => v[0])).toFixed(2)}]`,
        'Y Range': `[${Math.min(...vectors.map(v => v[1])).toFixed(2)}, ${Math.max(...vectors.map(v => v[1])).toFixed(2)}]`
    };

    if (metadata) {
        stats['Metadata Columns'] = Object.keys(metadata).length;
    }

    if (clusterInfo.clusters) {
        stats['Clusters'] = clusterInfo.clusters.length;
    }

    statsDiv.innerHTML = Object.entries(stats)
        .map(([key, value]) => `<strong>${key}:</strong> ${value}`)
        .join('<br>');
}

async function runCoreSelection() {
    showStatus('Running core selection...', 'info');
    const btn = document.getElementById('btn-core-selection');
    const originalText = btn.textContent;
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner"></span> Computing...';

    try {
        const use_rep = document.getElementById('core-use-rep')?.value || null;
        const cluster_resolution = parseFloat(document.getElementById('core-cluster-resolution')?.value) || 1.2;
        const fast_view = document.getElementById('core-fast-view')?.checked !== false; // default true
        
        // Ground truth labels from JSON file only
        const truth_json_file = document.getElementById('core-truth-json-file')?.files[0];
        
        // Prepare request - use FormData if JSON file is provided, otherwise JSON
        let requestBody;
        let headers;
        let body;
        
        if (truth_json_file) {
            // Use FormData for file upload
            const formData = new FormData();
            formData.append('session_id', sessionId);
            formData.append('use_rep', use_rep || '');
            formData.append('key_added', 'X_cplearn_coremap');
            formData.append('cluster_resolution', cluster_resolution.toString());
            formData.append('fast_view', fast_view.toString());
            formData.append('truth_json_file', truth_json_file);
            
            headers = {}; // FormData sets Content-Type automatically
            body = formData;
        } else {
            // Use JSON for regular request
            requestBody = {
                session_id: sessionId,
                use_rep: use_rep,
                key_added: 'X_cplearn_coremap',
                cluster_resolution: cluster_resolution,
                fast_view: fast_view
            };
            
            headers = { 'Content-Type': 'application/json' };
            body = JSON.stringify(requestBody);
        }

        const response = await fetch(`${API_BASE}/core-selection`, {
            method: 'POST',
            headers: headers,
            body: body
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Core selection failed');
        }

        let statusMessage = `✓ Core selection complete! Embedding stored in ${data.key_added}`;
        if (data.assigned_points !== undefined) {
            statusMessage += ` (${data.assigned_points}/${data.total_points} points assigned`;
            if (data.core_cells !== undefined) {
                statusMessage += `, ${data.core_cells} core cells`;
            }
            statusMessage += ')';
        }
        
        showStatus(statusMessage, 'success');
        console.log('Core selection result:', data);
        
        // Handle visualization result - load and display directly
        if (data.visualization) {
            if (data.visualization.success && data.visualization.file_url) {
                // Load visualization HTML and display in plot container
                try {
                    const vizUrl = data.visualization.file_url;
                    console.log('[CORE] Loading visualization from:', vizUrl);
                    
                    // Fetch the HTML file
                    const vizResponse = await fetch(vizUrl);
                    if (vizResponse.ok) {
                        const htmlContent = await vizResponse.text();
                        
                        // Create an iframe to display the Plotly HTML
                        const plotDiv = document.getElementById('plot');
                        plotDiv.innerHTML = ''; // Clear existing content
                        
                        const iframe = document.createElement('iframe');
                        iframe.style.width = '100%';
                        iframe.style.height = '100%';
                        iframe.style.border = 'none';
                        iframe.srcdoc = htmlContent;
                        plotDiv.appendChild(iframe);
                        
                        // Show visualization controls
                        document.getElementById('visualization-section').style.display = 'block';
                        document.getElementById('search-section').style.display = 'block';
                        document.getElementById('stats-section').style.display = 'block';
                        
                        showStatus('✓ Visualization loaded and displayed!', 'success');
                    } else {
                        throw new Error('Failed to load visualization file');
                    }
                } catch (error) {
                    console.error('[CORE] Error loading visualization:', error);
                    showStatus(`Visualization generated but could not be loaded: ${error.message}. File available at: ${data.visualization.file_url}`, 'info');
                    
                    // Fallback: try to load coremap embedding data instead
                    if (data.key_added && data.obsm_keys && data.obsm_keys.includes(data.key_added)) {
                        console.log('[CORE] Attempting to load coremap embedding data...');
                        // Trigger visualization with coremap data
                        loadCoremapVisualization(data.key_added, data.cluster_key);
                    }
                }
            } else if (data.visualization.error) {
                console.warn('Visualization warning:', data.visualization.error);
                // Try to load coremap embedding if available
                if (data.key_added && data.obsm_keys && data.obsm_keys.includes(data.key_added)) {
                    loadCoremapVisualization(data.key_added, data.cluster_key);
                }
            }
        } else if (data.key_added && data.obsm_keys && data.obsm_keys.includes(data.key_added)) {
            // If visualization wasn't requested but coremap was generated, load it
            loadCoremapVisualization(data.key_added, data.cluster_key);
        }
        
        // Update cluster info if clustering was auto-run
        if (data.clustering_auto_run && data.cluster_key) {
            clusterInfo = {
                cluster_key: data.cluster_key,
                clusters: data.clusters || [],
                cluster_counts: data.cluster_counts || {},
                has_model: true
            };
            
            // Update available cluster keys for DEG analysis
            await updateAvailableClusterKeys();
            
            // Update color-by select
            updateMetadataSelect([data.cluster_key]);
        }
        
        // Enable downstream steps
        document.getElementById('btn-visualize').disabled = false;
        document.getElementById('btn-marker-genes').disabled = false;
        document.getElementById('deg-step').style.display = 'block';
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Core selection error:', error);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

async function runMarkerGenes() {
    showStatus('Finding marker genes...', 'info');
    const btn = document.getElementById('btn-marker-genes');
    const originalText = btn.textContent;
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner"></span> Analyzing...';

    try {
        const cluster_key = document.getElementById('deg-cluster-key').value || clusterInfo.cluster_key;

        const response = await fetch(`${API_BASE}/marker-genes`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                session_id: sessionId,
                cluster_key: cluster_key
            })
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Marker genes analysis failed');
        }

        console.log('Marker genes response:', data);
        console.log('Markers type:', typeof data.markers);
        console.log('Markers keys:', data.markers ? Object.keys(data.markers).slice(0, 5) : 'null');
        console.log('n_markers:', data.n_markers);
        console.log('Stats:', data.stats);

        if (!data.markers || Object.keys(data.markers).length === 0) {
            showStatus('⚠ No marker genes found. This may be normal if clusters are too similar.', 'info');
            displayMarkerGenes({}, 0, data.stats);
            return;
        }

        // Show status with statistics matching lotus_workflow.py output
        let statusMsg = '✓ Marker genes analysis complete!';
        if (data.stats) {
            statusMsg += ` (Total: ${data.stats.total}, p_adj<0.05: ${data.stats.significant_p05}, p_adj<0.01: ${data.stats.significant_p01})`;
        }
        showStatus(statusMsg, 'success');
        console.log('Marker genes:', data.markers);
        
        // Display marker genes results in table
        displayMarkerGenes(data.markers, data.n_markers, data.stats);
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Marker genes error:', error);
    } finally {
        btn.disabled = false;
        btn.textContent = originalText;
    }
}

function displayMarkerGenes(markers, nMarkers, stats) {
    const container = document.getElementById('marker-genes-results');
    const tableBody = document.getElementById('marker-genes-table-body');
    
    if (!markers || Object.keys(markers).length === 0) {
        let message = 'No marker genes found.';
        if (stats && stats.total > 0) {
            message = `No marker genes in top results (Total: ${stats.total}, p_adj<0.05: ${stats.significant_p05}, p_adj<0.01: ${stats.significant_p01})`;
        }
        tableBody.innerHTML = `<tr><td colspan="4" style="padding: 16px; text-align: center; color: #64748b;">${message}</td></tr>`;
        container.style.display = 'block';
        return;
    }
    
    // Display statistics if available (matching lotus_workflow.py output)
    if (stats) {
        console.log(`DEG Statistics (matching lotus_workflow.py):`);
        console.log(`  - Total differentially expressed genes: ${stats.total}`);
        console.log(`  - Significant genes (p_adj < 0.05): ${stats.significant_p05}`);
        console.log(`  - Significant genes (p_adj < 0.01): ${stats.significant_p01}`);
    }
    
    // Clear existing rows
    tableBody.innerHTML = '';
    
    // Convert markers dict to array and sort by absolute log2fc (descending)
    const markersArray = Object.values(markers).map(m => {
        // Ensure all values are valid numbers
        const log2fc = typeof m.log2fc === 'number' ? m.log2fc : parseFloat(m.log2fc) || 0;
        const pval = typeof m.pval === 'number' ? m.pval : parseFloat(m.pval) || 1;
        const pval_adj = typeof m.pval_adj === 'number' ? m.pval_adj : parseFloat(m.pval_adj) || 1;
        return {
            gene: m.gene || 'unknown',
            log2fc: log2fc,
            pval: pval,
            pval_adj: pval_adj
        };
    }).filter(m => m.gene !== 'unknown'); // Filter out invalid entries
    
    if (markersArray.length === 0) {
        tableBody.innerHTML = '<tr><td colspan="4" style="padding: 16px; text-align: center; color: #64748b;">No valid marker genes found.</td></tr>';
        container.style.display = 'block';
        return;
    }
    
    markersArray.sort((a, b) => Math.abs(b.log2fc) - Math.abs(a.log2fc));
    
    // Display top markers
    markersArray.forEach(marker => {
        const row = document.createElement('tr');
        row.style.borderBottom = '1px solid var(--border)';
        
        // Color coding based on log2fc
        const log2fcColor = marker.log2fc > 0 ? '#10b981' : '#ef4444';
        const pvalColor = marker.pval_adj < 0.05 ? '#1e40af' : '#64748b';
        
        row.innerHTML = `
            <td style="padding: 8px; font-weight: 500;">${marker.gene}</td>
            <td style="padding: 8px; text-align: right; color: ${log2fcColor}; font-weight: 500;">
                ${marker.log2fc.toFixed(3)}
            </td>
            <td style="padding: 8px; text-align: right; color: ${pvalColor};">
                ${marker.pval.toExponential(2)}
            </td>
            <td style="padding: 8px; text-align: right; color: ${pvalColor};">
                ${marker.pval_adj.toExponential(2)}
            </td>
        `;
        
        tableBody.appendChild(row);
    });
    
    // Show the results section
    container.style.display = 'block';
    
    // Scroll to results
    container.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

async function loadCoremapVisualization(coremapKey, clusterKey) {
    // Load coremap embedding data and display it
    try {
        showStatus('Loading coremap visualization...', 'info');
        
        // Get metadata to find available embeddings
        const metadataResponse = await fetch(`${API_BASE}/metadata?session_id=${sessionId}`);
        const metadataData = await metadataResponse.json();
        
        if (!metadataResponse.ok) {
            throw new Error(metadataData.error || 'Failed to get metadata');
        }
        
        // Check if coremap embedding exists
        if (!metadataData.obsm_keys || !metadataData.obsm_keys.includes(coremapKey)) {
            throw new Error(`Coremap embedding '${coremapKey}' not found`);
        }
        
        // Use visualize API but it will use existing embeddings
        // We need to create a custom endpoint or modify the approach
        // For now, let's use the visualize endpoint which will compute UMAP
        // But we can also directly fetch and use the coremap coordinates
        
        // Use visualize API with coremap embedding
        const response = await fetch(`${API_BASE}/visualize`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                session_id: sessionId,
                cluster_key: clusterKey || 'cplearn',
                n_components: 2,
                use_rep: coremapKey  // Use coremap embedding directly
            })
        });
        
        const data = await response.json();
        
        if (!response.ok) {
            throw new Error(data.error || 'Failed to load visualization');
        }
        
        // Check if we can get coremap coordinates directly
        // For now, use UMAP coordinates but note that coremap is available
        vectors = data.coordinates;
        metadata = data.metadata;
        
        // Add note about coremap availability
        if (metadataData.obsm_keys.includes(coremapKey)) {
            showStatus(`✓ Visualization loaded! (Coremap '${coremapKey}' available in data)`, 'success');
        } else {
            showStatus('✓ Visualization loaded!', 'success');
        }
        
        updatePlot();
        
    } catch (error) {
        console.error('[CORE] Error loading coremap visualization:', error);
        showStatus(`Error loading visualization: ${error.message}`, 'error');
    }
}

function showStatus(message, type) {
    const statusDiv = document.getElementById('status');
    statusDiv.textContent = message;
    statusDiv.className = `status ${type}`;
    statusDiv.style.display = 'block';
    
    if (type === 'success') {
        setTimeout(() => {
            statusDiv.style.display = 'none';
        }, 5000);
    }
}
