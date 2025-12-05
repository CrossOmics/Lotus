<<<<<<< HEAD
// Lotus - Main Application Logic with Backend Integration

const API_BASE = '/api';
let vectors = null;
let metadata = null;
let currentPlot = null;
let colorMapping = {};
// Generate unique session ID for each user
let sessionId = generateSessionId();
let clusterInfo = {}; // Store cluster information
let heartbeatInterval = null; // Heartbeat interval ID
let sessionCleanupSent = false; // Flag to prevent multiple cleanup requests
let visualizationCache = {}; // Store visualization results: { cache_key: { vectors, metadata, params } }
let coreAnalysisVisualization = null; // Store Core Analysis visualization data: { file_url, htmlContent, iframe }
// Cache key format: `${cluster_key}_${n_components}_${min_dist}_${spread}`

// Function to generate cache key from cluster key, method, and visualization parameters
function getVisualizationCacheKey(cluster_key, method, n_components, min_dist, spread) {
    if (method === 'coremap') {
        return `${cluster_key}_coremap_${n_components}`;
    } else if (method === 'umap') {
        return `${cluster_key}_umap_${n_components}_${min_dist}_${spread}`;
    } else {
        // For tsne, diffmap, draw_graph - use method and n_components
        return `${cluster_key}_${method}_${n_components}`;
    }
}

// Function to get current visualization parameters
function getVisualizationParams() {
    const method = document.getElementById('viz-method')?.value || 'umap';
    const n_components = document.getElementById('dimension').value === '3d' ? 3 : 2;
    const min_dist = parseFloat(document.getElementById('umap-min-dist').value) || 0.5;
    const spread = parseFloat(document.getElementById('umap-spread').value) || 1.0;
    return { method, n_components, min_dist, spread };
}

// Function to update visualization method selector based on Core Analysis availability
function updateVisualizationMethodSelector() {
    const methodSelect = document.getElementById('viz-method');
    if (!methodSelect) return;
    
    // Check if Core Analysis has been run
    // Use clusterInfo.has_model which is set when Core Analysis or cplearn clustering is run
    // This is preserved even when running scanpy clustering after Core Analysis
    const hasCoreAnalysis = clusterInfo && clusterInfo.has_model;
    
    // Get current value
    const currentValue = methodSelect.value;
    
    // Clear options
    methodSelect.innerHTML = '';
    
    // Always add UMAP
    const umapOption = document.createElement('option');
    umapOption.value = 'umap';
    umapOption.textContent = 'UMAP';
    methodSelect.appendChild(umapOption);
    
    // Add Coremap if Core Analysis is available
    if (hasCoreAnalysis) {
        const coremapOption = document.createElement('option');
        coremapOption.value = 'coremap';
        coremapOption.textContent = 'Coremap';
        methodSelect.appendChild(coremapOption);
    }
    
    // Restore previous selection if still valid, otherwise default to UMAP
    if (currentValue && Array.from(methodSelect.options).some(opt => opt.value === currentValue)) {
        methodSelect.value = currentValue;
    } else {
        methodSelect.value = 'umap';
    }
    
    // Update parameters visibility
    updateUMAPParamsVisibility();
}

// Function to update parameters visibility based on selected method
function updateUMAPParamsVisibility() {
    const methodSelect = document.getElementById('viz-method');
    const umapParamsGroup = document.getElementById('umap-params-group');
    const umapSpreadGroup = document.getElementById('umap-spread-group');
    const vizClusterKeyGroup = document.getElementById('viz-cluster-key')?.closest('.form-group');
    const dimensionGroup = document.getElementById('dimension')?.closest('.form-group');
    const coremapGroundTruthGroup = document.getElementById('coremap-ground-truth-group');
    const coremapFastViewGroup = document.getElementById('coremap-fast-view-group');
    
    if (methodSelect) {
        const method = methodSelect.value;
        if (method === 'coremap') {
            // Hide UMAP parameters, Cluster Key, and Dimension for coremap
            if (umapParamsGroup) umapParamsGroup.style.display = 'none';
            if (umapSpreadGroup) umapSpreadGroup.style.display = 'none';
            if (vizClusterKeyGroup) vizClusterKeyGroup.style.display = 'none';
            if (dimensionGroup) dimensionGroup.style.display = 'none';
            // Show Fast View and Ground Truth upload for coremap
            if (coremapFastViewGroup) coremapFastViewGroup.style.display = 'block';
            if (coremapGroundTruthGroup) coremapGroundTruthGroup.style.display = 'block';
        } else {
            // Show parameters for standard methods (umap, tsne, diffmap, draw_graph)
            // UMAP-specific parameters only shown for UMAP
            if (method === 'umap') {
                if (umapParamsGroup) umapParamsGroup.style.display = 'block';
                if (umapSpreadGroup) umapSpreadGroup.style.display = 'block';
            } else {
                // Hide UMAP-specific parameters for other methods
                if (umapParamsGroup) umapParamsGroup.style.display = 'none';
                if (umapSpreadGroup) umapSpreadGroup.style.display = 'none';
            }
            // Show Cluster Key and Dimension for all standard methods
            if (vizClusterKeyGroup) vizClusterKeyGroup.style.display = 'block';
            if (dimensionGroup) dimensionGroup.style.display = 'block';
            // Hide Fast View and Ground Truth upload for standard methods
            if (coremapFastViewGroup) coremapFastViewGroup.style.display = 'none';
            if (coremapGroundTruthGroup) coremapGroundTruthGroup.style.display = 'none';
        }
    }
}

// Generate unique session ID using UUID v4
function generateSessionId() {
    // Use crypto.randomUUID if available (modern browsers)
    if (typeof crypto !== 'undefined' && crypto.randomUUID) {
        return crypto.randomUUID();
    }
    // Fallback: generate UUID v4 manually
    return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
        const r = Math.random() * 16 | 0;
        const v = c === 'x' ? r : (r & 0x3 | 0x8);
        return v.toString(16);
    });
}

// Initialize session on page load
async function initializeSession() {
    try {
        const response = await fetch(`${API_BASE}/session/create`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId })
        });
        if (response.ok) {
            console.log(`[SESSION] Session initialized: ${sessionId}`);
            // Start heartbeat
            startHeartbeat();
        } else {
            console.warn('[SESSION] Failed to initialize session');
        }
    } catch (error) {
        console.error('[SESSION] Error initializing session:', error);
    }
}

// Send heartbeat to keep session alive
async function sendHeartbeat() {
    try {
        await fetch(`${API_BASE}/session/heartbeat`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId })
        });
    } catch (error) {
        console.error('[SESSION] Heartbeat error:', error);
    }
}

// Start heartbeat (send every 30 seconds)
function startHeartbeat() {
    // Send initial heartbeat
    sendHeartbeat();
    // Set up interval
    heartbeatInterval = setInterval(sendHeartbeat, 30000); // 30 seconds
}

// Stop heartbeat
function stopHeartbeat() {
    if (heartbeatInterval) {
        clearInterval(heartbeatInterval);
        heartbeatInterval = null;
    }
}

// Clean up session when page is closed or hidden
async function cleanupSession() {
    if (sessionCleanupSent) return; // Prevent multiple cleanup requests
    sessionCleanupSent = true;
    
    stopHeartbeat();
    
    try {
        // Use sendBeacon for reliable cleanup on page unload
        const data = JSON.stringify({ session_id: sessionId });
        if (navigator.sendBeacon) {
            navigator.sendBeacon(`${API_BASE}/session/cleanup`, data);
            console.log(`[SESSION] Session cleanup sent via beacon: ${sessionId}`);
        } else {
            // Fallback: use fetch with keepalive
            await fetch(`${API_BASE}/session/cleanup`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: data,
                keepalive: true
            });
            console.log(`[SESSION] Session cleanup sent: ${sessionId}`);
        }
    } catch (error) {
        console.error('[SESSION] Error cleaning up session:', error);
    }
}

// Set up cleanup handlers
window.addEventListener('beforeunload', cleanupSession);
document.addEventListener('visibilitychange', function() {
    if (document.visibilityState === 'hidden') {
        // Page is hidden, but don't cleanup yet (user might come back)
        // Just stop heartbeat to save resources
        stopHeartbeat();
    } else if (document.visibilityState === 'visible') {
        // Page is visible again, restart heartbeat
        if (!heartbeatInterval) {
            startHeartbeat();
        }
    }
});

// Initialize session when page loads
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initializeSession);
} else {
    initializeSession();
}

// Function to disable all buttons during operations
// Note: File upload and default dataset buttons are NOT disabled to allow reloading data at any time
function disableAllButtons() {
    const buttons = [
        'btn-preprocess',
        'btn-cluster',
        'btn-core-selection',
        'btn-visualize',
        'btn-marker-genes'
    ];
    buttons.forEach(btnId => {
        const btn = document.getElementById(btnId);
        if (btn) {
            btn.disabled = true;
        }
    });
}

// Function to disable data loading buttons during loading
function disableDataLoadingButtons() {
    const buttons = [
        'load-demo-data',
        'load-pbmc3k',
        'data-file'  // File input (though it's hidden, we can disable it)
    ];
    buttons.forEach(btnId => {
        const btn = document.getElementById(btnId);
        if (btn) {
            btn.disabled = true;
        }
    });
}

// Function to enable data loading buttons
function enableDataLoadingButtons() {
    const buttons = [
        'load-demo-data',
        'load-pbmc3k',
        'data-file'
    ];
    buttons.forEach(btnId => {
        const btn = document.getElementById(btnId);
        if (btn) {
            btn.disabled = false;
        }
    });
}

// Function to update button states based on current data state
// Rule 2: Visualization and Marker Genes require either Core Analysis or Clustering to be completed
function updateButtonStates() {
    // Check if we have clustering or core analysis results
    const hasClustering = clusterInfo && clusterInfo.cluster_key;
    const hasCoreAnalysis = clusterInfo && clusterInfo.has_model;
    
    // Rule 2: Visualization and Marker Genes require clustering or core analysis
    const canVisualize = hasClustering || hasCoreAnalysis;
    
    const btnVisualize = document.getElementById('btn-visualize');
    const btnMarkerGenes = document.getElementById('btn-marker-genes');
    
    if (btnVisualize) {
        btnVisualize.disabled = !canVisualize;
    }
    if (btnMarkerGenes) {
        btnMarkerGenes.disabled = !canVisualize;
    }
    
    // Update cluster key select for Marker Genes
    const degClusterKeySelect = document.getElementById('deg-cluster-key');
    if (degClusterKeySelect && canVisualize) {
        // If we can visualize, the select should be enabled if cluster keys are available
        // This will be handled by updateAvailableClusterKeys
    } else if (degClusterKeySelect && !canVisualize) {
        degClusterKeySelect.disabled = true;
    }
}

// Function to reset all state when new data is loaded
function resetAllState() {
    console.log('[RESET] Resetting all state...');
    
    // Clear cluster info
    clusterInfo = {};
    
    // Clear visualization cache
    visualizationCache = {};
    
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
    
    // Reset cluster keys dropdowns
    const degClusterKeySelect = document.getElementById('deg-cluster-key');
    if (degClusterKeySelect) {
        degClusterKeySelect.innerHTML = '<option value="">No cluster keys available - Run clustering first</option>';
        degClusterKeySelect.disabled = true;
    }
    const vizClusterKeySelect = document.getElementById('viz-cluster-key');
    if (vizClusterKeySelect) {
        vizClusterKeySelect.innerHTML = '<option value="">No cluster keys available - Run clustering first</option>';
        vizClusterKeySelect.disabled = true;
    }
    
    // DEG step is always visible (same as visualization)
    
    // Clear coremap visualization if exists
    const coremapViz = document.getElementById('coremap-visualization');
    if (coremapViz) {
        coremapViz.innerHTML = '';
        coremapViz.style.display = 'none';
    }
    
    // Clear Core Analysis visualization data (if exists, for caching)
    coreAnalysisVisualization = null;
    
    // Clear Ground Truth file
    clearVizGroundTruthFile();
    
    console.log('[RESET] All state reset complete');
    
    // Update button states based on rules
    updateButtonStates();
}

// Function to update available cluster keys dropdown
async function updateAvailableClusterKeys() {
    try {
        const response = await fetch(`${API_BASE}/cluster-keys`, {
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

        // Update both cluster key selectors (DEG and Visualization)
        const degSelect = document.getElementById('deg-cluster-key');
        const vizSelect = document.getElementById('viz-cluster-key');

        // Update DEG cluster key selector
        if (degSelect) {
            degSelect.innerHTML = '';
        if (data.cluster_keys && data.cluster_keys.length > 0) {
            data.cluster_keys.forEach(clusterKeyInfo => {
                const option = document.createElement('option');
                option.value = clusterKeyInfo.key;
                const label = `${clusterKeyInfo.key} (${clusterKeyInfo.method}, ${clusterKeyInfo.n_clusters} clusters)`;
                option.textContent = label;
                option.selected = clusterKeyInfo.key === clusterInfo.cluster_key;
                    degSelect.appendChild(option);
                });
                degSelect.disabled = false;
            } else {
                const option = document.createElement('option');
                option.value = '';
                option.textContent = 'No cluster keys available - Run clustering first';
                option.disabled = true;
                degSelect.appendChild(option);
                degSelect.disabled = true;
            }
        }

        // Update Visualization cluster key selector
        if (vizSelect) {
            vizSelect.innerHTML = '';
            if (data.cluster_keys && data.cluster_keys.length > 0) {
                data.cluster_keys.forEach(clusterKeyInfo => {
                    const option = document.createElement('option');
                    option.value = clusterKeyInfo.key;
                    const label = `${clusterKeyInfo.key} (${clusterKeyInfo.method}, ${clusterKeyInfo.n_clusters} clusters)`;
                    option.textContent = label;
                    option.selected = clusterKeyInfo.key === clusterInfo.cluster_key;
                    vizSelect.appendChild(option);
                });
                vizSelect.disabled = false;
                // Do not auto-trigger visualization when updating cluster keys
                // User should manually select cluster key or click "Generate Visualization" button
        } else {
            const option = document.createElement('option');
            option.value = '';
            option.textContent = 'No cluster keys available - Run clustering first';
            option.disabled = true;
                vizSelect.appendChild(option);
                vizSelect.disabled = true;
            }
        }
        
        // Update button states based on rules after cluster keys are updated
        updateButtonStates();
        
        // Update visualization method selector (to show/hide coremap option)
        updateVisualizationMethodSelector();
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
    try {
        console.log('[INIT] Initializing application...');
    checkHealth();
    setupFileInputs();
        setupDefaultDatasets();
    setupControls();
    setupPipelineButtons();
    // Initialize cluster keys dropdown (will be empty until clustering is run)
    updateAvailableClusterKeys();
        console.log('[INIT] Application initialized successfully');
    } catch (error) {
        console.error('[INIT] Error during initialization:', error);
        showStatus('Error initializing application. Please refresh the page.', 'error');
    }
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
    
    // Setup Ground Truth file upload for Visualization (coremap)
    setupGroundTruthFileInput();
}

function setupGroundTruthFileInput() {
    try {
        // Setup for Visualization Ground Truth upload (coremap)
        const vizGroundTruthInput = document.getElementById('viz-truth-json-file');
        const vizGroundTruthUploadArea = document.getElementById('viz-ground-truth-upload-area');
        
        if (!vizGroundTruthInput || !vizGroundTruthUploadArea) {
            // Elements don't exist yet (might be hidden), that's okay
            return;
        }
        
        // Drag and drop for Ground Truth
        vizGroundTruthUploadArea.addEventListener('dragover', (e) => {
            e.preventDefault();
            e.stopPropagation();
            vizGroundTruthUploadArea.classList.add('dragover');
        });
        
        vizGroundTruthUploadArea.addEventListener('dragleave', () => {
            vizGroundTruthUploadArea.classList.remove('dragover');
        });
        
        vizGroundTruthUploadArea.addEventListener('drop', (e) => {
            e.preventDefault();
            e.stopPropagation();
            vizGroundTruthUploadArea.classList.remove('dragover');
            
            const files = e.dataTransfer.files;
            if (files.length > 0 && files[0].type === 'application/json') {
                vizGroundTruthInput.files = files;
                
                // Update UI to show selected file
                updateVizGroundTruthFileDisplay(files[0].name);
            }
        });
        
        vizGroundTruthUploadArea.addEventListener('click', () => {
            vizGroundTruthInput.click();
        });
        
        // Update UI when file is selected via click
        vizGroundTruthInput.addEventListener('change', (e) => {
            if (e.target.files && e.target.files.length > 0) {
                updateVizGroundTruthFileDisplay(e.target.files[0].name);
            }
        });
        
        // Setup remove button
        const removeBtn = document.getElementById('viz-ground-truth-remove-btn');
        if (removeBtn) {
            removeBtn.addEventListener('click', (e) => {
                e.stopPropagation(); // Prevent triggering file input click
                clearVizGroundTruthFile();
            });
        }
    } catch (error) {
        console.error('[SETUP] Error setting up Ground Truth file input:', error);
    }
}

function updateVizGroundTruthFileDisplay(filename) {
    const uploadArea = document.getElementById('viz-ground-truth-upload-area');
    if (!uploadArea) return;
    
    const textElement = uploadArea.querySelector('.file-upload-text');
    const hintElement = uploadArea.querySelector('.file-upload-hint');
    const removeBtn = document.getElementById('viz-ground-truth-remove-btn');
    
    if (textElement) {
        textElement.textContent = filename;
        textElement.style.color = 'var(--primary)';
    }
    if (hintElement) {
        hintElement.textContent = 'Click to change file';
    }
    if (removeBtn) {
        removeBtn.style.display = 'block';
    }
}

function clearVizGroundTruthFile() {
    const vizGroundTruthInput = document.getElementById('viz-truth-json-file');
    const uploadArea = document.getElementById('viz-ground-truth-upload-area');
    if (!uploadArea || !vizGroundTruthInput) return;
    
    // Clear file input
    vizGroundTruthInput.value = '';
    
    // Reset display
    const textElement = uploadArea.querySelector('.file-upload-text');
    const hintElement = uploadArea.querySelector('.file-upload-hint');
    const removeBtn = document.getElementById('viz-ground-truth-remove-btn');
    
    if (textElement) {
        textElement.textContent = 'Click or drag JSON file';
        textElement.style.color = '';
    }
    if (hintElement) {
        hintElement.textContent = '.json';
    }
    if (removeBtn) {
        removeBtn.style.display = 'none';
    }
}

function setupDefaultDatasets() {
    const loadDemoDataBtn = document.getElementById('load-demo-data');
    const loadPbmc3kBtn = document.getElementById('load-pbmc3k');
    
    if (loadDemoDataBtn) {
        loadDemoDataBtn.addEventListener('click', () => {
            loadDefaultDataset('demo_data.h5ad');
        });
    }
    
    if (loadPbmc3kBtn) {
        loadPbmc3kBtn.addEventListener('click', () => {
            loadDefaultDataset('pbmc3k_raw.h5ad');
        });
    }
}

async function loadDefaultDataset(filename) {
    console.log(`[LOAD] Loading default dataset: ${filename}`);
    showStatus(`Loading ${filename}...`, 'info');
    
    // Reset all state when loading new data
    resetAllState();
    
    // Disable pipeline buttons and data loading buttons during loading
    disableAllButtons();
    disableDataLoadingButtons();
    
    try {
        const response = await fetch(`${API_BASE}/load-default-dataset`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                filename: filename,
                session_id: sessionId
            })
        });

        let data;
        try {
            data = await response.json();
            console.log('[LOAD] Response data:', data);
        } catch (e) {
            const text = await response.text();
            console.error('[LOAD] Failed to parse JSON:', text);
            throw new Error(`Server error: ${response.status} ${response.statusText}. Response: ${text.substring(0, 200)}`);
        }

        if (!response.ok) {
            console.error('[LOAD] Error response:', data);
            throw new Error(data.error || `Load failed: ${response.status}`);
        }

        console.log('[LOAD] Success:', data);
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
        
        // Update available cluster keys (in case loaded data already has clustering results)
        await updateAvailableClusterKeys();
        
        // Update button states based on rules (check if loaded data has clustering/core analysis)
        updateButtonStates();
        
        // Update visualization method selector
        updateVisualizationMethodSelector();
        
    } catch (error) {
        console.error('[LOAD] Exception:', error);
        showStatus(`Error: ${error.message}`, 'error');
    } finally {
        // Always re-enable data loading buttons to allow reloading at any time
        enableDataLoadingButtons();
    }
}

function setupControls() {
    // Initialize visualization method selector
    updateVisualizationMethodSelector();
    
    // Add event listener for visualization method selector
    const vizMethodSelect = document.getElementById('viz-method');
    if (vizMethodSelect) {
        vizMethodSelect.addEventListener('change', (e) => {
            // Update UMAP parameters visibility
            updateUMAPParamsVisibility();
            
            // Clear cache for current cluster key when method changes
            const vizClusterKeySelect = document.getElementById('viz-cluster-key');
            const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
            if (selectedKey && selectedKey !== '') {
                Object.keys(visualizationCache).forEach(key => {
                    if (key.startsWith(selectedKey + '_')) {
                        delete visualizationCache[key];
                    }
                });
            }
            
            // Don't auto-generate visualization when switching to coremap
            // User needs to click "Generate Visualization" button
        });
    }
    
    document.getElementById('color-by').addEventListener('change', updatePlot);
    document.getElementById('dimension').addEventListener('change', (e) => {
        // When dimension changes, clear current visualization cache for the selected cluster key
        // This forces regeneration with new parameters
        const vizClusterKeySelect = document.getElementById('viz-cluster-key');
        const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
        if (selectedKey && selectedKey !== '') {
            // Clear cache entries for this cluster key (they will be regenerated with new parameters)
            Object.keys(visualizationCache).forEach(key => {
                if (key.startsWith(selectedKey + '_')) {
                    delete visualizationCache[key];
                }
            });
        }
        updatePlot();
    });
    
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
    
    // Update cplearn advanced parameters visibility based on method selection
    const clusterMethodSelect = document.getElementById('cluster-method');
    const cplearnAdvancedParams = document.getElementById('cplearn-advanced-params');
    
    function updateCplearnParamsVisibility() {
        if (clusterMethodSelect && cplearnAdvancedParams) {
            const method = clusterMethodSelect.value;
            if (method === 'cplearn') {
                cplearnAdvancedParams.style.display = 'block';
            } else {
                cplearnAdvancedParams.style.display = 'none';
            }
        }
    }
    
    if (clusterMethodSelect) {
        clusterMethodSelect.addEventListener('change', updateCplearnParamsVisibility);
        // Initialize visibility on page load
        updateCplearnParamsVisibility();
    }
    
    // Add event listeners for visualization parameters to clear cache when changed
    const umapMinDistInput = document.getElementById('umap-min-dist');
    const umapSpreadInput = document.getElementById('umap-spread');
    
    if (umapMinDistInput) {
        umapMinDistInput.addEventListener('change', () => {
            // Clear cache for current cluster key when min_dist changes
            const vizClusterKeySelect = document.getElementById('viz-cluster-key');
            const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
            if (selectedKey && selectedKey !== '') {
                Object.keys(visualizationCache).forEach(key => {
                    if (key.startsWith(selectedKey + '_')) {
                        delete visualizationCache[key];
                    }
                });
            }
        });
    }
    
    if (umapSpreadInput) {
        umapSpreadInput.addEventListener('change', () => {
            // Clear cache for current cluster key when spread changes
            const vizClusterKeySelect = document.getElementById('viz-cluster-key');
            const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
            if (selectedKey && selectedKey !== '') {
                Object.keys(visualizationCache).forEach(key => {
                    if (key.startsWith(selectedKey + '_')) {
                        delete visualizationCache[key];
                    }
                });
            }
        });
    }
    
    // Add event listener for Visualization cluster key selector
    const vizClusterKeySelect = document.getElementById('viz-cluster-key');
    if (vizClusterKeySelect) {
        vizClusterKeySelect.addEventListener('change', async (e) => {
            const selectedKey = e.target.value;
            if (selectedKey && selectedKey !== '') {
                // Automatically update Color By to match the selected cluster key
                const colorBySelect = document.getElementById('color-by');
                if (colorBySelect) {
                    // Check if the cluster key exists in the color-by options
                    const optionExists = Array.from(colorBySelect.options).some(opt => opt.value === selectedKey);
                    if (!optionExists) {
                        // Add the cluster key to color-by options if it doesn't exist
                        const option = document.createElement('option');
                        option.value = selectedKey;
                        option.textContent = selectedKey;
                        colorBySelect.appendChild(option);
                    }
                    // Set color-by to the selected cluster key
                    colorBySelect.value = selectedKey;
                }
                
                // Automatically generate new visualization when switching cluster key
                // Cache will be checked inside runVisualizationForClusterKey with current parameters
                await runVisualizationForClusterKey(selectedKey);
            } else {
                // Clear visualization if no cluster key selected
                vectors = null;
                metadata = null;
                const plotDiv = document.getElementById('plot');
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #64748b;">Select a cluster key to visualize</div>';
                }
            }
        });
    }
}

function setupPipelineButtons() {
    try {
        const btnPreprocess = document.getElementById('btn-preprocess');
        const btnCluster = document.getElementById('btn-cluster');
        const btnVisualize = document.getElementById('btn-visualize');
        const btnCoreSelection = document.getElementById('btn-core-selection');
        const btnMarkerGenes = document.getElementById('btn-marker-genes');
        
        if (btnPreprocess) {
            btnPreprocess.addEventListener('click', runPreprocess);
        } else {
            console.error('[SETUP] btn-preprocess not found');
        }
        
        if (btnCluster) {
            btnCluster.addEventListener('click', runClustering);
        } else {
            console.error('[SETUP] btn-cluster not found');
        }
        
        if (btnVisualize) {
            btnVisualize.addEventListener('click', runVisualization);
        } else {
            console.error('[SETUP] btn-visualize not found');
        }
        
        if (btnCoreSelection) {
            btnCoreSelection.addEventListener('click', runCoreSelection);
        } else {
            console.error('[SETUP] btn-core-selection not found');
        }
        
        if (btnMarkerGenes) {
            btnMarkerGenes.addEventListener('click', runMarkerGenes);
        } else {
            console.error('[SETUP] btn-marker-genes not found');
        }
    } catch (error) {
        console.error('[SETUP] Error setting up pipeline buttons:', error);
    }
}

async function handleDataFile(event) {
    const file = event.target.files[0];
    if (!file) {
        console.log('[UPLOAD] No file selected');
        return;
    }

    console.log('[UPLOAD] File selected:', file.name, 'Size:', file.size, 'bytes');
    showStatus('Uploading data...', 'info');
    
    // Reset all state when loading new data
    resetAllState();
    
    // Disable pipeline buttons and data loading buttons during upload
    disableAllButtons();
    disableDataLoadingButtons();
    
    // Determine file type from extension (only h5ad, csv, tsv allowed)
    const fileName = file.name.toLowerCase();
    let fileType = null;
    
    if (fileName.endsWith('.h5ad')) {
        fileType = 'h5ad';
    } else if (fileName.endsWith('.csv')) {
        fileType = 'csv';
    } else if (fileName.endsWith('.tsv')) {
        fileType = 'tsv';
    } else {
        showStatus('Unsupported file format. Please upload .h5ad, .csv, or .tsv files only.', 'error');
        enableDataLoadingButtons(); // Re-enable buttons if file format is invalid
        return;
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
        
        // Update button states based on rules (check if uploaded data has clustering/core analysis)
        updateButtonStates();
        
        // Update visualization method selector
        updateVisualizationMethodSelector();
        
    } catch (error) {
        console.error('[UPLOAD] Exception:', error);
        showStatus(`Error: ${error.message}`, 'error');
    } finally {
        // Always re-enable data loading buttons to allow reloading at any time
        enableDataLoadingButtons();
    }
}

async function runPreprocess() {
    showStatus('Running preprocessing... This may take a while.', 'info');
    const btn = document.getElementById('btn-preprocess');
    const originalText = btn.textContent;
    
    // Disable all buttons during preprocessing
    disableAllButtons();
    btn.innerHTML = '<span class="spinner"></span> Processing...';

    try {
        const n_pcs_value = document.getElementById('n-pcs').value?.trim();
        const n_pcs = n_pcs_value ? parseInt(n_pcs_value) : null;
        const n_neighbors = parseInt(document.getElementById('n-neighbors').value) || 15;
        const target_sum = parseFloat(document.getElementById('target-sum')?.value) || 1e4;
        const n_top_genes = document.getElementById('n-top-genes')?.value ? parseInt(document.getElementById('n-top-genes').value) : null;
        const use_rep = document.getElementById('preprocess-use-rep')?.value || 'X_pca';
        const save_raw = document.getElementById('save-raw')?.checked !== false;
        const min_genes = document.getElementById('min-genes')?.value ? parseInt(document.getElementById('min-genes').value) : null;
        const min_cells = document.getElementById('min-cells')?.value ? parseInt(document.getElementById('min-cells').value) : null;
        const min_counts = document.getElementById('min-counts')?.value ? parseInt(document.getElementById('min-counts').value) : null;
        const max_counts = document.getElementById('max-counts')?.value ? parseInt(document.getElementById('max-counts').value) : null;
        const max_genes = document.getElementById('max-genes')?.value ? parseInt(document.getElementById('max-genes').value) : null;
        const pct_mt_max = document.getElementById('pct-mt-max')?.value ? parseFloat(document.getElementById('pct-mt-max').value) : null;
        const hvg_flavor = document.getElementById('hvg-flavor')?.value || 'seurat_v3';
        const batch_key = document.getElementById('batch-key')?.value?.trim() || null;
        const regress_out_keys_str = document.getElementById('regress-out-keys')?.value?.trim() || null;
        const regress_out_keys = regress_out_keys_str ? regress_out_keys_str.split(',').map(k => k.trim()).filter(k => k) : null;
        const use_combat = document.getElementById('use-combat')?.checked || false;

        const requestBody = {
            session_id: sessionId,
            n_neighbors: n_neighbors,
            target_sum: target_sum,
            use_rep: use_rep,
            save_raw: save_raw,
            hvg_flavor: hvg_flavor
        };
        
        // Include n_pcs only if set (null means auto-select)
        if (n_pcs !== null) requestBody.n_pcs = n_pcs;
        
        // Only include optional parameters if they are set
        if (n_top_genes !== null) requestBody.n_top_genes = n_top_genes;
        if (min_genes !== null) requestBody.min_genes = min_genes;
        if (min_cells !== null) requestBody.min_cells = min_cells;
        if (min_counts !== null) requestBody.min_counts = min_counts;
        if (max_counts !== null) requestBody.max_counts = max_counts;
        if (max_genes !== null) requestBody.max_genes = max_genes;
        if (pct_mt_max !== null) requestBody.pct_mt_max = pct_mt_max;
        if (batch_key !== null) requestBody.batch_key = batch_key;
        if (regress_out_keys !== null && regress_out_keys.length > 0) requestBody.regress_out_keys = regress_out_keys;
        if (use_combat) requestBody.use_combat = use_combat;

        const response = await fetch(`${API_BASE}/preprocess`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(requestBody)
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Preprocessing failed');
        }

        showStatus('✓ Preprocessing complete!', 'success');
        
        // Display PCA variance ratio if available
        if (data.pca_variance_ratio && data.pca_variance_ratio.length > 0) {
            displayPCAVarianceRatio(data.pca_variance_ratio);
        }
        
        // Enable core analysis and cluster buttons (core analysis is independent)
        document.getElementById('btn-core-selection').disabled = false;
        document.getElementById('btn-cluster').disabled = false;
        
        // Update button states based on rules (visualization/marker genes still disabled until clustering/core analysis)
        updateButtonStates();
        
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
    
    // Disable all buttons during clustering
    disableAllButtons();
    btn.innerHTML = '<span class="spinner"></span> Clustering...';

    try {
        const method = document.getElementById('cluster-method').value;
        const resolution = parseFloat(document.getElementById('cluster-resolution').value) || 1.2;
        const use_rep = document.getElementById('use-rep').value || null;
        const key_added = null; // Use default (method name)
        
        // cplearn-specific parameters
        const stable_core_frac = parseFloat(document.getElementById('stable-core-frac')?.value) || 0.25;
        const stable_ng_num = parseInt(document.getElementById('stable-ng-num')?.value) || 8;
        const fine_grained = document.getElementById('fine-grained')?.checked || false;
        const propagate = document.getElementById('propagate')?.checked !== false;

        const requestBody = {
            session_id: sessionId,
            method: method,
            resolution: resolution
        };
        
        if (use_rep) requestBody.use_rep = use_rep;
        if (key_added) requestBody.key_added = key_added;
        
        // Add cplearn-specific parameters
        if (method === 'cplearn') {
            requestBody.stable_core_frac = stable_core_frac;
            requestBody.stable_ng_num = stable_ng_num;
            requestBody.fine_grained = fine_grained;
            requestBody.propagate = propagate;
        }

        const response = await fetch(`${API_BASE}/cluster`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(requestBody)
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Clustering failed');
        }

        showStatus(`✓ Clustering complete! Found ${data.n_clusters} clusters using ${data.cluster_key}.`, 'success');
        
        // Store cluster info
        // Preserve has_model from previous Core Analysis if it exists
        // This prevents coremap option from disappearing when running scanpy clustering
        const previousHasModel = clusterInfo && clusterInfo.has_model;
        clusterInfo = {
            cluster_key: data.cluster_key,
            clusters: data.clusters,
            cluster_counts: data.cluster_counts,
            has_model: data.has_model || previousHasModel || false
        };
        
        // If cplearn clustering, automatically run Core Analysis with same parameters (without visualization)
        if (method === 'cplearn' && data.has_model) {
            try {
                showStatus('Running Core Analysis with cplearn clustering parameters...', 'info');
                
                // Get parameters from Clustering form (same as used for clustering)
                const use_rep = document.getElementById('use-rep').value || null;
                const key_added = 'X_cplearn_coremap'; // Default coremap key
                const cluster_resolution = resolution; // Use same resolution as clustering
                const cluster_key = data.cluster_key; // Use the cluster key from clustering
                
                // Use same cplearn parameters from Clustering form
                const stable_core_frac = parseFloat(document.getElementById('stable-core-frac')?.value) || 0.25;
                const stable_ng_num = parseInt(document.getElementById('stable-ng-num')?.value) || 8;
                const fine_grained = document.getElementById('fine-grained')?.checked || false;
                const propagate = document.getElementById('propagate')?.checked !== false;
                
                // Call Core Analysis API with clustering parameters (without visualization)
                const coreResponse = await fetch(`${API_BASE}/core-analyze`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        session_id: sessionId,
                        use_rep: use_rep,
                        key_added: key_added,
                        cluster_resolution: cluster_resolution,
                        cluster_key: cluster_key,
                        stable_core_frac: stable_core_frac,
                        stable_ng_num: stable_ng_num,
                        fine_grained: fine_grained,
                        propagate: propagate,
                        visualize: false  // Don't generate visualization here, let Visualization module handle it
                    })
                });
                
                const coreData = await coreResponse.json();
                
                if (coreResponse.ok) {
                    // Update has_model flag
                    clusterInfo.has_model = true;
                    showStatus('✓ Core Analysis complete! Coremap embedding available for visualization.', 'success');
                } else {
                    console.warn('[CLUSTER] Core Analysis failed:', coreData.error || 'Unknown error');
                    showStatus('Clustering complete, but Core Analysis failed. You can run Core Analysis separately.', 'info');
                }
            } catch (coreError) {
                console.error('[CLUSTER] Error running Core Analysis after cplearn clustering:', coreError);
                // Don't fail the clustering operation if Core Analysis fails
                showStatus('Clustering complete, but Core Analysis failed. You can run Core Analysis separately.', 'info');
            }
        }
        
        // Update button states based on rules (now we have clustering, so visualization/marker genes can be enabled)
        updateButtonStates();
        
        // Update available cluster keys for DEG analysis (only existing ones)
        await updateAvailableClusterKeys();
        
        // Update visualization method selector (will show coremap if Core Analysis was run)
        updateVisualizationMethodSelector();
        
        // Update color-by select to include cluster key
        updateMetadataSelect([data.cluster_key]);
        document.getElementById('color-by').value = data.cluster_key;
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Clustering error:', error);
    } finally {
        // Restore button states
        btn.disabled = false;
        btn.textContent = originalText;
        
        // Re-enable buttons that should be available
        // Clustering doesn't prevent re-running preprocessing or core analysis
        // Users may want to compare different clustering parameters
        document.getElementById('btn-preprocess').disabled = false;
        document.getElementById('btn-core-selection').disabled = false;
        
        // Update button states based on rules (visualization/marker genes only if clustering/core analysis exists)
        updateButtonStates();
    }
}

async function runVisualization() {
    const params = getVisualizationParams();
    
    // For coremap, don't need cluster key
    if (params.method === 'coremap') {
        await runVisualizationForClusterKey(null);
        return;
    }
    
    // For UMAP, need cluster key
    const vizClusterKeySelect = document.getElementById('viz-cluster-key');
    const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
    
    if (selectedKey && selectedKey !== '') {
        // Use selected cluster key
        await runVisualizationForClusterKey(selectedKey);
    } else {
        // Use default cluster key
        const defaultKey = clusterInfo.cluster_key || 'cplearn';
        await runVisualizationForClusterKey(defaultKey);
    }
}

async function runVisualizationForClusterKey(cluster_key) {
    // Get current visualization parameters
    const params = getVisualizationParams();
    
    // Handle coremap method - generate visualization via API
    if (params.method === 'coremap') {
        // Check if coremap embedding exists (by checking if Core Analysis was run)
        if (!clusterInfo || !clusterInfo.has_model) {
            showStatus('Error: Core Analysis not available. Please run Core Analysis or cplearn clustering first.', 'error');
            return;
        }
        
    const btn = document.getElementById('btn-visualize');
    const originalText = btn.textContent;
        
        // Disable all buttons during visualization
        disableAllButtons();
        btn.innerHTML = '<span class="spinner"></span> Generating Coremap...';
        
        try {
            showStatus('Generating coremap visualization...', 'info');
            
            // Get Ground Truth file if provided
            const truthJsonFile = document.getElementById('viz-truth-json-file')?.files[0];
            let requestBody = {
                session_id: sessionId,
                cluster_key: cluster_key || 'cplearn',
                n_components: params.n_components,
                method: 'coremap',
                use_rep: 'X_cplearn_coremap'  // Default coremap key
            };
            
            let headers = { 'Content-Type': 'application/json' };
            let body;
            
            // If Ground Truth file is provided, use FormData
            if (truthJsonFile) {
                const formData = new FormData();
                formData.append('session_id', sessionId);
                formData.append('cluster_key', cluster_key || 'cplearn');
                formData.append('n_components', params.n_components.toString());
                formData.append('method', 'coremap');
                formData.append('use_rep', 'X_cplearn_coremap');
                formData.append('truth_json_file', truthJsonFile);
                
                headers = {}; // FormData sets Content-Type automatically
                body = formData;
            } else {
                body = JSON.stringify(requestBody);
            }
            
            // Call Visualization API to generate coremap visualization
        const response = await fetch(`${API_BASE}/visualize`, {
            method: 'POST',
                headers: headers,
                body: body
            });
            
            const data = await response.json();
            
            if (!response.ok) {
                throw new Error(data.error || 'Coremap visualization failed');
            }
            
            // Display the visualization in iframe
            const plotDiv = document.getElementById('plot');
            plotDiv.innerHTML = ''; // Clear existing content
            
            const iframe = document.createElement('iframe');
            iframe.style.width = '100%';
            iframe.style.height = '100%';
            iframe.style.border = 'none';
            iframe.srcdoc = data.html_content;
            plotDiv.appendChild(iframe);
            
            // Store visualization for future use
            if (data.file_url) {
                coreAnalysisVisualization = {
                    file_url: data.file_url,
                    htmlContent: data.html_content
                };
            }
            
            // Show visualization controls
            document.getElementById('visualization-section').style.display = 'block';
            document.getElementById('search-section').style.display = 'block';
            document.getElementById('stats-section').style.display = 'block';
            
            showStatus('✓ Coremap visualization generated!', 'success');
            
        } catch (error) {
            console.error('[VIZ] Error generating coremap visualization:', error);
            showStatus(`Error generating visualization: ${error.message}`, 'error');
        } finally {
            btn.disabled = false;
            btn.textContent = originalText;
            
            // Re-enable buttons that should be available
            document.getElementById('btn-marker-genes').disabled = false;
            document.getElementById('btn-cluster').disabled = false;
            document.getElementById('btn-core-selection').disabled = false;
            const btnPreprocess = document.getElementById('btn-preprocess');
            const pipelineSection = document.getElementById('pipeline-section');
            if (btnPreprocess && pipelineSection && pipelineSection.style.display !== 'none') {
                btnPreprocess.disabled = false;
            }
        }
        return;
    }
    
    // Handle standard visualization methods (umap, tsne, diffmap, draw_graph)
    const cacheKey = getVisualizationCacheKey(cluster_key, params.method, params.n_components, params.min_dist, params.spread);
    
    // Check if visualization is already cached with current parameters
    if (visualizationCache[cacheKey]) {
        vectors = visualizationCache[cacheKey].vectors;
        metadata = visualizationCache[cacheKey].metadata;
        updatePlot();
        showStatus(`✓ Loaded cached visualization for ${cluster_key}`, 'success');
        return;
    }
    
    const methodLabels = {
        'umap': 'UMAP',
        'tsne': 't-SNE',
        'diffmap': 'Diffusion Map',
        'draw_graph': 'Force-Directed Graph'
    };
    const methodLabel = methodLabels[params.method] || params.method.toUpperCase();
    
    showStatus(`Computing ${methodLabel} visualization for ${cluster_key}...`, 'info');
    const btn = document.getElementById('btn-visualize');
    const originalText = btn.textContent;
    
    // Disable all buttons during visualization
    disableAllButtons();
    btn.innerHTML = `<span class="spinner"></span> Computing ${methodLabel}...`;

    try {
        const requestBody = {
            session_id: sessionId,
            cluster_key: cluster_key,
            n_components: params.n_components,
            method: params.method
        };
        
        // Add method-specific parameters
        if (params.method === 'umap') {
            requestBody.min_dist = params.min_dist;
            requestBody.spread = params.spread;
        } else if (params.method === 'draw_graph') {
            requestBody.layout = 'fa'; // Default layout
        }
        
        const response = await fetch(`${API_BASE}/visualize`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(requestBody)
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Visualization failed');
        }

        // Store data in cache with parameters included
        visualizationCache[cacheKey] = {
            vectors: data.coordinates,
            metadata: data.metadata,
            params: params,
            method: params.method
        };
        
        // Update current vectors and metadata
        vectors = data.coordinates;
        metadata = data.metadata;
        
        // Update cluster key selector to show selected key
        const vizClusterKeySelect = document.getElementById('viz-cluster-key');
        if (vizClusterKeySelect) {
            vizClusterKeySelect.value = cluster_key;
        }
        
        // Automatically update Color By to match the cluster key
        const colorBySelect = document.getElementById('color-by');
        if (colorBySelect) {
            // Check if the cluster key exists in the color-by options
            const optionExists = Array.from(colorBySelect.options).some(opt => opt.value === cluster_key);
            if (!optionExists) {
                // Add the cluster key to color-by options if it doesn't exist
                const option = document.createElement('option');
                option.value = cluster_key;
                option.textContent = cluster_key;
                colorBySelect.appendChild(option);
            }
            // Set color-by to the cluster key
            colorBySelect.value = cluster_key;
        }
        
        showStatus(`✓ Visualization ready for ${cluster_key}!`, 'success');
        updatePlot();
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Visualization error:', error);
    } finally {
        // Restore button states
        btn.disabled = false;
        btn.textContent = originalText;
        
        // Re-enable buttons that should be available
        document.getElementById('btn-cluster').disabled = false;
        document.getElementById('btn-core-selection').disabled = false;
        // Preprocess button should remain enabled if data is loaded
        const btnPreprocess = document.getElementById('btn-preprocess');
        if (btnPreprocess && !btnPreprocess.disabled) {
            // If preprocess was enabled before, keep it enabled
        } else {
            // Otherwise, check if we have data (pipeline section is visible)
            const pipelineSection = document.getElementById('pipeline-section');
            if (pipelineSection && pipelineSection.style.display !== 'none') {
                btnPreprocess.disabled = false;
            }
        }
        
        // Update button states based on rules (visualization/marker genes only if clustering/core analysis exists)
        updateButtonStates();
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
            text: `Visualization (${vectors.length} points)`,
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

function displayPCAVarianceRatio(varianceRatio) {
    /**
     * Display PCA variance ratio plot in the visualization area
     */
    const plotDiv = document.getElementById('plot');
    if (!plotDiv) return;
    
    // Create PC numbers (1-indexed)
    const pcNumbers = varianceRatio.map((_, i) => i + 1);
    
    // Create trace for variance ratio
    const trace = {
        x: pcNumbers,
        y: varianceRatio,
        type: 'scatter',
        mode: 'lines+markers',
        name: 'Variance Ratio',
        line: {
            color: '#3b82f6',
            width: 2
        },
        marker: {
            color: '#3b82f6',
            size: 6
        }
    };
    
    // Calculate cumulative variance
    const cumulativeVariance = varianceRatio.reduce((acc, val, i) => {
        acc.push((i === 0 ? 0 : acc[i - 1]) + val);
        return acc;
    }, []);
    
    // Add cumulative variance trace
    const cumulativeTrace = {
        x: pcNumbers,
        y: cumulativeVariance,
        type: 'scatter',
        mode: 'lines',
        name: 'Cumulative Variance',
        yaxis: 'y2',
        line: {
            color: '#ef4444',
            width: 2,
            dash: 'dash'
        }
    };
    
    const layout = {
        title: {
            text: 'PCA Variance Ratio',
            font: { size: 18, color: '#1e293b' }
        },
        xaxis: {
            title: 'Principal Component',
            titlefont: { size: 14 },
            tickfont: { size: 12 }
        },
        yaxis: {
            title: 'Variance Ratio',
            titlefont: { size: 14 },
            tickfont: { size: 12 },
            type: 'log' // Log scale for better visualization
        },
        yaxis2: {
            title: 'Cumulative Variance',
            titlefont: { size: 14, color: '#ef4444' },
            tickfont: { size: 12, color: '#ef4444' },
            overlaying: 'y',
            side: 'right',
            range: [0, 1.1] // 0 to 110% for cumulative
        },
        hovermode: 'closest',
        showlegend: true,
        legend: {
            x: 0.7,
            y: 0.95
        },
        margin: { l: 60, r: 60, t: 60, b: 60 },
        plot_bgcolor: 'white',
        paper_bgcolor: 'white'
    };
    
    const config = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['lasso2d', 'select2d']
    };
    
    Plotly.newPlot('plot', [trace, cumulativeTrace], layout, config);
    
    // Show visualization section
    const visualizationSection = document.getElementById('visualization-section');
    if (visualizationSection) {
        visualizationSection.style.display = 'block';
    }
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
        showStatus('Running core analysis...', 'info');
    const btn = document.getElementById('btn-core-selection');
    const originalText = btn.textContent;
    
    // Disable all buttons during core analysis
    disableAllButtons();
    btn.innerHTML = '<span class="spinner"></span> Computing...';

    try {
        const use_rep = document.getElementById('core-use-rep')?.value || null;
        const cluster_resolution = parseFloat(document.getElementById('core-cluster-resolution')?.value) || 1.2;
        // Note: fast_view is now handled by Visualization module
        const key_added = 'X_cplearn_coremap'; // Use default
        const cluster_key = 'cplearn'; // Use default
        
        // cplearn clustering parameters (used when auto-running clustering)
        const stable_core_frac = parseFloat(document.getElementById('core-stable-core-frac')?.value) || 0.25;
        const stable_ng_num = parseInt(document.getElementById('core-stable-ng-num')?.value) || 8;
        const fine_grained = document.getElementById('core-fine-grained')?.checked || false;
        const propagate = document.getElementById('core-propagate')?.checked !== false;
        
        // Note: Ground Truth is now handled in Visualization module when using coremap
        
        // Prepare JSON request
        const requestBody = {
                session_id: sessionId,
                use_rep: use_rep,
                key_added: key_added,
                cluster_key: cluster_key,
                cluster_resolution: cluster_resolution,
            // Note: fast_view is now handled by Visualization module
                stable_core_frac: stable_core_frac,
                stable_ng_num: stable_ng_num,
                fine_grained: fine_grained,
                propagate: propagate
            };
            
        const headers = { 'Content-Type': 'application/json' };
        const body = JSON.stringify(requestBody);

        // Add timeout and retry logic for long-running operations
        const controller = new AbortController();
        const timeoutId = setTimeout(() => controller.abort(), 600000); // 10 minutes timeout
        
        let response;
        let data;
        let retryCount = 0;
        const maxRetries = 2;
        
        while (retryCount <= maxRetries) {
            try {
                response = await fetch(`${API_BASE}/core-analyze`, {
                    method: 'POST',
                    headers: headers,
                    body: body,
                    signal: controller.signal
                });
                
                clearTimeout(timeoutId);
                
                // Check if response is actually JSON
                const contentType = response.headers.get('content-type');
                const isJson = contentType && contentType.includes('application/json');
                
                // Handle non-JSON responses (502 Bad Gateway, 504 Gateway Timeout, etc.)
                if (!isJson) {
                    // Clone response to read text without consuming the body
                    const text = await response.clone().text();
                    console.error('[CORE] Non-JSON response:', text.substring(0, 200));
                    
                    if (response.status === 502 || response.status === 504) {
                        // Bad Gateway or Gateway Timeout - retry
                        if (retryCount < maxRetries) {
                            retryCount++;
                            const waitTime = Math.min(1000 * Math.pow(2, retryCount), 10000); // Exponential backoff, max 10s
                            showStatus(`Server timeout (${response.status}). Retrying in ${waitTime/1000}s... (attempt ${retryCount + 1}/${maxRetries + 1})`, 'info');
                            await new Promise(resolve => setTimeout(resolve, waitTime));
                            continue;
                        } else {
                            throw new Error(`Server error (${response.status}): The request took too long or the server is overloaded. Please try again later or use a smaller dataset.`);
                        }
                    } else {
                        throw new Error(`Server error (${response.status}): ${text.substring(0, 100)}`);
                    }
                }
                
                // Parse JSON response
                try {
                    data = await response.json();
                } catch (jsonError) {
                    // If JSON parsing fails, clone response to get text
                    const text = await response.clone().text();
                    console.error('[CORE] JSON parse error. Response text:', text.substring(0, 500));
                    throw new Error(`Invalid response from server: ${text.substring(0, 100)}`);
                }
                
                if (!response.ok) {
                    throw new Error(data.error || `Core analysis failed with status ${response.status}`);
                }
                
                // Success - break out of retry loop
                break;
                
            } catch (error) {
                clearTimeout(timeoutId);
                
                if (error.name === 'AbortError') {
                    throw new Error('Request timeout: The operation took longer than 10 minutes. Please try with a smaller dataset or contact support.');
                }
                
                // Network error or other fetch error
                if (retryCount < maxRetries && (error.message.includes('fetch') || error.message.includes('network'))) {
                    retryCount++;
                    const waitTime = Math.min(1000 * Math.pow(2, retryCount), 10000);
                    showStatus(`Network error. Retrying in ${waitTime/1000}s... (attempt ${retryCount + 1}/${maxRetries + 1})`, 'info');
                    await new Promise(resolve => setTimeout(resolve, waitTime));
                    continue;
                }
                
                // Re-throw if no more retries or non-retryable error
                throw error;
            }
        }

        let statusMessage = `✓ Core analysis complete! Embedding stored in ${data.key_added}`;
        if (data.assigned_points !== undefined) {
            statusMessage += ` (${data.assigned_points}/${data.total_points} points assigned`;
            if (data.core_cells !== undefined) {
                statusMessage += `, ${data.core_cells} core cells`;
            }
            statusMessage += ')';
        }
        
        showStatus(statusMessage, 'success');
        console.log('Core analysis result:', data);
        
        // Note: Visualization is now handled by Visualization module
        // No need to store visualization data here
        
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
            
            // Update color-by select to include cluster key and set it as selected
            if (data.cluster_key) {
                updateMetadataSelect([data.cluster_key]);
                document.getElementById('color-by').value = data.cluster_key;
            }
        }
        
        // Update button states based on rules (now we have core analysis, so visualization/marker genes can be enabled)
        updateButtonStates();
        
        // Update visualization method selector to show coremap option
        updateVisualizationMethodSelector();
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Core analysis error:', error);
    } finally {
        // Restore button states
        btn.disabled = false;
        btn.textContent = originalText;
        
        // Re-enable buttons that should be available
        // Core Analysis can run independently
        document.getElementById('btn-cluster').disabled = false;
        // Preprocess button should remain enabled if data is loaded
        const btnPreprocess = document.getElementById('btn-preprocess');
        if (btnPreprocess && !btnPreprocess.disabled) {
            // If preprocess was enabled before, keep it enabled
        } else {
            // Otherwise, check if we have data (pipeline section is visible)
            const pipelineSection = document.getElementById('pipeline-section');
            if (pipelineSection && pipelineSection.style.display !== 'none') {
                btnPreprocess.disabled = false;
            }
        }
        
        // Update button states based on rules (visualization/marker genes only if clustering/core analysis exists)
        updateButtonStates();
    }
}

async function runMarkerGenes() {
    showStatus('Finding marker genes...', 'info');
    const btn = document.getElementById('btn-marker-genes');
    const originalText = btn.textContent;
    
    // Disable all buttons during marker genes analysis
    disableAllButtons();
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
        // Restore button states
        btn.disabled = false;
        btn.textContent = originalText;
        
        // Re-enable buttons that should be available
        document.getElementById('btn-cluster').disabled = false;
        document.getElementById('btn-core-selection').disabled = false;
        // Preprocess button should remain enabled if data is loaded
        const btnPreprocess = document.getElementById('btn-preprocess');
        if (btnPreprocess && !btnPreprocess.disabled) {
            // If preprocess was enabled before, keep it enabled
        } else {
            // Otherwise, check if we have data (pipeline section is visible)
            const pipelineSection = document.getElementById('pipeline-section');
            if (pipelineSection && pipelineSection.style.display !== 'none') {
                btnPreprocess.disabled = false;
            }
        }
        
        // Update button states based on rules (visualization/marker genes only if clustering/core analysis exists)
        updateButtonStates();
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
    // Create toast container if it doesn't exist
    let toastContainer = document.getElementById('toast-container');
    if (!toastContainer) {
        toastContainer = document.createElement('div');
        toastContainer.id = 'toast-container';
        toastContainer.className = 'toast-container';
        document.body.appendChild(toastContainer);
    }

    // Create toast element
    const toast = document.createElement('div');
    toast.className = `toast ${type}`;

    // Add icon based on type
    const icon = document.createElement('span');
    icon.className = 'toast-icon';
    if (type === 'success') {
        icon.textContent = '✓';
    } else if (type === 'error') {
        icon.textContent = '✗';
    } else {
        icon.textContent = 'ℹ';
    }

    // Add message
    const messageSpan = document.createElement('span');
    messageSpan.className = 'toast-message';
    messageSpan.textContent = message;

    // Add close button
    const closeBtn = document.createElement('button');
    closeBtn.className = 'toast-close';
    closeBtn.innerHTML = '×';
    closeBtn.onclick = () => removeToast(toast);

    // Assemble toast
    toast.appendChild(icon);
    toast.appendChild(messageSpan);
    toast.appendChild(closeBtn);

    // Add to container
    toastContainer.appendChild(toast);

    // Auto-remove after delay (errors don't auto-remove, only success/info)
    if (type !== 'error') {
        const delay = type === 'info' ? 5000 : 4000;
        setTimeout(() => {
            removeToast(toast);
        }, delay);
    }
}

function removeToast(toast) {
    if (toast && toast.parentNode) {
        toast.classList.add('slide-out');
        setTimeout(() => {
            if (toast.parentNode) {
                toast.parentNode.removeChild(toast);
            }
        }, 300);
    }
}
=======
// Lotus - Main Application Logic with Backend Integration

const API_BASE = '/api';
let vectors = null;
let metadata = null;
let currentPlot = null;
let colorMapping = {};
// Generate unique session ID for each user
let sessionId = generateSessionId();
let clusterInfo = {}; // Store cluster information
let heartbeatInterval = null; // Heartbeat interval ID
let sessionCleanupSent = false; // Flag to prevent multiple cleanup requests
let visualizationCache = {}; // Store visualization results: { cache_key: { vectors, metadata, params } }
let coreAnalysisVisualization = null; // Store Core Analysis visualization data: { file_url, htmlContent, iframe }
// Cache key format: `${cluster_key}_${n_components}_${min_dist}_${spread}`

// Function to generate cache key from cluster key, method, and visualization parameters
function getVisualizationCacheKey(cluster_key, method, n_components, min_dist, spread) {
    if (method === 'coremap') {
        return `${cluster_key}_coremap_${n_components}`;
    } else if (method === 'umap') {
        return `${cluster_key}_umap_${n_components}_${min_dist}_${spread}`;
    } else {
        // For tsne, diffmap, draw_graph - use method and n_components
        return `${cluster_key}_${method}_${n_components}`;
    }
}

// Function to get current visualization parameters
function getVisualizationParams() {
    const method = document.getElementById('viz-method')?.value || 'umap';
    const n_components = document.getElementById('dimension').value === '3d' ? 3 : 2;
    const min_dist = parseFloat(document.getElementById('umap-min-dist').value) || 0.5;
    const spread = parseFloat(document.getElementById('umap-spread').value) || 1.0;
    return { method, n_components, min_dist, spread };
}

// Function to update visualization method selector based on Core Analysis availability
function updateVisualizationMethodSelector() {
    const methodSelect = document.getElementById('viz-method');
    if (!methodSelect) return;
    
    // Check if Core Analysis has been run
    // Use clusterInfo.has_model which is set when Core Analysis or cplearn clustering is run
    // This is preserved even when running scanpy clustering after Core Analysis
    const hasCoreAnalysis = clusterInfo && clusterInfo.has_model;
    
    // Get current value
    const currentValue = methodSelect.value;
    
    // Clear options
    methodSelect.innerHTML = '';
    
    // Always add UMAP
    const umapOption = document.createElement('option');
    umapOption.value = 'umap';
    umapOption.textContent = 'UMAP';
    methodSelect.appendChild(umapOption);
    
    // Add Coremap if Core Analysis is available
    if (hasCoreAnalysis) {
        const coremapOption = document.createElement('option');
        coremapOption.value = 'coremap';
        coremapOption.textContent = 'Coremap';
        methodSelect.appendChild(coremapOption);
    }
    
    // Restore previous selection if still valid, otherwise default to UMAP
    if (currentValue && Array.from(methodSelect.options).some(opt => opt.value === currentValue)) {
        methodSelect.value = currentValue;
    } else {
        methodSelect.value = 'umap';
    }
    
    // Update parameters visibility
    updateUMAPParamsVisibility();
}

// Function to update parameters visibility based on selected method
function updateUMAPParamsVisibility() {
    const methodSelect = document.getElementById('viz-method');
    const umapParamsGroup = document.getElementById('umap-params-group');
    const umapSpreadGroup = document.getElementById('umap-spread-group');
    const vizClusterKeyGroup = document.getElementById('viz-cluster-key')?.closest('.form-group');
    const dimensionGroup = document.getElementById('dimension')?.closest('.form-group');
    const coremapGroundTruthGroup = document.getElementById('coremap-ground-truth-group');
    const coremapFastViewGroup = document.getElementById('coremap-fast-view-group');
    
    if (methodSelect) {
        const method = methodSelect.value;
        if (method === 'coremap') {
            // Hide UMAP parameters, Cluster Key, and Dimension for coremap
            if (umapParamsGroup) umapParamsGroup.style.display = 'none';
            if (umapSpreadGroup) umapSpreadGroup.style.display = 'none';
            if (vizClusterKeyGroup) vizClusterKeyGroup.style.display = 'none';
            if (dimensionGroup) dimensionGroup.style.display = 'none';
            // Show Fast View and Ground Truth upload for coremap
            if (coremapFastViewGroup) coremapFastViewGroup.style.display = 'block';
            if (coremapGroundTruthGroup) coremapGroundTruthGroup.style.display = 'block';
        } else {
            // Show parameters for standard methods (umap, tsne, diffmap, draw_graph)
            // UMAP-specific parameters only shown for UMAP
            if (method === 'umap') {
                if (umapParamsGroup) umapParamsGroup.style.display = 'block';
                if (umapSpreadGroup) umapSpreadGroup.style.display = 'block';
            } else {
                // Hide UMAP-specific parameters for other methods
                if (umapParamsGroup) umapParamsGroup.style.display = 'none';
                if (umapSpreadGroup) umapSpreadGroup.style.display = 'none';
            }
            // Show Cluster Key and Dimension for all standard methods
            if (vizClusterKeyGroup) vizClusterKeyGroup.style.display = 'block';
            if (dimensionGroup) dimensionGroup.style.display = 'block';
            // Hide Fast View and Ground Truth upload for standard methods
            if (coremapFastViewGroup) coremapFastViewGroup.style.display = 'none';
            if (coremapGroundTruthGroup) coremapGroundTruthGroup.style.display = 'none';
        }
    }
}

// Generate unique session ID using UUID v4
function generateSessionId() {
    // Use crypto.randomUUID if available (modern browsers)
    if (typeof crypto !== 'undefined' && crypto.randomUUID) {
        return crypto.randomUUID();
    }
    // Fallback: generate UUID v4 manually
    return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
        const r = Math.random() * 16 | 0;
        const v = c === 'x' ? r : (r & 0x3 | 0x8);
        return v.toString(16);
    });
}

// Initialize session on page load
async function initializeSession() {
    try {
        const response = await fetch(`${API_BASE}/session/create`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId })
        });
        if (response.ok) {
            console.log(`[SESSION] Session initialized: ${sessionId}`);
            // Start heartbeat
            startHeartbeat();
        } else {
            console.warn('[SESSION] Failed to initialize session');
        }
    } catch (error) {
        console.error('[SESSION] Error initializing session:', error);
    }
}

// Send heartbeat to keep session alive
async function sendHeartbeat() {
    try {
        await fetch(`${API_BASE}/session/heartbeat`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId })
        });
    } catch (error) {
        console.error('[SESSION] Heartbeat error:', error);
    }
}

// Start heartbeat (send every 30 seconds)
function startHeartbeat() {
    // Send initial heartbeat
    sendHeartbeat();
    // Set up interval
    heartbeatInterval = setInterval(sendHeartbeat, 30000); // 30 seconds
}

// Stop heartbeat
function stopHeartbeat() {
    if (heartbeatInterval) {
        clearInterval(heartbeatInterval);
        heartbeatInterval = null;
    }
}

// Clean up session when page is closed or hidden
async function cleanupSession() {
    if (sessionCleanupSent) return; // Prevent multiple cleanup requests
    sessionCleanupSent = true;
    
    stopHeartbeat();
    
    try {
        // Use sendBeacon for reliable cleanup on page unload
        const data = JSON.stringify({ session_id: sessionId });
        if (navigator.sendBeacon) {
            navigator.sendBeacon(`${API_BASE}/session/cleanup`, data);
            console.log(`[SESSION] Session cleanup sent via beacon: ${sessionId}`);
        } else {
            // Fallback: use fetch with keepalive
            await fetch(`${API_BASE}/session/cleanup`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: data,
                keepalive: true
            });
            console.log(`[SESSION] Session cleanup sent: ${sessionId}`);
        }
    } catch (error) {
        console.error('[SESSION] Error cleaning up session:', error);
    }
}

// Set up cleanup handlers
window.addEventListener('beforeunload', cleanupSession);
document.addEventListener('visibilitychange', function() {
    if (document.visibilityState === 'hidden') {
        // Page is hidden, but don't cleanup yet (user might come back)
        // Just stop heartbeat to save resources
        stopHeartbeat();
    } else if (document.visibilityState === 'visible') {
        // Page is visible again, restart heartbeat
        if (!heartbeatInterval) {
            startHeartbeat();
        }
    }
});

// Initialize session when page loads
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initializeSession);
} else {
    initializeSession();
}

// Function to disable all buttons during operations
// Note: File upload and default dataset buttons are NOT disabled to allow reloading data at any time
function disableAllButtons() {
    const buttons = [
        'btn-preprocess',
        'btn-cluster',
        'btn-core-selection',
        'btn-visualize',
        'btn-marker-genes'
    ];
    buttons.forEach(btnId => {
        const btn = document.getElementById(btnId);
        if (btn) {
            btn.disabled = true;
        }
    });
}

// Function to disable data loading buttons during loading
function disableDataLoadingButtons() {
    const buttons = [
        'load-demo-data',
        'load-pbmc3k',
        'data-file'  // File input (though it's hidden, we can disable it)
    ];
    buttons.forEach(btnId => {
        const btn = document.getElementById(btnId);
        if (btn) {
            btn.disabled = true;
        }
    });
}

// Function to enable data loading buttons
function enableDataLoadingButtons() {
    const buttons = [
        'load-demo-data',
        'load-pbmc3k',
        'data-file'
    ];
    buttons.forEach(btnId => {
        const btn = document.getElementById(btnId);
        if (btn) {
            btn.disabled = false;
        }
    });
}

// Function to update button states based on current data state
// Rule 2: Visualization and Marker Genes require either Core Analysis or Clustering to be completed
function updateButtonStates() {
    // Check if we have clustering or core analysis results
    const hasClustering = clusterInfo && clusterInfo.cluster_key;
    const hasCoreAnalysis = clusterInfo && clusterInfo.has_model;
    
    // Rule 2: Visualization and Marker Genes require clustering or core analysis
    const canVisualize = hasClustering || hasCoreAnalysis;
    
    const btnVisualize = document.getElementById('btn-visualize');
    const btnMarkerGenes = document.getElementById('btn-marker-genes');
    
    if (btnVisualize) {
        btnVisualize.disabled = !canVisualize;
    }
    if (btnMarkerGenes) {
        btnMarkerGenes.disabled = !canVisualize;
    }
    
    // Update cluster key select for Marker Genes
    const degClusterKeySelect = document.getElementById('deg-cluster-key');
    if (degClusterKeySelect && canVisualize) {
        // If we can visualize, the select should be enabled if cluster keys are available
        // This will be handled by updateAvailableClusterKeys
    } else if (degClusterKeySelect && !canVisualize) {
        degClusterKeySelect.disabled = true;
    }
}

// Function to reset all state when new data is loaded
function resetAllState() {
    console.log('[RESET] Resetting all state...');
    
    // Clear cluster info
    clusterInfo = {};
    
    // Clear visualization cache
    visualizationCache = {};
    
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
    
    // Reset cluster keys dropdowns
    const degClusterKeySelect = document.getElementById('deg-cluster-key');
    if (degClusterKeySelect) {
        degClusterKeySelect.innerHTML = '<option value="">No cluster keys available - Run clustering first</option>';
        degClusterKeySelect.disabled = true;
    }
    const vizClusterKeySelect = document.getElementById('viz-cluster-key');
    if (vizClusterKeySelect) {
        vizClusterKeySelect.innerHTML = '<option value="">No cluster keys available - Run clustering first</option>';
        vizClusterKeySelect.disabled = true;
    }
    
    // DEG step is always visible (same as visualization)
    
    // Clear coremap visualization if exists
    const coremapViz = document.getElementById('coremap-visualization');
    if (coremapViz) {
        coremapViz.innerHTML = '';
        coremapViz.style.display = 'none';
    }
    
    // Clear Core Analysis visualization data (if exists, for caching)
    coreAnalysisVisualization = null;
    
    // Clear Ground Truth file
    clearVizGroundTruthFile();
    
    console.log('[RESET] All state reset complete');
    
    // Update button states based on rules
    updateButtonStates();
}

// Function to update available cluster keys dropdown
async function updateAvailableClusterKeys() {
    try {
        const response = await fetch(`${API_BASE}/cluster-keys`, {
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

        // Update both cluster key selectors (DEG and Visualization)
        const degSelect = document.getElementById('deg-cluster-key');
        const vizSelect = document.getElementById('viz-cluster-key');

        // Update DEG cluster key selector
        if (degSelect) {
            degSelect.innerHTML = '';
        if (data.cluster_keys && data.cluster_keys.length > 0) {
            data.cluster_keys.forEach(clusterKeyInfo => {
                const option = document.createElement('option');
                option.value = clusterKeyInfo.key;
                const label = `${clusterKeyInfo.key} (${clusterKeyInfo.method}, ${clusterKeyInfo.n_clusters} clusters)`;
                option.textContent = label;
                option.selected = clusterKeyInfo.key === clusterInfo.cluster_key;
                    degSelect.appendChild(option);
                });
                degSelect.disabled = false;
            } else {
                const option = document.createElement('option');
                option.value = '';
                option.textContent = 'No cluster keys available - Run clustering first';
                option.disabled = true;
                degSelect.appendChild(option);
                degSelect.disabled = true;
            }
        }

        // Update Visualization cluster key selector
        if (vizSelect) {
            vizSelect.innerHTML = '';
            if (data.cluster_keys && data.cluster_keys.length > 0) {
                data.cluster_keys.forEach(clusterKeyInfo => {
                    const option = document.createElement('option');
                    option.value = clusterKeyInfo.key;
                    const label = `${clusterKeyInfo.key} (${clusterKeyInfo.method}, ${clusterKeyInfo.n_clusters} clusters)`;
                    option.textContent = label;
                    option.selected = clusterKeyInfo.key === clusterInfo.cluster_key;
                    vizSelect.appendChild(option);
                });
                vizSelect.disabled = false;
                // Do not auto-trigger visualization when updating cluster keys
                // User should manually select cluster key or click "Generate Visualization" button
        } else {
            const option = document.createElement('option');
            option.value = '';
            option.textContent = 'No cluster keys available - Run clustering first';
            option.disabled = true;
                vizSelect.appendChild(option);
                vizSelect.disabled = true;
            }
        }
        
        // Update button states based on rules after cluster keys are updated
        updateButtonStates();
        
        // Update visualization method selector (to show/hide coremap option)
        updateVisualizationMethodSelector();
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
    try {
        console.log('[INIT] Initializing application...');
    checkHealth();
    setupFileInputs();
        setupDefaultDatasets();
    setupControls();
    setupPipelineButtons();
    // Initialize cluster keys dropdown (will be empty until clustering is run)
    updateAvailableClusterKeys();
        console.log('[INIT] Application initialized successfully');
    } catch (error) {
        console.error('[INIT] Error during initialization:', error);
        showStatus('Error initializing application. Please refresh the page.', 'error');
    }
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

// Store MTX files temporarily
let mtxFiles = {
    matrix: null,
    features: null,
    barcodes: null
};

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
    
    // Setup MTX file inputs - click buttons to upload
    const matrixBtn = document.getElementById('mtx-matrix-btn');
    const featuresBtn = document.getElementById('mtx-features-btn');
    const barcodesBtn = document.getElementById('mtx-barcodes-btn');
    const matrixInput = document.getElementById('mtx-matrix-file');
    const featuresInput = document.getElementById('mtx-features-file');
    const barcodesInput = document.getElementById('mtx-barcodes-file');
    
    if (matrixBtn && matrixInput) {
        matrixBtn.addEventListener('click', () => matrixInput.click());
        matrixInput.addEventListener('change', (e) => handleMtxFileUpload(e, 'matrix'));
    }
    if (featuresBtn && featuresInput) {
        featuresBtn.addEventListener('click', () => featuresInput.click());
        featuresInput.addEventListener('change', (e) => handleMtxFileUpload(e, 'features'));
    }
    if (barcodesBtn && barcodesInput) {
        barcodesBtn.addEventListener('click', () => barcodesInput.click());
        barcodesInput.addEventListener('change', (e) => handleMtxFileUpload(e, 'barcodes'));
    }
    
    // Setup Ground Truth file upload for Visualization (coremap)
    setupGroundTruthFileInput();
}

function setupGroundTruthFileInput() {
    try {
        // Setup for Visualization Ground Truth upload (coremap)
        const vizGroundTruthInput = document.getElementById('viz-truth-json-file');
        const vizGroundTruthUploadArea = document.getElementById('viz-ground-truth-upload-area');
        
        if (!vizGroundTruthInput || !vizGroundTruthUploadArea) {
            // Elements don't exist yet (might be hidden), that's okay
            return;
        }
        
        // Drag and drop for Ground Truth
        vizGroundTruthUploadArea.addEventListener('dragover', (e) => {
            e.preventDefault();
            e.stopPropagation();
            vizGroundTruthUploadArea.classList.add('dragover');
        });
        
        vizGroundTruthUploadArea.addEventListener('dragleave', () => {
            vizGroundTruthUploadArea.classList.remove('dragover');
        });
        
        vizGroundTruthUploadArea.addEventListener('drop', (e) => {
            e.preventDefault();
            e.stopPropagation();
            vizGroundTruthUploadArea.classList.remove('dragover');
            
            const files = e.dataTransfer.files;
            if (files.length > 0 && files[0].type === 'application/json') {
                vizGroundTruthInput.files = files;
                
                // Update UI to show selected file
                updateVizGroundTruthFileDisplay(files[0].name);
            }
        });
        
        vizGroundTruthUploadArea.addEventListener('click', () => {
            vizGroundTruthInput.click();
        });
        
        // Update UI when file is selected via click
        vizGroundTruthInput.addEventListener('change', (e) => {
            if (e.target.files && e.target.files.length > 0) {
                updateVizGroundTruthFileDisplay(e.target.files[0].name);
            }
        });
        
        // Setup remove button
        const removeBtn = document.getElementById('viz-ground-truth-remove-btn');
        if (removeBtn) {
            removeBtn.addEventListener('click', (e) => {
                e.stopPropagation(); // Prevent triggering file input click
                clearVizGroundTruthFile();
            });
        }
    } catch (error) {
        console.error('[SETUP] Error setting up Ground Truth file input:', error);
    }
}

function updateVizGroundTruthFileDisplay(filename) {
    const uploadArea = document.getElementById('viz-ground-truth-upload-area');
    if (!uploadArea) return;
    
    const textElement = uploadArea.querySelector('.file-upload-text');
    const hintElement = uploadArea.querySelector('.file-upload-hint');
    const removeBtn = document.getElementById('viz-ground-truth-remove-btn');
    
    if (textElement) {
        textElement.textContent = filename;
        textElement.style.color = 'var(--primary)';
    }
    if (hintElement) {
        hintElement.textContent = 'Click to change file';
    }
    if (removeBtn) {
        removeBtn.style.display = 'block';
    }
}

function clearVizGroundTruthFile() {
    const vizGroundTruthInput = document.getElementById('viz-truth-json-file');
    const uploadArea = document.getElementById('viz-ground-truth-upload-area');
    if (!uploadArea || !vizGroundTruthInput) return;
    
    // Clear file input
    vizGroundTruthInput.value = '';
    
    // Reset display
    const textElement = uploadArea.querySelector('.file-upload-text');
    const hintElement = uploadArea.querySelector('.file-upload-hint');
    const removeBtn = document.getElementById('viz-ground-truth-remove-btn');
    
    if (textElement) {
        textElement.textContent = 'Click or drag JSON file';
        textElement.style.color = '';
    }
    if (hintElement) {
        hintElement.textContent = '.json';
    }
    if (removeBtn) {
        removeBtn.style.display = 'none';
    }
}

function setupDefaultDatasets() {
    const loadDemoDataBtn = document.getElementById('load-demo-data');
    const loadPbmc3kBtn = document.getElementById('load-pbmc3k');
    
    if (loadDemoDataBtn) {
        loadDemoDataBtn.addEventListener('click', () => {
            loadDefaultDataset('demo_data.h5ad');
        });
    }
    
    if (loadPbmc3kBtn) {
        loadPbmc3kBtn.addEventListener('click', () => {
            loadDefaultDataset('pbmc3k_raw.h5ad');
        });
    }
}

async function loadDefaultDataset(filename) {
    console.log(`[LOAD] Loading default dataset: ${filename}`);
    showStatus(`Loading ${filename}...`, 'info');
    
    // Reset all state when loading new data
    resetAllState();
    
    // Disable pipeline buttons and data loading buttons during loading
    disableAllButtons();
    disableDataLoadingButtons();
    
    try {
        const response = await fetch(`${API_BASE}/load-default-dataset`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                filename: filename,
                session_id: sessionId
            })
        });

        let data;
        try {
            data = await response.json();
            console.log('[LOAD] Response data:', data);
        } catch (e) {
            const text = await response.text();
            console.error('[LOAD] Failed to parse JSON:', text);
            throw new Error(`Server error: ${response.status} ${response.statusText}. Response: ${text.substring(0, 200)}`);
        }

        if (!response.ok) {
            console.error('[LOAD] Error response:', data);
            throw new Error(data.error || `Load failed: ${response.status}`);
        }

        console.log('[LOAD] Success:', data);
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
        
        // Update available cluster keys (in case loaded data already has clustering results)
        await updateAvailableClusterKeys();
        
        // Update button states based on rules (check if loaded data has clustering/core analysis)
        updateButtonStates();
        
        // Update visualization method selector
        updateVisualizationMethodSelector();
        
    } catch (error) {
        console.error('[LOAD] Exception:', error);
        showStatus(`Error: ${error.message}`, 'error');
    } finally {
        // Always re-enable data loading buttons to allow reloading at any time
        enableDataLoadingButtons();
    }
}

function setupControls() {
    // Initialize visualization method selector
    updateVisualizationMethodSelector();
    
    // Add event listener for visualization method selector
    const vizMethodSelect = document.getElementById('viz-method');
    if (vizMethodSelect) {
        vizMethodSelect.addEventListener('change', (e) => {
            // Update UMAP parameters visibility
            updateUMAPParamsVisibility();
            
            // Clear cache for current cluster key when method changes
            const vizClusterKeySelect = document.getElementById('viz-cluster-key');
            const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
            if (selectedKey && selectedKey !== '') {
                Object.keys(visualizationCache).forEach(key => {
                    if (key.startsWith(selectedKey + '_')) {
                        delete visualizationCache[key];
                    }
                });
            }
            
            // Don't auto-generate visualization when switching to coremap
            // User needs to click "Generate Visualization" button
        });
    }
    
    document.getElementById('color-by').addEventListener('change', updatePlot);
    document.getElementById('dimension').addEventListener('change', (e) => {
        // When dimension changes, clear current visualization cache for the selected cluster key
        // This forces regeneration with new parameters
        const vizClusterKeySelect = document.getElementById('viz-cluster-key');
        const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
        if (selectedKey && selectedKey !== '') {
            // Clear cache entries for this cluster key (they will be regenerated with new parameters)
            Object.keys(visualizationCache).forEach(key => {
                if (key.startsWith(selectedKey + '_')) {
                    delete visualizationCache[key];
                }
            });
        }
        updatePlot();
    });
    
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
    
    // Update cplearn advanced parameters visibility based on method selection
    const clusterMethodSelect = document.getElementById('cluster-method');
    const cplearnAdvancedParams = document.getElementById('cplearn-advanced-params');
    
    function updateCplearnParamsVisibility() {
        if (clusterMethodSelect && cplearnAdvancedParams) {
            const method = clusterMethodSelect.value;
            if (method === 'cplearn') {
                cplearnAdvancedParams.style.display = 'block';
            } else {
                cplearnAdvancedParams.style.display = 'none';
            }
        }
    }
    
    if (clusterMethodSelect) {
        clusterMethodSelect.addEventListener('change', updateCplearnParamsVisibility);
        // Initialize visibility on page load
        updateCplearnParamsVisibility();
    }
    
    // Add event listeners for visualization parameters to clear cache when changed
    const umapMinDistInput = document.getElementById('umap-min-dist');
    const umapSpreadInput = document.getElementById('umap-spread');
    
    if (umapMinDistInput) {
        umapMinDistInput.addEventListener('change', () => {
            // Clear cache for current cluster key when min_dist changes
            const vizClusterKeySelect = document.getElementById('viz-cluster-key');
            const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
            if (selectedKey && selectedKey !== '') {
                Object.keys(visualizationCache).forEach(key => {
                    if (key.startsWith(selectedKey + '_')) {
                        delete visualizationCache[key];
                    }
                });
            }
        });
    }
    
    if (umapSpreadInput) {
        umapSpreadInput.addEventListener('change', () => {
            // Clear cache for current cluster key when spread changes
            const vizClusterKeySelect = document.getElementById('viz-cluster-key');
            const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
            if (selectedKey && selectedKey !== '') {
                Object.keys(visualizationCache).forEach(key => {
                    if (key.startsWith(selectedKey + '_')) {
                        delete visualizationCache[key];
                    }
                });
            }
        });
    }
    
    // Add event listener for Visualization cluster key selector
    const vizClusterKeySelect = document.getElementById('viz-cluster-key');
    if (vizClusterKeySelect) {
        vizClusterKeySelect.addEventListener('change', async (e) => {
            const selectedKey = e.target.value;
            if (selectedKey && selectedKey !== '') {
                // Automatically update Color By to match the selected cluster key
                const colorBySelect = document.getElementById('color-by');
                if (colorBySelect) {
                    // Check if the cluster key exists in the color-by options
                    const optionExists = Array.from(colorBySelect.options).some(opt => opt.value === selectedKey);
                    if (!optionExists) {
                        // Add the cluster key to color-by options if it doesn't exist
                        const option = document.createElement('option');
                        option.value = selectedKey;
                        option.textContent = selectedKey;
                        colorBySelect.appendChild(option);
                    }
                    // Set color-by to the selected cluster key
                    colorBySelect.value = selectedKey;
                }
                
                // Automatically generate new visualization when switching cluster key
                // Cache will be checked inside runVisualizationForClusterKey with current parameters
                await runVisualizationForClusterKey(selectedKey);
            } else {
                // Clear visualization if no cluster key selected
                vectors = null;
                metadata = null;
                const plotDiv = document.getElementById('plot');
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #64748b;">Select a cluster key to visualize</div>';
                }
            }
        });
    }
}

function setupPipelineButtons() {
    try {
        const btnPreprocess = document.getElementById('btn-preprocess');
        const btnCluster = document.getElementById('btn-cluster');
        const btnVisualize = document.getElementById('btn-visualize');
        const btnCoreSelection = document.getElementById('btn-core-selection');
        const btnMarkerGenes = document.getElementById('btn-marker-genes');
        
        if (btnPreprocess) {
            btnPreprocess.addEventListener('click', runPreprocess);
        } else {
            console.error('[SETUP] btn-preprocess not found');
        }
        
        if (btnCluster) {
            btnCluster.addEventListener('click', runClustering);
        } else {
            console.error('[SETUP] btn-cluster not found');
        }
        
        if (btnVisualize) {
            btnVisualize.addEventListener('click', runVisualization);
        } else {
            console.error('[SETUP] btn-visualize not found');
        }
        
        if (btnCoreSelection) {
            btnCoreSelection.addEventListener('click', runCoreSelection);
        } else {
            console.error('[SETUP] btn-core-selection not found');
        }
        
        if (btnMarkerGenes) {
            btnMarkerGenes.addEventListener('click', runMarkerGenes);
        } else {
            console.error('[SETUP] btn-marker-genes not found');
        }
    } catch (error) {
        console.error('[SETUP] Error setting up pipeline buttons:', error);
    }
}

async function handleDataFile(event) {
    const file = event.target.files[0];
    if (!file) {
        console.log('[UPLOAD] No file selected');
        return;
    }

    console.log('[UPLOAD] File selected:', file.name, 'Size:', file.size, 'bytes');
    showStatus('Uploading data...', 'info');
    
    // Reset all state when loading new data
    resetAllState();
    
    // Disable pipeline buttons and data loading buttons during upload
    disableAllButtons();
    disableDataLoadingButtons();
    
    // Determine file type from extension (h5ad, csv, tsv allowed; mtx uses separate upload)
    const fileName = file.name.toLowerCase();
    let fileType = null;
    
    if (fileName.endsWith('.h5ad')) {
        fileType = 'h5ad';
    } else if (fileName.endsWith('.csv')) {
        fileType = 'csv';
    } else if (fileName.endsWith('.tsv')) {
        fileType = 'tsv';
    } else {
        showStatus('Unsupported file format. Please upload .h5ad, .csv, or .tsv files only. For MTX format, use the separate file upload option below.', 'error');
        enableDataLoadingButtons(); // Re-enable buttons if file format is invalid
        return;
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
        
        // Update button states based on rules (check if uploaded data has clustering/core analysis)
        updateButtonStates();
        
        // Update visualization method selector
        updateVisualizationMethodSelector();
        
    } catch (error) {
        console.error('[UPLOAD] Exception:', error);
        showStatus(`Error: ${error.message}`, 'error');
    } finally {
        // Always re-enable data loading buttons to allow reloading at any time
        enableDataLoadingButtons();
    }
}

// Handle individual MTX file upload
async function handleMtxFileUpload(event, fileType) {
    const file = event.target.files[0];
    if (!file) {
        return;
    }
    
    console.log(`[UPLOAD] MTX ${fileType} file selected:`, file.name);
    
    // Update UI to show file is selected
    const btnId = `mtx-${fileType}-btn`;
    const statusId = `mtx-${fileType}-status`;
    const btn = document.getElementById(btnId);
    const status = document.getElementById(statusId);
    
    if (btn) {
        btn.classList.add('uploaded');
        btn.innerHTML = `✓ ${fileType.charAt(0).toUpperCase() + fileType.slice(1)}:`;
    }
    if (status) {
        status.textContent = file.name;
        status.style.color = 'var(--success)';
    }
    
    // Store file reference
    mtxFiles[fileType] = file;
    
    // Check if all three files are ready
    if (mtxFiles.matrix && mtxFiles.features && mtxFiles.barcodes) {
        // All files ready, combine and load
        showStatus('All files ready, loading data...', 'info');
        await combineAndLoadMtxFiles();
    } else {
        // Show which files are still needed
        const missing = [];
        if (!mtxFiles.matrix) missing.push('Matrix');
        if (!mtxFiles.features) missing.push('Features');
        if (!mtxFiles.barcodes) missing.push('Barcodes');
        showStatus(`✓ ${fileType} uploaded. Still need: ${missing.join(', ')}`, 'info');
    }
}

// Combine and load all MTX files
async function combineAndLoadMtxFiles() {
    if (!mtxFiles.matrix || !mtxFiles.features || !mtxFiles.barcodes) {
        showStatus('Please upload all three files: Matrix, Features, and Barcodes', 'error');
        return;
    }
    
    console.log('[UPLOAD] All MTX files ready, combining and loading...', {
        matrix: mtxFiles.matrix.name,
        features: mtxFiles.features.name,
        barcodes: mtxFiles.barcodes.name
    });
    
    showStatus('Combining and loading MTX files...', 'info');
    
    // Reset all state when loading new data
    resetAllState();
    
    // Disable pipeline buttons and data loading buttons during upload
    disableAllButtons();
    disableDataLoadingButtons();
    
    try {
        const formData = new FormData();
        formData.append('matrix_file', mtxFiles.matrix);
        formData.append('features_file', mtxFiles.features);
        formData.append('barcodes_file', mtxFiles.barcodes);
        formData.append('type', 'mtx_multi');
        formData.append('session_id', sessionId);

        console.log('[UPLOAD] Sending MTX files to:', `${API_BASE}/upload`);
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
        
        // Show warning if Seurat format was detected
        if (data.warning) {
            setTimeout(() => {
                showStatus(`⚠ ${data.warning}`, 'warning');
            }, 500); // Show warning after success message
        }
        
        // Reset all state when new data is loaded
        resetAllState();
        
        // Clear MTX files and reset UI
        mtxFiles = { matrix: null, features: null, barcodes: null };
        
        // Reset MTX upload buttons and status
        const matrixBtn = document.getElementById('mtx-matrix-btn');
        const featuresBtn = document.getElementById('mtx-features-btn');
        const barcodesBtn = document.getElementById('mtx-barcodes-btn');
        const matrixStatus = document.getElementById('mtx-matrix-status');
        const featuresStatus = document.getElementById('mtx-features-status');
        const barcodesStatus = document.getElementById('mtx-barcodes-status');
        
        if (matrixBtn) {
            matrixBtn.classList.remove('uploaded');
            matrixBtn.textContent = 'Matrix:';
        }
        if (featuresBtn) {
            featuresBtn.classList.remove('uploaded');
            featuresBtn.textContent = 'Features:';
        }
        if (barcodesBtn) {
            barcodesBtn.classList.remove('uploaded');
            barcodesBtn.textContent = 'Barcodes:';
        }
        if (matrixStatus) matrixStatus.textContent = '';
        if (featuresStatus) featuresStatus.textContent = '';
        if (barcodesStatus) barcodesStatus.textContent = '';
        
        // Clear file inputs
        document.getElementById('mtx-matrix-file').value = '';
        document.getElementById('mtx-features-file').value = '';
        document.getElementById('mtx-barcodes-file').value = '';
        
        // Enable pipeline buttons
        document.getElementById('btn-preprocess').disabled = false;
        document.getElementById('pipeline-section').style.display = 'block';
        
        // Update metadata select
        if (data.obs_columns) {
            updateMetadataSelect(data.obs_columns);
        }
        
        // Update available cluster keys
        await updateAvailableClusterKeys();
        
        // Update button states
        updateButtonStates();
        
        // Update visualization method selector
        updateVisualizationMethodSelector();
        
    } catch (error) {
        console.error('[UPLOAD] Exception:', error);
        showStatus(`Error: ${error.message}`, 'error');
        // Clear files on error and reset UI
        mtxFiles = { matrix: null, features: null, barcodes: null };
        
        // Reset MTX upload buttons and status
        const matrixBtn = document.getElementById('mtx-matrix-btn');
        const featuresBtn = document.getElementById('mtx-features-btn');
        const barcodesBtn = document.getElementById('mtx-barcodes-btn');
        const matrixStatus = document.getElementById('mtx-matrix-status');
        const featuresStatus = document.getElementById('mtx-features-status');
        const barcodesStatus = document.getElementById('mtx-barcodes-status');
        
        if (matrixBtn) {
            matrixBtn.classList.remove('uploaded');
            matrixBtn.textContent = 'Matrix:';
        }
        if (featuresBtn) {
            featuresBtn.classList.remove('uploaded');
            featuresBtn.textContent = 'Features:';
        }
        if (barcodesBtn) {
            barcodesBtn.classList.remove('uploaded');
            barcodesBtn.textContent = 'Barcodes:';
        }
        if (matrixStatus) matrixStatus.textContent = '';
        if (featuresStatus) featuresStatus.textContent = '';
        if (barcodesStatus) barcodesStatus.textContent = '';
    } finally {
        // Always re-enable data loading buttons to allow reloading at any time
        enableDataLoadingButtons();
    }
}

async function runPreprocess() {
    showStatus('Running preprocessing... This may take a while.', 'info');
    const btn = document.getElementById('btn-preprocess');
    const originalText = btn.textContent;
    
    // Disable all buttons during preprocessing
    disableAllButtons();
    btn.innerHTML = '<span class="spinner"></span> Processing...';

    try {
        const n_pcs_value = document.getElementById('n-pcs').value?.trim();
        const n_pcs = n_pcs_value ? parseInt(n_pcs_value) : null;
        const n_neighbors = parseInt(document.getElementById('n-neighbors').value) || 15;
        const target_sum = parseFloat(document.getElementById('target-sum')?.value) || 1e4;
        const n_top_genes = document.getElementById('n-top-genes')?.value ? parseInt(document.getElementById('n-top-genes').value) : null;
        const use_rep = document.getElementById('preprocess-use-rep')?.value || 'X_pca';
        const save_raw = document.getElementById('save-raw')?.checked !== false;
        const min_genes = document.getElementById('min-genes')?.value ? parseInt(document.getElementById('min-genes').value) : null;
        const min_cells = document.getElementById('min-cells')?.value ? parseInt(document.getElementById('min-cells').value) : null;
        const min_counts = document.getElementById('min-counts')?.value ? parseInt(document.getElementById('min-counts').value) : null;
        const max_counts = document.getElementById('max-counts')?.value ? parseInt(document.getElementById('max-counts').value) : null;
        const max_genes = document.getElementById('max-genes')?.value ? parseInt(document.getElementById('max-genes').value) : null;
        const pct_mt_max = document.getElementById('pct-mt-max')?.value ? parseFloat(document.getElementById('pct-mt-max').value) : null;
        const hvg_flavor = document.getElementById('hvg-flavor')?.value || 'seurat_v3';
        const batch_key = document.getElementById('batch-key')?.value?.trim() || null;
        const regress_out_keys_str = document.getElementById('regress-out-keys')?.value?.trim() || null;
        const regress_out_keys = regress_out_keys_str ? regress_out_keys_str.split(',').map(k => k.trim()).filter(k => k) : null;
        const use_combat = document.getElementById('use-combat')?.checked || false;

        const requestBody = {
            session_id: sessionId,
            n_neighbors: n_neighbors,
            target_sum: target_sum,
            use_rep: use_rep,
            save_raw: save_raw,
            hvg_flavor: hvg_flavor
        };
        
        // Include n_pcs only if set (null means auto-select)
        if (n_pcs !== null) requestBody.n_pcs = n_pcs;
        
        // Only include optional parameters if they are set
        if (n_top_genes !== null) requestBody.n_top_genes = n_top_genes;
        if (min_genes !== null) requestBody.min_genes = min_genes;
        if (min_cells !== null) requestBody.min_cells = min_cells;
        if (min_counts !== null) requestBody.min_counts = min_counts;
        if (max_counts !== null) requestBody.max_counts = max_counts;
        if (max_genes !== null) requestBody.max_genes = max_genes;
        if (pct_mt_max !== null) requestBody.pct_mt_max = pct_mt_max;
        if (batch_key !== null) requestBody.batch_key = batch_key;
        if (regress_out_keys !== null && regress_out_keys.length > 0) requestBody.regress_out_keys = regress_out_keys;
        if (use_combat) requestBody.use_combat = use_combat;

        const response = await fetch(`${API_BASE}/preprocess`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(requestBody)
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Preprocessing failed');
        }

        showStatus('✓ Preprocessing complete!', 'success');
        
        // Display PCA variance ratio if available
        if (data.pca_variance_ratio && data.pca_variance_ratio.length > 0) {
            displayPCAVarianceRatio(data.pca_variance_ratio);
        }
        
        // Enable core analysis and cluster buttons (core analysis is independent)
        document.getElementById('btn-core-selection').disabled = false;
        document.getElementById('btn-cluster').disabled = false;
        
        // Update button states based on rules (visualization/marker genes still disabled until clustering/core analysis)
        updateButtonStates();
        
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
    
    // Disable all buttons during clustering
    disableAllButtons();
    btn.innerHTML = '<span class="spinner"></span> Clustering...';

    try {
        const method = document.getElementById('cluster-method').value;
        const resolution = parseFloat(document.getElementById('cluster-resolution').value) || 1.2;
        const use_rep = document.getElementById('use-rep').value || null;
        const key_added = null; // Use default (method name)
        
        // cplearn-specific parameters
        const stable_core_frac = parseFloat(document.getElementById('stable-core-frac')?.value) || 0.25;
        const stable_ng_num = parseInt(document.getElementById('stable-ng-num')?.value) || 8;
        const fine_grained = document.getElementById('fine-grained')?.checked || false;
        const propagate = document.getElementById('propagate')?.checked !== false;

        const requestBody = {
            session_id: sessionId,
            method: method,
            resolution: resolution
        };
        
        if (use_rep) requestBody.use_rep = use_rep;
        if (key_added) requestBody.key_added = key_added;
        
        // Add cplearn-specific parameters
        if (method === 'cplearn') {
            requestBody.stable_core_frac = stable_core_frac;
            requestBody.stable_ng_num = stable_ng_num;
            requestBody.fine_grained = fine_grained;
            requestBody.propagate = propagate;
        }

        const response = await fetch(`${API_BASE}/cluster`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(requestBody)
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Clustering failed');
        }

        showStatus(`✓ Clustering complete! Found ${data.n_clusters} clusters using ${data.cluster_key}.`, 'success');
        
        // Store cluster info
        // Preserve has_model from previous Core Analysis if it exists
        // This prevents coremap option from disappearing when running scanpy clustering
        const previousHasModel = clusterInfo && clusterInfo.has_model;
        clusterInfo = {
            cluster_key: data.cluster_key,
            clusters: data.clusters,
            cluster_counts: data.cluster_counts,
            has_model: data.has_model || previousHasModel || false
        };
        
        // If cplearn clustering, automatically run Core Analysis with same parameters (without visualization)
        if (method === 'cplearn' && data.has_model) {
            try {
                showStatus('Running Core Analysis with cplearn clustering parameters...', 'info');
                
                // Get parameters from Clustering form (same as used for clustering)
                const use_rep = document.getElementById('use-rep').value || null;
                const key_added = 'X_cplearn_coremap'; // Default coremap key
                const cluster_resolution = resolution; // Use same resolution as clustering
                const cluster_key = data.cluster_key; // Use the cluster key from clustering
                
                // Use same cplearn parameters from Clustering form
                const stable_core_frac = parseFloat(document.getElementById('stable-core-frac')?.value) || 0.25;
                const stable_ng_num = parseInt(document.getElementById('stable-ng-num')?.value) || 8;
                const fine_grained = document.getElementById('fine-grained')?.checked || false;
                const propagate = document.getElementById('propagate')?.checked !== false;
                
                // Call Core Analysis API with clustering parameters (without visualization)
                const coreResponse = await fetch(`${API_BASE}/core-analyze`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        session_id: sessionId,
                        use_rep: use_rep,
                        key_added: key_added,
                        cluster_resolution: cluster_resolution,
                        cluster_key: cluster_key,
                        stable_core_frac: stable_core_frac,
                        stable_ng_num: stable_ng_num,
                        fine_grained: fine_grained,
                        propagate: propagate,
                        visualize: false  // Don't generate visualization here, let Visualization module handle it
                    })
                });
                
                const coreData = await coreResponse.json();
                
                if (coreResponse.ok) {
                    // Update has_model flag
                    clusterInfo.has_model = true;
                    showStatus('✓ Core Analysis complete! Coremap embedding available for visualization.', 'success');
                } else {
                    console.warn('[CLUSTER] Core Analysis failed:', coreData.error || 'Unknown error');
                    showStatus('Clustering complete, but Core Analysis failed. You can run Core Analysis separately.', 'info');
                }
            } catch (coreError) {
                console.error('[CLUSTER] Error running Core Analysis after cplearn clustering:', coreError);
                // Don't fail the clustering operation if Core Analysis fails
                showStatus('Clustering complete, but Core Analysis failed. You can run Core Analysis separately.', 'info');
            }
        }
        
        // Update button states based on rules (now we have clustering, so visualization/marker genes can be enabled)
        updateButtonStates();
        
        // Update available cluster keys for DEG analysis (only existing ones)
        await updateAvailableClusterKeys();
        
        // Update visualization method selector (will show coremap if Core Analysis was run)
        updateVisualizationMethodSelector();
        
        // Update color-by select to include cluster key
        updateMetadataSelect([data.cluster_key]);
        document.getElementById('color-by').value = data.cluster_key;
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Clustering error:', error);
    } finally {
        // Restore button states
        btn.disabled = false;
        btn.textContent = originalText;
        
        // Re-enable buttons that should be available
        // Clustering doesn't prevent re-running preprocessing or core analysis
        // Users may want to compare different clustering parameters
        document.getElementById('btn-preprocess').disabled = false;
        document.getElementById('btn-core-selection').disabled = false;
        
        // Update button states based on rules (visualization/marker genes only if clustering/core analysis exists)
        updateButtonStates();
    }
}

async function runVisualization() {
    const params = getVisualizationParams();
    
    // For coremap, don't need cluster key
    if (params.method === 'coremap') {
        await runVisualizationForClusterKey(null);
        return;
    }
    
    // For UMAP, need cluster key
    const vizClusterKeySelect = document.getElementById('viz-cluster-key');
    const selectedKey = vizClusterKeySelect ? vizClusterKeySelect.value : null;
    
    if (selectedKey && selectedKey !== '') {
        // Use selected cluster key
        await runVisualizationForClusterKey(selectedKey);
    } else {
        // Use default cluster key
        const defaultKey = clusterInfo.cluster_key || 'cplearn';
        await runVisualizationForClusterKey(defaultKey);
    }
}

async function runVisualizationForClusterKey(cluster_key) {
    // Get current visualization parameters
    const params = getVisualizationParams();
    
    // Handle coremap method - generate visualization via API
    if (params.method === 'coremap') {
        // Check if coremap embedding exists (by checking if Core Analysis was run)
        if (!clusterInfo || !clusterInfo.has_model) {
            showStatus('Error: Core Analysis not available. Please run Core Analysis or cplearn clustering first.', 'error');
            return;
        }
        
    const btn = document.getElementById('btn-visualize');
    const originalText = btn.textContent;
        
        // Disable all buttons during visualization
        disableAllButtons();
        btn.innerHTML = '<span class="spinner"></span> Generating Coremap...';
        
        try {
            showStatus('Generating coremap visualization...', 'info');
            
            // Get Ground Truth file if provided
            const truthJsonFile = document.getElementById('viz-truth-json-file')?.files[0];
            let requestBody = {
                session_id: sessionId,
                cluster_key: cluster_key || 'cplearn',
                n_components: params.n_components,
                method: 'coremap',
                use_rep: 'X_cplearn_coremap'  // Default coremap key
            };
            
            let headers = { 'Content-Type': 'application/json' };
            let body;
            
            // If Ground Truth file is provided, use FormData
            if (truthJsonFile) {
                const formData = new FormData();
                formData.append('session_id', sessionId);
                formData.append('cluster_key', cluster_key || 'cplearn');
                formData.append('n_components', params.n_components.toString());
                formData.append('method', 'coremap');
                formData.append('use_rep', 'X_cplearn_coremap');
                formData.append('truth_json_file', truthJsonFile);
                
                headers = {}; // FormData sets Content-Type automatically
                body = formData;
            } else {
                body = JSON.stringify(requestBody);
            }
            
            // Call Visualization API to generate coremap visualization
        const response = await fetch(`${API_BASE}/visualize`, {
            method: 'POST',
                headers: headers,
                body: body
            });
            
            const data = await response.json();
            
            if (!response.ok) {
                throw new Error(data.error || 'Coremap visualization failed');
            }
            
            // Display the visualization in iframe
            const plotDiv = document.getElementById('plot');
            plotDiv.innerHTML = ''; // Clear existing content
            
            const iframe = document.createElement('iframe');
            iframe.style.width = '100%';
            iframe.style.height = '100%';
            iframe.style.border = 'none';
            iframe.srcdoc = data.html_content;
            plotDiv.appendChild(iframe);
            
            // Store visualization for future use
            if (data.file_url) {
                coreAnalysisVisualization = {
                    file_url: data.file_url,
                    htmlContent: data.html_content
                };
            }
            
            // Show visualization controls
            document.getElementById('visualization-section').style.display = 'block';
            document.getElementById('search-section').style.display = 'block';
            document.getElementById('stats-section').style.display = 'block';
            
            showStatus('✓ Coremap visualization generated!', 'success');
            
        } catch (error) {
            console.error('[VIZ] Error generating coremap visualization:', error);
            showStatus(`Error generating visualization: ${error.message}`, 'error');
        } finally {
            btn.disabled = false;
            btn.textContent = originalText;
            
            // Re-enable buttons that should be available
            document.getElementById('btn-marker-genes').disabled = false;
            document.getElementById('btn-cluster').disabled = false;
            document.getElementById('btn-core-selection').disabled = false;
            const btnPreprocess = document.getElementById('btn-preprocess');
            const pipelineSection = document.getElementById('pipeline-section');
            if (btnPreprocess && pipelineSection && pipelineSection.style.display !== 'none') {
                btnPreprocess.disabled = false;
            }
        }
        return;
    }
    
    // Handle standard visualization methods (umap, tsne, diffmap, draw_graph)
    const cacheKey = getVisualizationCacheKey(cluster_key, params.method, params.n_components, params.min_dist, params.spread);
    
    // Check if visualization is already cached with current parameters
    if (visualizationCache[cacheKey]) {
        vectors = visualizationCache[cacheKey].vectors;
        metadata = visualizationCache[cacheKey].metadata;
        updatePlot();
        showStatus(`✓ Loaded cached visualization for ${cluster_key}`, 'success');
        return;
    }
    
    const methodLabels = {
        'umap': 'UMAP',
        'tsne': 't-SNE',
        'diffmap': 'Diffusion Map',
        'draw_graph': 'Force-Directed Graph'
    };
    const methodLabel = methodLabels[params.method] || params.method.toUpperCase();
    
    showStatus(`Computing ${methodLabel} visualization for ${cluster_key}...`, 'info');
    const btn = document.getElementById('btn-visualize');
    const originalText = btn.textContent;
    
    // Disable all buttons during visualization
    disableAllButtons();
    btn.innerHTML = `<span class="spinner"></span> Computing ${methodLabel}...`;

    try {
        const requestBody = {
            session_id: sessionId,
            cluster_key: cluster_key,
            n_components: params.n_components,
            method: params.method
        };
        
        // Add method-specific parameters
        if (params.method === 'umap') {
            requestBody.min_dist = params.min_dist;
            requestBody.spread = params.spread;
        } else if (params.method === 'draw_graph') {
            requestBody.layout = 'fa'; // Default layout
        }
        
        const response = await fetch(`${API_BASE}/visualize`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(requestBody)
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Visualization failed');
        }

        // Store data in cache with parameters included
        visualizationCache[cacheKey] = {
            vectors: data.coordinates,
            metadata: data.metadata,
            params: params,
            method: params.method
        };
        
        // Update current vectors and metadata
        vectors = data.coordinates;
        metadata = data.metadata;
        
        // Update cluster key selector to show selected key
        const vizClusterKeySelect = document.getElementById('viz-cluster-key');
        if (vizClusterKeySelect) {
            vizClusterKeySelect.value = cluster_key;
        }
        
        // Automatically update Color By to match the cluster key
        const colorBySelect = document.getElementById('color-by');
        if (colorBySelect) {
            // Check if the cluster key exists in the color-by options
            const optionExists = Array.from(colorBySelect.options).some(opt => opt.value === cluster_key);
            if (!optionExists) {
                // Add the cluster key to color-by options if it doesn't exist
                const option = document.createElement('option');
                option.value = cluster_key;
                option.textContent = cluster_key;
                colorBySelect.appendChild(option);
            }
            // Set color-by to the cluster key
            colorBySelect.value = cluster_key;
        }
        
        showStatus(`✓ Visualization ready for ${cluster_key}!`, 'success');
        updatePlot();
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Visualization error:', error);
    } finally {
        // Restore button states
        btn.disabled = false;
        btn.textContent = originalText;
        
        // Re-enable buttons that should be available
        document.getElementById('btn-cluster').disabled = false;
        document.getElementById('btn-core-selection').disabled = false;
        // Preprocess button should remain enabled if data is loaded
        const btnPreprocess = document.getElementById('btn-preprocess');
        if (btnPreprocess && !btnPreprocess.disabled) {
            // If preprocess was enabled before, keep it enabled
        } else {
            // Otherwise, check if we have data (pipeline section is visible)
            const pipelineSection = document.getElementById('pipeline-section');
            if (pipelineSection && pipelineSection.style.display !== 'none') {
                btnPreprocess.disabled = false;
            }
        }
        
        // Update button states based on rules (visualization/marker genes only if clustering/core analysis exists)
        updateButtonStates();
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
            text: `Visualization (${vectors.length} points)`,
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

function displayPCAVarianceRatio(varianceRatio) {
    /**
     * Display PCA variance ratio plot in the visualization area
     */
    const plotDiv = document.getElementById('plot');
    if (!plotDiv) return;
    
    // Create PC numbers (1-indexed)
    const pcNumbers = varianceRatio.map((_, i) => i + 1);
    
    // Create trace for variance ratio
    const trace = {
        x: pcNumbers,
        y: varianceRatio,
        type: 'scatter',
        mode: 'lines+markers',
        name: 'Variance Ratio',
        line: {
            color: '#3b82f6',
            width: 2
        },
        marker: {
            color: '#3b82f6',
            size: 6
        }
    };
    
    // Calculate cumulative variance
    const cumulativeVariance = varianceRatio.reduce((acc, val, i) => {
        acc.push((i === 0 ? 0 : acc[i - 1]) + val);
        return acc;
    }, []);
    
    // Add cumulative variance trace
    const cumulativeTrace = {
        x: pcNumbers,
        y: cumulativeVariance,
        type: 'scatter',
        mode: 'lines',
        name: 'Cumulative Variance',
        yaxis: 'y2',
        line: {
            color: '#ef4444',
            width: 2,
            dash: 'dash'
        }
    };
    
    const layout = {
        title: {
            text: 'PCA Variance Ratio',
            font: { size: 18, color: '#1e293b' }
        },
        xaxis: {
            title: 'Principal Component',
            titlefont: { size: 14 },
            tickfont: { size: 12 }
        },
        yaxis: {
            title: 'Variance Ratio',
            titlefont: { size: 14 },
            tickfont: { size: 12 },
            type: 'log' // Log scale for better visualization
        },
        yaxis2: {
            title: 'Cumulative Variance',
            titlefont: { size: 14, color: '#ef4444' },
            tickfont: { size: 12, color: '#ef4444' },
            overlaying: 'y',
            side: 'right',
            range: [0, 1.1] // 0 to 110% for cumulative
        },
        hovermode: 'closest',
        showlegend: true,
        legend: {
            x: 0.7,
            y: 0.95
        },
        margin: { l: 60, r: 60, t: 60, b: 60 },
        plot_bgcolor: 'white',
        paper_bgcolor: 'white'
    };
    
    const config = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['lasso2d', 'select2d']
    };
    
    Plotly.newPlot('plot', [trace, cumulativeTrace], layout, config);
    
    // Show visualization section
    const visualizationSection = document.getElementById('visualization-section');
    if (visualizationSection) {
        visualizationSection.style.display = 'block';
    }
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
        showStatus('Running core analysis...', 'info');
    const btn = document.getElementById('btn-core-selection');
    const originalText = btn.textContent;
    
    // Disable all buttons during core analysis
    disableAllButtons();
    btn.innerHTML = '<span class="spinner"></span> Computing...';

    try {
        const use_rep = document.getElementById('core-use-rep')?.value || null;
        const cluster_resolution = parseFloat(document.getElementById('core-cluster-resolution')?.value) || 1.2;
        // Note: fast_view is now handled by Visualization module
        const key_added = 'X_cplearn_coremap'; // Use default
        const cluster_key = 'cplearn'; // Use default
        
        // cplearn clustering parameters (used when auto-running clustering)
        const stable_core_frac = parseFloat(document.getElementById('core-stable-core-frac')?.value) || 0.25;
        const stable_ng_num = parseInt(document.getElementById('core-stable-ng-num')?.value) || 8;
        const fine_grained = document.getElementById('core-fine-grained')?.checked || false;
        const propagate = document.getElementById('core-propagate')?.checked !== false;
        
        // Note: Ground Truth is now handled in Visualization module when using coremap
        
        // Prepare JSON request
        const requestBody = {
                session_id: sessionId,
                use_rep: use_rep,
                key_added: key_added,
                cluster_key: cluster_key,
                cluster_resolution: cluster_resolution,
            // Note: fast_view is now handled by Visualization module
                stable_core_frac: stable_core_frac,
                stable_ng_num: stable_ng_num,
                fine_grained: fine_grained,
                propagate: propagate
            };
            
        const headers = { 'Content-Type': 'application/json' };
        const body = JSON.stringify(requestBody);

        // Add timeout and retry logic for long-running operations
        const controller = new AbortController();
        const timeoutId = setTimeout(() => controller.abort(), 600000); // 10 minutes timeout
        
        let response;
        let data;
        let retryCount = 0;
        const maxRetries = 2;
        
        while (retryCount <= maxRetries) {
            try {
                response = await fetch(`${API_BASE}/core-analyze`, {
                    method: 'POST',
                    headers: headers,
                    body: body,
                    signal: controller.signal
                });
                
                clearTimeout(timeoutId);
                
                // Check if response is actually JSON
                const contentType = response.headers.get('content-type');
                const isJson = contentType && contentType.includes('application/json');
                
                // Handle non-JSON responses (502 Bad Gateway, 504 Gateway Timeout, etc.)
                if (!isJson) {
                    // Clone response to read text without consuming the body
                    const text = await response.clone().text();
                    console.error('[CORE] Non-JSON response:', text.substring(0, 200));
                    
                    if (response.status === 502 || response.status === 504) {
                        // Bad Gateway or Gateway Timeout - retry
                        if (retryCount < maxRetries) {
                            retryCount++;
                            const waitTime = Math.min(1000 * Math.pow(2, retryCount), 10000); // Exponential backoff, max 10s
                            showStatus(`Server timeout (${response.status}). Retrying in ${waitTime/1000}s... (attempt ${retryCount + 1}/${maxRetries + 1})`, 'info');
                            await new Promise(resolve => setTimeout(resolve, waitTime));
                            continue;
                        } else {
                            throw new Error(`Server error (${response.status}): The request took too long or the server is overloaded. Please try again later or use a smaller dataset.`);
                        }
                    } else {
                        throw new Error(`Server error (${response.status}): ${text.substring(0, 100)}`);
                    }
                }
                
                // Parse JSON response
                try {
                    data = await response.json();
                } catch (jsonError) {
                    // If JSON parsing fails, clone response to get text
                    const text = await response.clone().text();
                    console.error('[CORE] JSON parse error. Response text:', text.substring(0, 500));
                    throw new Error(`Invalid response from server: ${text.substring(0, 100)}`);
                }
                
                if (!response.ok) {
                    throw new Error(data.error || `Core analysis failed with status ${response.status}`);
                }
                
                // Success - break out of retry loop
                break;
                
            } catch (error) {
                clearTimeout(timeoutId);
                
                if (error.name === 'AbortError') {
                    throw new Error('Request timeout: The operation took longer than 10 minutes. Please try with a smaller dataset or contact support.');
                }
                
                // Network error or other fetch error
                if (retryCount < maxRetries && (error.message.includes('fetch') || error.message.includes('network'))) {
                    retryCount++;
                    const waitTime = Math.min(1000 * Math.pow(2, retryCount), 10000);
                    showStatus(`Network error. Retrying in ${waitTime/1000}s... (attempt ${retryCount + 1}/${maxRetries + 1})`, 'info');
                    await new Promise(resolve => setTimeout(resolve, waitTime));
                    continue;
                }
                
                // Re-throw if no more retries or non-retryable error
                throw error;
            }
        }

        let statusMessage = `✓ Core analysis complete! Embedding stored in ${data.key_added}`;
        if (data.assigned_points !== undefined) {
            statusMessage += ` (${data.assigned_points}/${data.total_points} points assigned`;
            if (data.core_cells !== undefined) {
                statusMessage += `, ${data.core_cells} core cells`;
            }
            statusMessage += ')';
        }
        
        showStatus(statusMessage, 'success');
        console.log('Core analysis result:', data);
        
        // Note: Visualization is now handled by Visualization module
        // No need to store visualization data here
        
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
            
            // Update color-by select to include cluster key and set it as selected
            if (data.cluster_key) {
                updateMetadataSelect([data.cluster_key]);
                document.getElementById('color-by').value = data.cluster_key;
            }
        }
        
        // Update button states based on rules (now we have core analysis, so visualization/marker genes can be enabled)
        updateButtonStates();
        
        // Update visualization method selector to show coremap option
        updateVisualizationMethodSelector();
        
    } catch (error) {
        showStatus(`Error: ${error.message}`, 'error');
        console.error('Core analysis error:', error);
    } finally {
        // Restore button states
        btn.disabled = false;
        btn.textContent = originalText;
        
        // Re-enable buttons that should be available
        // Core Analysis can run independently
        document.getElementById('btn-cluster').disabled = false;
        // Preprocess button should remain enabled if data is loaded
        const btnPreprocess = document.getElementById('btn-preprocess');
        if (btnPreprocess && !btnPreprocess.disabled) {
            // If preprocess was enabled before, keep it enabled
        } else {
            // Otherwise, check if we have data (pipeline section is visible)
            const pipelineSection = document.getElementById('pipeline-section');
            if (pipelineSection && pipelineSection.style.display !== 'none') {
                btnPreprocess.disabled = false;
            }
        }
        
        // Update button states based on rules (visualization/marker genes only if clustering/core analysis exists)
        updateButtonStates();
    }
}

async function runMarkerGenes() {
    showStatus('Finding marker genes...', 'info');
    const btn = document.getElementById('btn-marker-genes');
    const originalText = btn.textContent;
    
    // Disable all buttons during marker genes analysis
    disableAllButtons();
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
        // Restore button states
        btn.disabled = false;
        btn.textContent = originalText;
        
        // Re-enable buttons that should be available
        document.getElementById('btn-cluster').disabled = false;
        document.getElementById('btn-core-selection').disabled = false;
        // Preprocess button should remain enabled if data is loaded
        const btnPreprocess = document.getElementById('btn-preprocess');
        if (btnPreprocess && !btnPreprocess.disabled) {
            // If preprocess was enabled before, keep it enabled
        } else {
            // Otherwise, check if we have data (pipeline section is visible)
            const pipelineSection = document.getElementById('pipeline-section');
            if (pipelineSection && pipelineSection.style.display !== 'none') {
                btnPreprocess.disabled = false;
            }
        }
        
        // Update button states based on rules (visualization/marker genes only if clustering/core analysis exists)
        updateButtonStates();
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
    // Create toast container if it doesn't exist
    let toastContainer = document.getElementById('toast-container');
    if (!toastContainer) {
        toastContainer = document.createElement('div');
        toastContainer.id = 'toast-container';
        toastContainer.className = 'toast-container';
        document.body.appendChild(toastContainer);
    }

    // Create toast element
    const toast = document.createElement('div');
    toast.className = `toast ${type}`;

    // Add icon based on type
    const icon = document.createElement('span');
    icon.className = 'toast-icon';
    if (type === 'success') {
        icon.textContent = '✓';
    } else if (type === 'error') {
        icon.textContent = '✗';
    } else if (type === 'warning') {
        icon.textContent = '⚠';
    } else {
        icon.textContent = 'ℹ';
    }

    // Add message
    const messageSpan = document.createElement('span');
    messageSpan.className = 'toast-message';
    messageSpan.textContent = message;

    // Add close button
    const closeBtn = document.createElement('button');
    closeBtn.className = 'toast-close';
    closeBtn.innerHTML = '×';
    closeBtn.onclick = () => removeToast(toast);

    // Assemble toast
    toast.appendChild(icon);
    toast.appendChild(messageSpan);
    toast.appendChild(closeBtn);

    // Add to container
    toastContainer.appendChild(toast);

    // Auto-remove after delay (errors don't auto-remove, only success/info)
    if (type !== 'error') {
        const delay = type === 'info' ? 5000 : 4000;
        setTimeout(() => {
            removeToast(toast);
        }, delay);
    }
}

function removeToast(toast) {
    if (toast && toast.parentNode) {
        toast.classList.add('slide-out');
        setTimeout(() => {
            if (toast.parentNode) {
                toast.parentNode.removeChild(toast);
            }
        }, 300);
    }
}
>>>>>>> 8b02f6559a88ee0942e4fa185dfd86d6e1447c40
