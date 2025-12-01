// Sidebar dropdown menu functionality for API Reference
(function() {
    'use strict';
    
    function initSidebarDropdown() {
        // Try multiple selectors for Furo theme sidebar
        const sidebarTree = document.querySelector('.sidebar-tree') ||
                           document.querySelector('.sidebar-drawer .sidebar-tree') ||
                           document.querySelector('nav .sidebar-tree');
        
        if (!sidebarTree) {
            console.log('Sidebar tree not found, retrying...');
            setTimeout(initSidebarDropdown, 200);
            return;
        }

        // Find all links in sidebar
        const allLinks = sidebarTree.querySelectorAll('a');
        let apiRefLink = null;
        let apiRefItem = null;
        
        // Find the API Reference link and its parent list item
        for (const link of allLinks) {
            const href = link.getAttribute('href') || '';
            const text = link.textContent.trim();
            
            // Match API Reference
            if (href.includes('api/index') || 
                text === 'API Reference' ||
                (href.includes('api/') && text.toLowerCase().includes('api'))) {
                apiRefLink = link;
                apiRefItem = link.closest('li');
                break;
            }
        }

        if (!apiRefItem || !apiRefLink) {
            console.log('API Reference item not found');
            return;
        }

        // Check if API Reference already has a submenu (from Sphinx toctree)
        let submenu = apiRefItem.querySelector('ul');
        
        // If no submenu from toctree, the toctree might not be rendering correctly
        // In that case, we should not interfere with Furo's native toctree functionality
        if (!submenu) {
            console.log('No submenu found - toctree may not be rendering correctly');
            return;
        }

        // Check if Furo's native toggle button already exists
        // Furo uses a button with class 'toctree-toggle' or similar
        const existingToggle = apiRefItem.querySelector('button[aria-expanded]') ||
                              apiRefItem.querySelector('.toctree-toggle') ||
                              apiRefItem.querySelector('button[class*="toggle"]');
        
        if (existingToggle) {
            // Furo's native toggle exists, ensure it works correctly
            console.log('Furo native toggle found, ensuring it works');
            
            // Make sure the toggle button is clickable
            existingToggle.addEventListener('click', function(e) {
                e.stopPropagation();
                // Let Furo handle the toggle
            }, true);
            
            return;
        }

        // If no native toggle exists, create one that works with Furo's structure
        // Find the label element (Furo wraps links in labels)
        const label = apiRefLink.closest('label') || apiRefLink.parentElement;
        
        // Create toggle button compatible with Furo
        const toggleBtn = document.createElement('button');
        toggleBtn.className = 'toctree-toggle';
        toggleBtn.setAttribute('aria-expanded', 'false');
        toggleBtn.setAttribute('aria-label', 'Toggle API Reference submenu');
        toggleBtn.setAttribute('type', 'button');
        toggleBtn.innerHTML = '▼';
        
        // Insert toggle button before the label/link
        if (label && label.parentNode === apiRefItem) {
            apiRefItem.insertBefore(toggleBtn, label);
        } else {
            apiRefItem.insertBefore(toggleBtn, apiRefLink);
        }

        // Initially hide the submenu if not expanded
        const isExpanded = apiRefItem.classList.contains('current') || 
                          apiRefItem.classList.contains('expanded');
        if (!isExpanded) {
            submenu.style.display = 'none';
        }
        toggleBtn.setAttribute('aria-expanded', isExpanded ? 'true' : 'false');
        toggleBtn.textContent = isExpanded ? '▲' : '▼';

        // Toggle function
        function toggleSubmenu() {
            const isExpanded = toggleBtn.getAttribute('aria-expanded') === 'true';
            const newState = !isExpanded;
            
            toggleBtn.setAttribute('aria-expanded', newState);
            submenu.style.display = newState ? 'block' : 'none';
            toggleBtn.textContent = newState ? '▲' : '▼';
            apiRefItem.classList.toggle('expanded', newState);
        }

        // Click handler for toggle button
        toggleBtn.addEventListener('click', function(e) {
            e.preventDefault();
            e.stopPropagation();
            toggleSubmenu();
        });

        // Expand by default if we're on an API page
        const currentPath = window.location.pathname;
        if (currentPath.includes('/api/') && !currentPath.includes('/api/index')) {
            if (!isExpanded) {
                toggleSubmenu();
            }
        }
        
        console.log('Sidebar dropdown initialized successfully');
    }

    // Initialize when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', initSidebarDropdown);
    } else {
        // DOM is already ready
        initSidebarDropdown();
    }
    
    // Also try after a short delay in case Sphinx is still rendering
    setTimeout(initSidebarDropdown, 500);
})();

// Style Return type sections with grey background
(function() {
    'use strict';
    
    function styleReturnType() {
        // Find all field-list dt elements
        const fieldListDts = document.querySelectorAll('.field-list dt');
        
        fieldListDts.forEach(dt => {
            // Check if this dt contains "Return type" text
            const text = dt.textContent.trim();
            if (text === 'Return type' || text.startsWith('Return type')) {
                dt.style.backgroundColor = 'rgba(128, 128, 128, 0.15)';
                dt.style.borderLeft = '4px solid rgba(128, 128, 128, 0.4)';
            }
        });
    }
    
    // Run when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', styleReturnType);
    } else {
        styleReturnType();
    }
    
    // Also try after a delay in case content is loaded dynamically
    setTimeout(styleReturnType, 500);
    setTimeout(styleReturnType, 1000);
})();
