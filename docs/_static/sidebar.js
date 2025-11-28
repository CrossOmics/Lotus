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
        
        // Find the API Reference link
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

        // Check if API Reference already has a submenu
        let submenu = apiRefItem.querySelector('ul');
        
        // If no submenu, try to find API module links that should be grouped
        if (!submenu) {
            // Look for API module links in the sidebar
            const apiModules = [];
            const allListItems = sidebarTree.querySelectorAll('li');
            let foundApiRef = false;
            
            for (const item of allListItems) {
                const link = item.querySelector('a');
                if (!link) continue;
                
                const href = link.getAttribute('href') || '';
                
                if (item === apiRefItem) {
                    foundApiRef = true;
                    continue;
                }
                
                if (foundApiRef) {
                    // Check if this is an API module
                    if (href.includes('api/') && 
                        !href.includes('api/index') &&
                        (href.includes('workflows') ||
                         href.includes('preprocess') ||
                         href.includes('core_selection') ||
                         href.includes('clustering') ||
                         href.includes('visualization') ||
                         href.includes('deg') ||
                         href.includes('tools') ||
                         href.includes('preprocessing_compat'))) {
                        apiModules.push(item);
                    } else if (!href.includes('api/')) {
                        // Stop when we hit a non-API link
                        break;
                    }
                }
            }
            
            // Create submenu if we found API modules
            if (apiModules.length > 0) {
                submenu = document.createElement('ul');
                submenu.className = 'api-ref-submenu';
                
                apiModules.forEach(li => {
                    const clonedLi = li.cloneNode(true);
                    submenu.appendChild(clonedLi);
                });
                
                apiRefItem.appendChild(submenu);
            }
        }

        // Only proceed if we have a submenu with items
        if (!submenu || submenu.children.length === 0) {
            console.log('No submenu found or submenu is empty');
            return;
        }

        // Check if toggle button already exists
        if (apiRefItem.querySelector('.api-ref-toggle')) {
            return;
        }

        // Create toggle button
        const toggleBtn = document.createElement('button');
        toggleBtn.className = 'api-ref-toggle';
        toggleBtn.setAttribute('aria-expanded', 'false');
        toggleBtn.setAttribute('aria-label', 'Toggle API Reference submenu');
        toggleBtn.setAttribute('type', 'button');
        toggleBtn.innerHTML = '<span class="toggle-icon">▼</span>';
        
        // Insert toggle button before the link
        if (apiRefLink.parentNode === apiRefItem) {
            apiRefItem.insertBefore(toggleBtn, apiRefLink);
        } else {
            // If link is wrapped in another element, insert before that wrapper
            const linkWrapper = apiRefLink.parentNode;
            linkWrapper.parentNode.insertBefore(toggleBtn, linkWrapper);
        }

        // Initially hide the submenu
        submenu.style.display = 'none';
        submenu.classList.add('api-ref-submenu');

        // Toggle function
        function toggleSubmenu() {
            const isExpanded = toggleBtn.getAttribute('aria-expanded') === 'true';
            const newState = !isExpanded;
            
            toggleBtn.setAttribute('aria-expanded', newState);
            submenu.style.display = newState ? 'block' : 'none';
            toggleBtn.querySelector('.toggle-icon').textContent = newState ? '▲' : '▼';
            apiRefItem.classList.toggle('expanded', newState);
        }

        // Click handler for toggle button
        toggleBtn.addEventListener('click', function(e) {
            e.preventDefault();
            e.stopPropagation();
            toggleSubmenu();
        });

        // Also allow clicking on the API Reference link to toggle (when on API index page)
        const currentPath = window.location.pathname;
        if (currentPath.includes('/api/index') || currentPath.endsWith('/api/')) {
            apiRefLink.addEventListener('click', function(e) {
                e.preventDefault();
                toggleSubmenu();
            });
        }

        // Expand by default if we're on an API page
        if (currentPath.includes('/api/') && !currentPath.includes('/api/index')) {
            toggleSubmenu();
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
