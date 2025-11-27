// Sidebar dropdown menu functionality for API Reference
document.addEventListener('DOMContentLoaded', function() {
    // Wait a bit for Sphinx to render the sidebar
    setTimeout(function() {
        // Find the sidebar tree - Furo theme uses .sidebar-tree
        const sidebarTree = document.querySelector('.sidebar-tree');
        if (!sidebarTree) return;

        // Find all list items
        const allItems = sidebarTree.querySelectorAll('li');
        let apiRefItem = null;
        
        // Find the API Reference item
        for (const item of allItems) {
            const link = item.querySelector('a');
            if (link) {
                const href = link.getAttribute('href') || '';
                const text = link.textContent.trim();
                // Check if this is the API Reference link
                if (href.includes('api/index') || 
                    (text === 'API Reference' && href.includes('api/'))) {
                    apiRefItem = item;
                    break;
                }
            }
        }

        if (!apiRefItem) return;

        // Check if API Reference already has a submenu (ul)
        let submenu = apiRefItem.querySelector('ul');
        
        // If no submenu, check if there are API module links as siblings
        if (!submenu) {
            // Look for the next sibling items that are API modules
            const apiModules = [];
            let current = apiRefItem.nextElementSibling;
            
            while (current && current.tagName === 'LI') {
                const link = current.querySelector('a');
                if (link) {
                    const href = link.getAttribute('href') || '';
                    if (href.includes('api/') && !href.includes('api/index')) {
                        apiModules.push(current);
                        current = current.nextElementSibling;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }
            
            // If we found API modules, create a submenu and move them
            if (apiModules.length > 0) {
                submenu = document.createElement('ul');
                submenu.className = 'api-ref-submenu';
                apiModules.forEach(li => {
                    submenu.appendChild(li.cloneNode(true));
                    li.remove();
                });
                apiRefItem.appendChild(submenu);
            }
        }

        if (!submenu || submenu.children.length === 0) return;

        // Create toggle button
        const toggleBtn = document.createElement('button');
        toggleBtn.className = 'api-ref-toggle';
        toggleBtn.setAttribute('aria-expanded', 'false');
        toggleBtn.setAttribute('aria-label', 'Toggle API Reference submenu');
        toggleBtn.innerHTML = '<span class="toggle-icon">▼</span>';
        
        // Insert toggle button
        const link = apiRefItem.querySelector('a');
        if (link && link.parentNode === apiRefItem) {
            apiRefItem.insertBefore(toggleBtn, link);
        } else if (link) {
            link.parentNode.insertBefore(toggleBtn, link);
        } else {
            apiRefItem.insertBefore(toggleBtn, apiRefItem.firstChild);
        }

        // Initially hide the submenu
        submenu.style.display = 'none';
        submenu.classList.add('api-ref-submenu');

        // Toggle function
        function toggleSubmenu() {
            const isExpanded = toggleBtn.getAttribute('aria-expanded') === 'true';
            toggleBtn.setAttribute('aria-expanded', !isExpanded);
            submenu.style.display = isExpanded ? 'none' : 'block';
            toggleBtn.querySelector('.toggle-icon').textContent = isExpanded ? '▼' : '▲';
            apiRefItem.classList.toggle('expanded', !isExpanded);
        }

        // Click handler for toggle button
        toggleBtn.addEventListener('click', function(e) {
            e.preventDefault();
            e.stopPropagation();
            toggleSubmenu();
        });

        // Also allow clicking on the API Reference link to toggle
        if (link) {
            const originalClick = link.onclick;
            link.addEventListener('click', function(e) {
                // Only prevent default if we're on the API index page
                if (window.location.pathname.includes('api/index') || 
                    link.getAttribute('href') === 'api/index.html') {
                    e.preventDefault();
                    toggleSubmenu();
                }
            });
        }

        // Check if we're currently on an API page and expand if so
        const currentPath = window.location.pathname;
        if (currentPath.includes('/api/')) {
            toggleSubmenu(); // Expand by default if on API page
        }
    }, 100);
});
