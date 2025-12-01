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
