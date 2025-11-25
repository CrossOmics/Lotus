"""
Session management API endpoints
Handles session creation, heartbeat, and cleanup
"""

import gc
from flask import Blueprint, request, jsonify
from pathlib import Path
import shutil
import time
import threading
from .utils import get_session_dir
from .config import UPLOAD_FOLDER

bp = Blueprint('session', __name__, url_prefix='/api')

# Track session last access times
# Format: {session_id: timestamp}
_session_access_times = {}
_session_lock = threading.Lock()

# Session timeout: 30 minutes of inactivity
SESSION_TIMEOUT = 30 * 60  # 30 minutes in seconds


@bp.route('/session/create', methods=['POST'])
def create_session():
    """Create a new session"""
    try:
        data = request.json
        session_id = data.get('session_id')
        
        if not session_id:
            return jsonify({'error': 'session_id is required'}), 400
        
        # Create session directory
        session_dir = get_session_dir(session_id)
        session_dir.mkdir(exist_ok=True, parents=True)
        
        # Update access time
        with _session_lock:
            _session_access_times[session_id] = time.time()
        
        print(f"[SESSION] Created session: {session_id}")
        return jsonify({
            'success': True,
            'session_id': session_id,
            'message': 'Session created'
        })
    
    except Exception as e:
        print(f"[SESSION] Error creating session: {e}")
        return jsonify({'error': str(e)}), 500


@bp.route('/session/heartbeat', methods=['POST'])
def heartbeat():
    """Update session last access time"""
    try:
        data = request.json
        session_id = data.get('session_id')
        
        if not session_id:
            return jsonify({'error': 'session_id is required'}), 400
        
        # Update access time
        with _session_lock:
            _session_access_times[session_id] = time.time()
        
        return jsonify({
            'success': True,
            'session_id': session_id,
            'message': 'Heartbeat received'
        })
    
    except Exception as e:
        print(f"[SESSION] Error updating heartbeat: {e}")
        return jsonify({'error': str(e)}), 500


@bp.route('/session/cleanup', methods=['POST'])
def cleanup_session():
    """Clean up a session (delete all files)"""
    try:
        data = request.json
        session_id = data.get('session_id')
        
        if not session_id:
            return jsonify({'error': 'session_id is required'}), 400
        
        # Remove from access times
        with _session_lock:
            if session_id in _session_access_times:
                del _session_access_times[session_id]
        
        # Delete session directory
        session_dir = get_session_dir(session_id)
        if session_dir.exists():
            try:
                shutil.rmtree(session_dir)
                print(f"[SESSION] Cleaned up session: {session_id}")
            except Exception as e:
                print(f"[SESSION] Warning: Could not fully delete session directory {session_id}: {e}")
        
        return jsonify({
            'success': True,
            'session_id': session_id,
            'message': 'Session cleaned up'
        })
    
    except Exception as e:
        print(f"[SESSION] Error cleaning up session: {e}")
        return jsonify({'error': str(e)}), 500


def cleanup_expired_sessions():
    """Clean up sessions that have been inactive for more than SESSION_TIMEOUT"""
    current_time = time.time()
    expired_sessions = []
    
    with _session_lock:
        for session_id, last_access in list(_session_access_times.items()):
            if current_time - last_access > SESSION_TIMEOUT:
                expired_sessions.append(session_id)
        
        # Remove expired sessions from tracking
        for session_id in expired_sessions:
            del _session_access_times[session_id]
    
    # Delete expired session directories
    for session_id in expired_sessions:
        session_dir = get_session_dir(session_id)
        if session_dir.exists():
            try:
                shutil.rmtree(session_dir)
                print(f"[SESSION] Auto-cleaned expired session: {session_id}")
            except Exception as e:
                print(f"[SESSION] Warning: Could not delete expired session {session_id}: {e}")
    
    # Force garbage collection after cleanup to free memory
    if expired_sessions:
        gc.collect()


# Background thread to periodically clean up expired sessions
def start_cleanup_thread():
    """Start background thread for periodic session cleanup"""
    def cleanup_loop():
        while True:
            try:
                cleanup_expired_sessions()
            except Exception as e:
                print(f"[SESSION] Error in cleanup thread: {e}")
            time.sleep(60)  # Check every minute
    
    cleanup_thread = threading.Thread(target=cleanup_loop, daemon=True)
    cleanup_thread.start()
    print("[SESSION] Background cleanup thread started")


# Start cleanup thread when module is imported
start_cleanup_thread()

