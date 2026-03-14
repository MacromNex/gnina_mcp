"""Shared utilities for Gnina MCP Server."""

from pathlib import Path
from typing import Dict, Any, Optional
import os
import sys
import shutil

def setup_paths():
    """Setup paths for script imports and environment."""
    script_dir = Path(__file__).parent.resolve()
    mcp_root = script_dir.parent
    scripts_dir = mcp_root / "scripts"

    # Add to Python path
    sys.path.insert(0, str(script_dir))
    sys.path.insert(0, str(scripts_dir))

    return {
        "script_dir": script_dir,
        "mcp_root": mcp_root,
        "scripts_dir": scripts_dir
    }

def resolve_gnina_executable():
    """Find gnina executable, checking PATH, env/bin, and mock fallback."""
    if shutil.which('gnina'):
        return 'gnina'
    script_dir = Path(__file__).parent.resolve()
    mcp_root = script_dir.parent
    # Check the conda env bundled with the project
    env_path = mcp_root / 'env' / 'bin' / 'gnina'
    if env_path.exists():
        return str(env_path)
    # Fall back to mock_gnina.py for testing
    mock_path = mcp_root / 'mock_gnina.py'
    if mock_path.exists():
        return str(mock_path)
    return 'gnina'

def gnina_subprocess_env():
    """Get environment dict with LD_LIBRARY_PATH set for gnina."""
    env = os.environ.copy()
    mcp_root = Path(__file__).resolve().parent.parent
    env_lib = mcp_root / 'env' / 'lib'
    if env_lib.exists():
        ld_path = env.get('LD_LIBRARY_PATH', '')
        if str(env_lib) not in ld_path:
            env['LD_LIBRARY_PATH'] = f"{env_lib}:{ld_path}" if ld_path else str(env_lib)
    return env

def resolve_python_executable():
    """Find python executable in the project's conda env."""
    script_dir = Path(__file__).parent.resolve()
    mcp_root = script_dir.parent
    env_python = mcp_root / 'env' / 'bin' / 'python'
    if env_python.exists():
        return str(env_python)
    return sys.executable

def validate_file_path(file_path: str, must_exist: bool = True) -> Optional[str]:
    """Validate and return absolute file path."""
    if not file_path:
        return None

    path = Path(file_path).resolve()

    if must_exist and not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    return str(path)

def ensure_output_dir(output_path: str) -> str:
    """Ensure output directory exists and return absolute path."""
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return str(output_path)

def format_error_response(error_msg: str, error_type: str = "error") -> Dict[str, Any]:
    """Format standardized error response."""
    return {
        "status": "error",
        "error_type": error_type,
        "error": str(error_msg)
    }

def format_success_response(result: Any, **kwargs) -> Dict[str, Any]:
    """Format standardized success response."""
    response = {
        "status": "success",
        "result": result
    }
    response.update(kwargs)
    return response