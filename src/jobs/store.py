"""Job state persistence utilities."""

from pathlib import Path
from typing import Dict, Any, Optional
import json
from datetime import datetime

class JobStore:
    """Simple file-based job state persistence."""

    def __init__(self, store_dir: Path = None):
        self.store_dir = store_dir or Path(__file__).parent.parent.parent / "jobs"
        self.store_dir.mkdir(parents=True, exist_ok=True)

    def save_job_data(self, job_id: str, data: Dict[str, Any]) -> bool:
        """Save arbitrary job data to disk."""
        try:
            job_dir = self.store_dir / job_id
            job_dir.mkdir(parents=True, exist_ok=True)

            data_file = job_dir / "job_data.json"
            data["last_updated"] = datetime.now().isoformat()

            with open(data_file, 'w') as f:
                json.dump(data, f, indent=2)

            return True
        except Exception:
            return False

    def load_job_data(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Load arbitrary job data from disk."""
        try:
            data_file = self.store_dir / job_id / "job_data.json"
            if data_file.exists():
                with open(data_file) as f:
                    return json.load(f)
        except Exception:
            pass
        return None

    def delete_job_data(self, job_id: str) -> bool:
        """Delete job data from disk."""
        try:
            job_dir = self.store_dir / job_id
            if job_dir.exists():
                import shutil
                shutil.rmtree(job_dir)
            return True
        except Exception:
            return False

# Global store instance
job_store = JobStore()