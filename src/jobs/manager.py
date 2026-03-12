"""Job management for long-running molecular docking tasks."""

import uuid
import json
import subprocess
import threading
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any
from enum import Enum
import fcntl
from loguru import logger

class JobStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"

class JobManager:
    """Manages asynchronous job execution for molecular docking computations."""

    def __init__(self, jobs_dir: Path = None):
        self.jobs_dir = jobs_dir or Path(__file__).parent.parent.parent / "jobs"
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self._running_jobs: Dict[str, subprocess.Popen] = {}
        self._lock = threading.Lock()

    def submit_job(
        self,
        script_path: str,
        args: Dict[str, Any],
        job_name: str = None
    ) -> Dict[str, Any]:
        """Submit a new job for background execution.

        Args:
            script_path: Path to the script to run
            args: Arguments to pass to the script
            job_name: Optional name for the job

        Returns:
            Dict with job_id and status
        """
        job_id = str(uuid.uuid4())[:8]
        job_dir = self.jobs_dir / job_id
        job_dir.mkdir(parents=True, exist_ok=True)

        # Save job metadata
        metadata = {
            "job_id": job_id,
            "job_name": job_name or f"job_{job_id}",
            "script": script_path,
            "args": args,
            "status": JobStatus.PENDING.value,
            "submitted_at": datetime.now().isoformat(),
            "started_at": None,
            "completed_at": None,
            "error": None
        }

        self._save_metadata(job_id, metadata)

        # Start job in background
        self._start_job(job_id, script_path, args, job_dir)

        return {
            "status": "submitted",
            "job_id": job_id,
            "message": f"Job submitted. Use get_job_status('{job_id}') to check progress."
        }

    def _start_job(self, job_id: str, script_path: str, args: Dict, job_dir: Path):
        """Start job execution in background thread."""
        def run_job():
            with self._lock:
                metadata = self._load_metadata(job_id)
                metadata["status"] = JobStatus.RUNNING.value
                metadata["started_at"] = datetime.now().isoformat()
                self._save_metadata(job_id, metadata)

            try:
                # Build command
                cmd = ["python", script_path]
                for key, value in args.items():
                    if value is not None and key != "output_dir":
                        if isinstance(value, bool):
                            if value:  # Only add flag if True
                                cmd.append(f"--{key}")
                        else:
                            cmd.extend([f"--{key}", str(value)])

                # Set output paths
                output_file = job_dir / "output.json"
                cmd.extend(["--output", str(output_file)])

                # Run script in the correct environment
                log_file = job_dir / "job.log"
                with open(log_file, 'w') as log:
                    # Use mamba run to ensure correct environment
                    full_cmd = ["mamba", "run", "-p", str(Path.cwd() / "env")] + cmd
                    process = subprocess.Popen(
                        full_cmd,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        cwd=str(Path(script_path).parent.parent)
                    )

                    with self._lock:
                        self._running_jobs[job_id] = process

                    process.wait()

                # Update status
                with self._lock:
                    metadata = self._load_metadata(job_id)
                    if process.returncode == 0:
                        metadata["status"] = JobStatus.COMPLETED.value
                    else:
                        metadata["status"] = JobStatus.FAILED.value
                        metadata["error"] = f"Process exited with code {process.returncode}"

            except Exception as e:
                with self._lock:
                    metadata = self._load_metadata(job_id)
                    metadata["status"] = JobStatus.FAILED.value
                    metadata["error"] = str(e)
                    logger.error(f"Job {job_id} failed: {e}")

            finally:
                with self._lock:
                    metadata["completed_at"] = datetime.now().isoformat()
                    self._save_metadata(job_id, metadata)
                    self._running_jobs.pop(job_id, None)

        thread = threading.Thread(target=run_job, daemon=True)
        thread.start()

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """Get status of a submitted job."""
        metadata = self._load_metadata(job_id)
        if not metadata:
            return {"status": "error", "error": f"Job {job_id} not found"}

        result = {
            "job_id": job_id,
            "job_name": metadata.get("job_name"),
            "status": metadata["status"],
            "submitted_at": metadata.get("submitted_at"),
            "started_at": metadata.get("started_at"),
            "completed_at": metadata.get("completed_at")
        }

        if metadata["status"] == JobStatus.FAILED.value:
            result["error"] = metadata.get("error")

        # Add progress info for running jobs
        if metadata["status"] == JobStatus.RUNNING.value:
            log_file = self.jobs_dir / job_id / "job.log"
            if log_file.exists():
                # Get last few lines for progress indication
                with open(log_file) as f:
                    lines = f.readlines()
                    if lines:
                        result["last_log_line"] = lines[-1].strip()

        return result

    def get_job_result(self, job_id: str) -> Dict[str, Any]:
        """Get results of a completed job."""
        metadata = self._load_metadata(job_id)
        if not metadata:
            return {"status": "error", "error": f"Job {job_id} not found"}

        if metadata["status"] != JobStatus.COMPLETED.value:
            return {
                "status": "error",
                "error": f"Job not completed. Current status: {metadata['status']}"
            }

        # Load output
        job_dir = self.jobs_dir / job_id
        output_file = job_dir / "output.json"

        if output_file.exists():
            with open(output_file) as f:
                content = f.read().strip()

                # Try to parse as JSON first
                try:
                    result = json.loads(content)
                    return {"status": "success", "result": result}
                except json.JSONDecodeError:
                    # If JSON parsing fails, treat as CSV and convert
                    import csv
                    import io

                    try:
                        # Parse CSV content
                        reader = csv.DictReader(io.StringIO(content))
                        data = list(reader)

                        # Convert to structured result
                        result = {
                            "format": "csv",
                            "data": data,
                            "total_records": len(data),
                            "output_file": str(output_file)
                        }

                        # Add summary if it's scoring data
                        if data and 'affinity' in data[0]:
                            affinities = [float(row['affinity']) for row in data if row['affinity'] and row['affinity'] != 'None']
                            if affinities:
                                result["summary"] = {
                                    "best_affinity": min(affinities),
                                    "mean_affinity": sum(affinities) / len(affinities),
                                    "total_compounds": len(affinities)
                                }

                        return {"status": "success", "result": result}
                    except Exception as e:
                        # If CSV parsing also fails, return raw content
                        return {
                            "status": "success",
                            "result": {
                                "format": "raw",
                                "content": content,
                                "output_file": str(output_file)
                            }
                        }
        else:
            return {"status": "error", "error": "Output file not found"}

    def get_job_log(self, job_id: str, tail: int = 50) -> Dict[str, Any]:
        """Get log output from a job."""
        job_dir = self.jobs_dir / job_id
        log_file = job_dir / "job.log"

        if not log_file.exists():
            return {"status": "error", "error": f"Log not found for job {job_id}"}

        with open(log_file) as f:
            lines = f.readlines()

        return {
            "status": "success",
            "job_id": job_id,
            "log_lines": lines[-tail:] if tail else lines,
            "total_lines": len(lines)
        }

    def cancel_job(self, job_id: str) -> Dict[str, Any]:
        """Cancel a running job."""
        with self._lock:
            if job_id in self._running_jobs:
                self._running_jobs[job_id].terminate()
                metadata = self._load_metadata(job_id)
                metadata["status"] = JobStatus.CANCELLED.value
                metadata["completed_at"] = datetime.now().isoformat()
                self._save_metadata(job_id, metadata)
                return {"status": "success", "message": f"Job {job_id} cancelled"}

        return {"status": "error", "error": f"Job {job_id} not running"}

    def list_jobs(self, status: Optional[str] = None) -> Dict[str, Any]:
        """List all jobs, optionally filtered by status."""
        jobs = []
        for job_dir in self.jobs_dir.iterdir():
            if job_dir.is_dir():
                metadata = self._load_metadata(job_dir.name)
                if metadata:
                    if status is None or metadata["status"] == status:
                        jobs.append({
                            "job_id": metadata["job_id"],
                            "job_name": metadata.get("job_name"),
                            "status": metadata["status"],
                            "submitted_at": metadata.get("submitted_at"),
                            "script": metadata.get("script")
                        })

        # Sort by submission time (most recent first)
        jobs.sort(key=lambda x: x.get("submitted_at", ""), reverse=True)

        return {"status": "success", "jobs": jobs, "total": len(jobs)}

    def cleanup_old_jobs(self, older_than_days: int = 7) -> Dict[str, Any]:
        """Clean up completed jobs older than specified days."""
        cutoff = datetime.now()
        cleaned = 0

        for job_dir in self.jobs_dir.iterdir():
            if job_dir.is_dir():
                metadata = self._load_metadata(job_dir.name)
                if metadata and metadata["status"] in [JobStatus.COMPLETED.value, JobStatus.FAILED.value]:
                    completed_at = metadata.get("completed_at")
                    if completed_at:
                        completed_time = datetime.fromisoformat(completed_at)
                        if (cutoff - completed_time).days > older_than_days:
                            # Remove job directory
                            import shutil
                            shutil.rmtree(job_dir)
                            cleaned += 1

        return {"status": "success", "cleaned_jobs": cleaned}

    def _save_metadata(self, job_id: str, metadata: Dict):
        """Save job metadata to disk with file locking."""
        meta_file = self.jobs_dir / job_id / "metadata.json"
        meta_file.parent.mkdir(parents=True, exist_ok=True)

        # Atomic write with file locking
        temp_file = meta_file.with_suffix('.tmp')
        with open(temp_file, 'w') as f:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            json.dump(metadata, f, indent=2)
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)

        temp_file.replace(meta_file)

    def _load_metadata(self, job_id: str) -> Optional[Dict]:
        """Load job metadata from disk."""
        meta_file = self.jobs_dir / job_id / "metadata.json"
        if meta_file.exists():
            try:
                with open(meta_file) as f:
                    fcntl.flock(f.fileno(), fcntl.LOCK_SH)
                    data = json.load(f)
                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)
                    return data
            except (json.JSONDecodeError, IOError) as e:
                logger.error(f"Error loading metadata for job {job_id}: {e}")
                return None
        return None

    def get_queue_info(self) -> Dict[str, Any]:
        """Get information about the job queue."""
        with self._lock:
            running_count = len(self._running_jobs)

        all_jobs = self.list_jobs()["jobs"]
        pending_jobs = [j for j in all_jobs if j["status"] == JobStatus.PENDING.value]

        return {
            "status": "success",
            "running_jobs": running_count,
            "pending_jobs": len(pending_jobs),
            "total_jobs": len(all_jobs),
            "queue_info": "Single-threaded execution - one job at a time"
        }

# Global job manager instance
job_manager = JobManager()