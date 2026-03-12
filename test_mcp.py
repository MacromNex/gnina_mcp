#!/usr/bin/env python3
"""Test script for Gnina MCP server functionality."""

import sys
import os
import json
from pathlib import Path

# Add src to path
sys.path.insert(0, 'src')

from jobs.manager import job_manager

def test_job_manager():
    """Test basic job manager functionality."""
    print("Testing job manager...")

    # Test job submission
    print("1. Testing job submission")
    result = job_manager.submit_job(
        script_path="scripts/score_protein_ligand.py",
        args={
            "receptor": "examples/data/10gs_rec.pdb",
            "ligand": "examples/data/10gs_lig.sdf"
        },
        job_name="test_scoring_job"
    )

    if result["status"] == "submitted":
        job_id = result["job_id"]
        print(f"✓ Job submitted successfully: {job_id}")

        # Test job status
        print("2. Testing job status check")
        import time
        for i in range(10):  # Wait up to 10 seconds
            status = job_manager.get_job_status(job_id)
            print(f"   Status: {status['status']}")

            if status["status"] == "completed":
                print("✓ Job completed successfully")

                # Test job result
                print("3. Testing job result retrieval")
                result_data = job_manager.get_job_result(job_id)
                if result_data["status"] == "success":
                    print("✓ Job results retrieved successfully")
                    print(f"   Result keys: {list(result_data['result'].keys())}")
                else:
                    print(f"✗ Failed to get results: {result_data}")

                break
            elif status["status"] == "failed":
                print(f"✗ Job failed: {status.get('error', 'Unknown error')}")
                # Get logs for debugging
                log_result = job_manager.get_job_log(job_id, tail=10)
                if log_result["status"] == "success":
                    print("Last log lines:")
                    for line in log_result["log_lines"]:
                        print(f"   {line.strip()}")
                break

            time.sleep(1)
        else:
            print(f"✗ Job did not complete within timeout")

    else:
        print(f"✗ Job submission failed: {result}")

    # Test job listing
    print("4. Testing job listing")
    jobs = job_manager.list_jobs()
    if jobs["status"] == "success":
        print(f"✓ Found {jobs['total']} jobs")
        for job in jobs["jobs"]:
            print(f"   Job {job['job_id']}: {job['status']} - {job.get('job_name', 'unnamed')}")
    else:
        print(f"✗ Failed to list jobs: {jobs}")

    print("\nJob manager testing completed!")

def test_sync_tools():
    """Test synchronous tools."""
    print("\nTesting synchronous tools...")

    # Set PATH to include current directory for gnina mock
    os.environ['PATH'] = f".:{os.environ.get('PATH', '')}"

    try:
        from server import score_protein_ligand, analyze_molecules

        # Test scoring tool
        print("1. Testing score_protein_ligand")
        result = score_protein_ligand(
            receptor_file="examples/data/10gs_rec.pdb",
            ligand_file="examples/data/10gs_lig.sdf"
        )

        if result["status"] == "success":
            print("✓ Scoring tool works successfully")
            print(f"   Found {len(result['results'])} scoring results")
        else:
            print(f"✗ Scoring tool failed: {result['error']}")

    except ImportError as e:
        print(f"✗ Import error: {e}")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")

    print("Synchronous tools testing completed!")

if __name__ == "__main__":
    # Set PATH to include mock gnina
    os.environ['PATH'] = f".:{os.environ.get('PATH', '')}"

    test_job_manager()
    test_sync_tools()