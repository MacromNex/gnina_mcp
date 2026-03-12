#!/usr/bin/env python3
"""Direct MCP tool testing script."""

import sys
import asyncio
from pathlib import Path

# Add src to path
sys.path.append('src')
from server import mcp

async def test_get_server_info():
    """Test server info tool."""
    print("=== Testing get_server_info ===")
    try:
        result = await mcp.call_tool("get_server_info", {})
        # Handle ToolResult object
        if hasattr(result, 'content'):
            data = result.content
        elif hasattr(result, 'result'):
            data = result.result
        else:
            data = result

        print(f"✓ Server info retrieved")
        print(f"  Result type: {type(data)}")
        print(f"  Content: {str(data)[:200]}...")
        return True
    except Exception as e:
        print(f"✗ Server info failed: {e}")
        return False

async def test_molecular_analysis():
    """Test molecular analysis tool."""
    print("\n=== Testing analyze_molecules ===")
    try:
        ligand_file = "examples/data/184l_lig.sdf"
        if not Path(ligand_file).exists():
            print(f"✗ Test skipped - file not found: {ligand_file}")
            return False

        result = await mcp.call_tool("analyze_molecules", {
            "ligand_file": ligand_file,
            "descriptors": ["molecular_weight", "logp", "tpsa"],
            "output_file": "test_output.csv"
        })

        # Handle ToolResult object
        if hasattr(result, 'content'):
            data = result.content
        elif hasattr(result, 'result'):
            data = result.result
        else:
            data = result

        print(f"✓ Molecular analysis completed")
        print(f"  Result: {str(data)[:200]}...")
        return True
    except Exception as e:
        print(f"✗ Molecular analysis failed: {e}")
        return False

async def test_protein_ligand_scoring():
    """Test protein-ligand scoring tool."""
    print("\n=== Testing score_protein_ligand ===")
    try:
        receptor_file = "examples/data/184l_rec_ref.pdb"
        ligand_file = "examples/data/184l_lig.sdf"

        if not Path(receptor_file).exists() or not Path(ligand_file).exists():
            print(f"✗ Test skipped - files not found")
            return False

        result = await mcp.call_tool("score_protein_ligand", {
            "receptor_file": receptor_file,
            "ligand_file": ligand_file,
            "scoring_function": "default"
        })

        if result.get("status") == "success":
            print(f"✓ Protein-ligand scoring completed")
            return True
        else:
            print(f"✗ Scoring failed: {result.get('error', 'Unknown error')}")
            return False
    except Exception as e:
        print(f"✗ Protein-ligand scoring failed: {e}")
        return False

async def test_job_submission():
    """Test job submission."""
    print("\n=== Testing submit_molecular_docking ===")
    try:
        receptor_file = "examples/data/184l_rec_ref.pdb"
        ligand_file = "examples/data/184l_lig.sdf"

        if not Path(receptor_file).exists() or not Path(ligand_file).exists():
            print(f"✗ Test skipped - files not found")
            return False

        result = await mcp.call_tool("submit_molecular_docking", {
            "receptor_file": receptor_file,
            "ligand_file": ligand_file,
            "num_modes": 3,
            "exhaustiveness": 4,
            "job_name": "test_docking"
        })

        if result.get("status") == "submitted":
            job_id = result.get("job_id")
            print(f"✓ Job submitted with ID: {job_id}")

            # Test job status
            status_result = await mcp.call_tool("get_job_status", {"job_id": job_id})
            print(f"✓ Job status: {status_result.get('status', 'unknown')}")
            return True
        else:
            print(f"✗ Job submission failed: {result.get('error', 'Unknown error')}")
            return False
    except Exception as e:
        print(f"✗ Job submission failed: {e}")
        return False

async def test_list_jobs():
    """Test list jobs functionality."""
    print("\n=== Testing list_jobs ===")
    try:
        result = await mcp.call_tool("list_jobs", {})

        if result.get("status") == "success":
            jobs = result.get("jobs", [])
            print(f"✓ Listed {len(jobs)} jobs")
            for job in jobs[:3]:  # Show first 3 jobs
                print(f"  Job {job.get('job_id')}: {job.get('status')} - {job.get('job_name')}")
            return True
        else:
            print(f"✗ List jobs failed: {result.get('error', 'Unknown error')}")
            return False
    except Exception as e:
        print(f"✗ List jobs failed: {e}")
        return False

async def test_get_queue_info():
    """Test queue info."""
    print("\n=== Testing get_queue_info ===")
    try:
        result = await mcp.call_tool("get_queue_info", {})

        if result.get("status") == "success":
            print(f"✓ Queue info retrieved")
            print(f"  Running jobs: {result.get('running_jobs', 0)}")
            print(f"  Pending jobs: {result.get('pending_jobs', 0)}")
            print(f"  Total jobs: {result.get('total_jobs', 0)}")
            return True
        else:
            print(f"✗ Queue info failed: {result.get('error', 'Unknown error')}")
            return False
    except Exception as e:
        print(f"✗ Queue info failed: {e}")
        return False

async def main():
    """Run all tests."""
    print("=== Gnina MCP Direct Tool Tests ===")
    print(f"Working directory: {Path.cwd()}")

    tests = [
        test_get_server_info,
        test_molecular_analysis,
        test_protein_ligand_scoring,
        test_job_submission,
        test_list_jobs,
        test_get_queue_info
    ]

    passed = 0
    total = len(tests)

    for test in tests:
        try:
            if await test():
                passed += 1
        except Exception as e:
            print(f"✗ Test failed with exception: {e}")

    print(f"\n=== TEST SUMMARY ===")
    print(f"Passed: {passed}/{total}")
    print(f"Pass rate: {passed/total*100:.1f}%")

    return passed == total

if __name__ == "__main__":
    asyncio.run(main())