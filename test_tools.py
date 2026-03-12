#!/usr/bin/env python3
"""Proper MCP tool testing with correct result handling."""

import sys
import asyncio
import json
from pathlib import Path

# Add src to path
sys.path.append('src')
from server import mcp

def extract_result(tool_result):
    """Extract JSON data from ToolResult."""
    if hasattr(tool_result, 'content') and tool_result.content:
        text_content = tool_result.content[0].text
        try:
            return json.loads(text_content)
        except json.JSONDecodeError:
            return {"error": "Failed to parse JSON", "raw": text_content}
    return {"error": "No content found"}

async def test_server_info():
    """Test server info."""
    print("=== Testing get_server_info ===")
    try:
        result = await mcp.call_tool("get_server_info", {})
        data = extract_result(result)

        if data.get("status") == "success":
            print(f"✓ Server: {data['server_name']} v{data['version']}")
            print(f"  Sync tools: {len(data['sync_tools'])}")
            print(f"  Submit tools: {len(data['submit_tools'])}")
            print(f"  Job management: {len(data['job_management'])}")
            return True
        else:
            print(f"✗ Server info failed: {data}")
            return False
    except Exception as e:
        print(f"✗ Server info failed: {e}")
        return False

async def test_molecular_analysis():
    """Test molecular analysis."""
    print("\n=== Testing analyze_molecules ===")
    try:
        ligand_file = "examples/data/184l_lig.sdf"
        if not Path(ligand_file).exists():
            print(f"✗ Test skipped - file not found: {ligand_file}")
            return False

        result = await mcp.call_tool("analyze_molecules", {
            "ligand_file": ligand_file,
            "descriptors": ["molecular_weight", "logp", "tpsa"]
        })
        data = extract_result(result)

        if data.get("status") == "success":
            print(f"✓ Molecular analysis completed")
            if "total_molecules" in data:
                print(f"  Analyzed molecules: {data['total_molecules']}")
            return True
        else:
            print(f"✗ Molecular analysis failed: {data.get('error', 'Unknown error')}")
            return False
    except Exception as e:
        print(f"✗ Molecular analysis failed: {e}")
        return False

async def test_protein_ligand_scoring():
    """Test protein-ligand scoring."""
    print("\n=== Testing score_protein_ligand ===")
    try:
        receptor_file = "examples/data/184l_rec_ref.pdb"
        ligand_file = "examples/data/184l_lig.sdf"

        if not Path(receptor_file).exists() or not Path(ligand_file).exists():
            print(f"✗ Test skipped - files not found")
            print(f"  Receptor: {Path(receptor_file).exists()}")
            print(f"  Ligand: {Path(ligand_file).exists()}")
            return False

        result = await mcp.call_tool("score_protein_ligand", {
            "receptor_file": receptor_file,
            "ligand_file": ligand_file,
            "scoring_function": "default"
        })
        data = extract_result(result)

        if data.get("status") == "success":
            print(f"✓ Protein-ligand scoring completed")
            return True
        else:
            print(f"✗ Scoring failed: {data.get('error', 'Unknown error')}")
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
        data = extract_result(result)

        if data.get("status") == "submitted":
            job_id = data.get("job_id")
            print(f"✓ Job submitted with ID: {job_id}")

            # Test job status
            status_result = await mcp.call_tool("get_job_status", {"job_id": job_id})
            status_data = extract_result(status_result)
            print(f"✓ Job status: {status_data.get('status', 'unknown')}")
            return True
        else:
            print(f"✗ Job submission failed: {data.get('error', 'Unknown error')}")
            return False
    except Exception as e:
        print(f"✗ Job submission failed: {e}")
        return False

async def test_list_jobs():
    """Test list jobs."""
    print("\n=== Testing list_jobs ===")
    try:
        result = await mcp.call_tool("list_jobs", {})
        data = extract_result(result)

        if data.get("status") == "success":
            jobs = data.get("jobs", [])
            print(f"✓ Listed {len(jobs)} jobs")
            for job in jobs[:3]:  # Show first 3 jobs
                print(f"  Job {job.get('job_id')}: {job.get('status')} - {job.get('job_name')}")
            return True
        else:
            print(f"✗ List jobs failed: {data.get('error', 'Unknown error')}")
            return False
    except Exception as e:
        print(f"✗ List jobs failed: {e}")
        return False

async def main():
    """Run all tests."""
    print("=== Gnina MCP Tool Tests ===")
    print(f"Working directory: {Path.cwd()}")
    print()

    tests = [
        ("Server Info", test_server_info),
        ("Molecular Analysis", test_molecular_analysis),
        ("Protein-Ligand Scoring", test_protein_ligand_scoring),
        ("Job Submission", test_job_submission),
        ("List Jobs", test_list_jobs)
    ]

    passed = 0
    total = len(tests)

    for test_name, test_func in tests:
        try:
            if await test_func():
                passed += 1
        except Exception as e:
            print(f"✗ {test_name} failed with exception: {e}")

    print(f"\n=== TEST SUMMARY ===")
    print(f"Passed: {passed}/{total}")
    print(f"Pass rate: {passed/total*100:.1f}%")

    return passed == total

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1)