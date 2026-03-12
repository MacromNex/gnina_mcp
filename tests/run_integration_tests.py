#!/usr/bin/env python3
"""Automated integration test runner for Gnina MCP server."""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

class MCPTestRunner:
    def __init__(self, server_path: str):
        self.server_path = Path(server_path)
        self.results = {
            "test_date": datetime.now().isoformat(),
            "server_name": "gnina-tools",
            "server_path": str(server_path),
            "tests": {},
            "issues": [],
            "summary": {}
        }

    def test_server_startup(self) -> bool:
        """Test that server starts without errors."""
        print("Testing server startup...")
        try:
            # Test import in the correct directory
            result = subprocess.run(
                ["python", "-c", "from server import mcp; print('Server imported successfully')"],
                capture_output=True, text=True, timeout=30,
                cwd=self.server_path.parent / "src"
            )
            success = result.returncode == 0
            self.results["tests"]["server_startup"] = {
                "status": "passed" if success else "failed",
                "output": result.stdout,
                "error": result.stderr if result.stderr else None
            }
            print(f"✓ Server startup: {'PASSED' if success else 'FAILED'}")
            return success
        except Exception as e:
            self.results["tests"]["server_startup"] = {"status": "error", "error": str(e)}
            print(f"✗ Server startup: ERROR - {e}")
            return False

    def test_rdkit_import(self) -> bool:
        """Test that RDKit is available."""
        print("Testing RDKit availability...")
        try:
            result = subprocess.run(
                ["python", "-c", "from rdkit import Chem; print('RDKit OK')"],
                capture_output=True, text=True, timeout=30
            )
            success = result.returncode == 0
            self.results["tests"]["rdkit_import"] = {
                "status": "passed" if success else "failed",
                "output": result.stdout,
                "error": result.stderr if result.stderr else None
            }
            print(f"✓ RDKit import: {'PASSED' if success else 'FAILED'}")
            return success
        except Exception as e:
            self.results["tests"]["rdkit_import"] = {"status": "error", "error": str(e)}
            print(f"✗ RDKit import: ERROR - {e}")
            return False

    def test_claude_mcp_registration(self) -> bool:
        """Test that MCP server is registered with Claude Code."""
        print("Testing Claude MCP registration...")
        try:
            result = subprocess.run(
                ["claude", "mcp", "list"],
                capture_output=True, text=True, timeout=30
            )
            # Check if gnina-tools is mentioned in the config
            config_check = subprocess.run(
                ["cat", str(Path.home() / ".claude.json")],
                capture_output=True, text=True
            )
            success = "gnina-tools" in config_check.stdout
            self.results["tests"]["claude_mcp_registration"] = {
                "status": "passed" if success else "failed",
                "output": f"Found gnina-tools in config: {success}",
                "error": None
            }
            print(f"✓ Claude MCP registration: {'PASSED' if success else 'FAILED'}")
            return success
        except Exception as e:
            self.results["tests"]["claude_mcp_registration"] = {"status": "error", "error": str(e)}
            print(f"✗ Claude MCP registration: ERROR - {e}")
            return False

    def test_script_files_exist(self) -> bool:
        """Test that required script files exist."""
        print("Testing script files existence...")
        scripts_dir = self.server_path.parent / "scripts"
        required_scripts = [
            "dock_ligand.py",
            "virtual_screening.py",
            "flexible_docking.py",
            "compare_cnn_models.py",
            "score_protein_ligand.py",
            "molecular_analysis.py"
        ]

        missing_scripts = []
        for script in required_scripts:
            if not (scripts_dir / script).exists():
                missing_scripts.append(script)

        success = len(missing_scripts) == 0
        self.results["tests"]["script_files_exist"] = {
            "status": "passed" if success else "failed",
            "output": f"Checked {len(required_scripts)} scripts",
            "error": f"Missing scripts: {missing_scripts}" if missing_scripts else None
        }
        print(f"✓ Script files: {'PASSED' if success else 'FAILED'}")
        if missing_scripts:
            print(f"  Missing: {', '.join(missing_scripts)}")
        return success

    def test_job_manager_import(self) -> bool:
        """Test that job manager can be imported."""
        print("Testing job manager import...")
        try:
            result = subprocess.run(
                ["python", "-c", "from jobs.manager import job_manager; print('Job manager imported successfully')"],
                capture_output=True, text=True, timeout=30,
                cwd=self.server_path.parent / "src"
            )
            success = result.returncode == 0
            self.results["tests"]["job_manager_import"] = {
                "status": "passed" if success else "failed",
                "output": result.stdout,
                "error": result.stderr if result.stderr else None
            }
            print(f"✓ Job manager import: {'PASSED' if success else 'FAILED'}")
            return success
        except Exception as e:
            self.results["tests"]["job_manager_import"] = {"status": "error", "error": str(e)}
            print(f"✗ Job manager import: ERROR - {e}")
            return False

    def test_environment_activation(self) -> bool:
        """Test that the conda environment can be activated."""
        print("Testing environment activation...")
        try:
            env_path = self.server_path.parent / "env"
            if env_path.exists():
                # Test if environment has the required packages
                result = subprocess.run(
                    ["conda", "run", "-p", str(env_path), "python", "-c", "import sys; print(sys.executable)"],
                    capture_output=True, text=True, timeout=30
                )
                success = result.returncode == 0 and str(env_path) in result.stdout
                self.results["tests"]["environment_activation"] = {
                    "status": "passed" if success else "failed",
                    "output": result.stdout,
                    "error": result.stderr if result.stderr else None
                }
                print(f"✓ Environment activation: {'PASSED' if success else 'FAILED'}")
                return success
            else:
                self.results["tests"]["environment_activation"] = {
                    "status": "failed",
                    "error": "Environment directory not found"
                }
                print("✗ Environment activation: FAILED - env directory not found")
                return False
        except Exception as e:
            self.results["tests"]["environment_activation"] = {"status": "error", "error": str(e)}
            print(f"✗ Environment activation: ERROR - {e}")
            return False

    def generate_report(self) -> str:
        """Generate JSON report."""
        total = len(self.results["tests"])
        passed = sum(1 for t in self.results["tests"].values() if t.get("status") == "passed")
        failed = sum(1 for t in self.results["tests"].values() if t.get("status") == "failed")
        errors = sum(1 for t in self.results["tests"].values() if t.get("status") == "error")

        self.results["summary"] = {
            "total_tests": total,
            "passed": passed,
            "failed": failed,
            "errors": errors,
            "pass_rate": f"{passed/total*100:.1f}%" if total > 0 else "N/A"
        }

        print(f"\n=== TEST SUMMARY ===")
        print(f"Total tests: {total}")
        print(f"Passed: {passed}")
        print(f"Failed: {failed}")
        print(f"Errors: {errors}")
        print(f"Pass rate: {self.results['summary']['pass_rate']}")

        return json.dumps(self.results, indent=2)

if __name__ == "__main__":
    print("=== Gnina MCP Integration Tests ===")
    print(f"Test started at: {datetime.now().isoformat()}")
    print()

    runner = MCPTestRunner("src/server.py")

    # Run all tests
    runner.test_server_startup()
    runner.test_rdkit_import()
    runner.test_claude_mcp_registration()
    runner.test_script_files_exist()
    runner.test_job_manager_import()
    runner.test_environment_activation()

    # Generate report
    report = runner.generate_report()

    # Save report
    report_path = Path("reports/step7_integration_test.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report)

    print(f"\nReport saved to: {report_path}")