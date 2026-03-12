# Gnina MCP Integration Test Prompts

## Tool Discovery Tests

### Prompt 1: List All Tools
"What MCP tools are available for molecular docking? Give me a brief description of each."

**Expected Response**: Should list all gnina-tools including:
- Sync tools: score_protein_ligand, analyze_molecules
- Submit tools: submit_molecular_docking, submit_virtual_screening, submit_flexible_docking, submit_cnn_comparison
- Job management: get_job_status, get_job_result, get_job_log, cancel_job, list_jobs, etc.

### Prompt 2: Server Information
"Get information about the gnina-tools MCP server capabilities."

**Expected Response**: Should return server info with version, supported formats, and tool lists.

## Sync Tool Tests

### Prompt 3: Score Protein-Ligand Complex
"Score a protein-ligand complex using gnina. Use receptor file 'examples/data/protein.pdb' and ligand file 'examples/data/ligand.sdf'"

**Expected Response**: Should execute scoring and return affinity scores and CNN model results.

### Prompt 4: Molecular Analysis
"Analyze molecular descriptors for a ligand file 'examples/data/ligand.sdf'. Calculate molecular weight, logP, TPSA, HBD, HBA, and drug-likeness."

**Expected Response**: Should return molecular properties and drug-likeness analysis.

### Prompt 5: Error Handling - Missing File
"Score a protein-ligand complex with receptor 'nonexistent.pdb' and ligand 'missing.sdf'"

**Expected Response**: Should return error with "File not found" message.

## Submit API Tests (Long-running Jobs)

### Prompt 6: Submit Molecular Docking
"Submit a molecular docking job with receptor 'examples/data/protein.pdb' and ligand 'examples/data/ligand.sdf'. Use 9 modes and exhaustiveness of 8."

**Expected Response**:
```json
{
  "status": "submitted",
  "job_id": "abc12345",
  "message": "Job submitted. Use get_job_status('abc12345') to check progress."
}
```

### Prompt 7: Check Job Status
"Check the status of job [insert job_id from previous test]"

**Expected Response**:
```json
{
  "job_id": "abc12345",
  "status": "running|completed|failed|pending",
  "submitted_at": "2024-01-01T12:00:00",
  "started_at": "2024-01-01T12:01:00",
  ...
}
```

### Prompt 8: Get Job Results
"Get the results for completed job [job_id]"

**Expected Response**: Should return docking results with poses, affinities, and output files.

### Prompt 9: View Job Logs
"Show me the last 30 lines of logs for job [job_id]"

**Expected Response**: Should return recent log lines from the docking process.

### Prompt 10: Submit Virtual Screening
"Submit a virtual screening job with receptor 'examples/data/protein.pdb' and ligand directory 'examples/data/ligands/'. Set top_n to 50 and affinity cutoff to -7.0 kcal/mol."

**Expected Response**: Job submission confirmation with job_id.

### Prompt 11: Submit Flexible Docking
"Submit a flexible docking job with receptor 'examples/data/protein.pdb' and ligand 'examples/data/ligand.sdf'. Use flexdist of 3.5 and compare with rigid docking."

**Expected Response**: Job submission confirmation for flexible docking.

## Job Management Tests

### Prompt 12: List All Jobs
"List all submitted molecular docking jobs"

**Expected Response**: Should return list of jobs with their status, submission times, and names.

### Prompt 13: List Completed Jobs Only
"List all jobs with status 'completed'"

**Expected Response**: Should return filtered list of only completed jobs.

### Prompt 14: Get Queue Information
"Get information about the job queue status"

**Expected Response**: Should return queue statistics including running and pending job counts.

### Prompt 15: Cancel Running Job
"Cancel job [running_job_id]"

**Expected Response**: Should confirm cancellation or return error if job not running.

## Batch Processing Tests

### Prompt 16: Multiple Scoring Jobs
"Score these protein-ligand pairs:
1. receptor: 'examples/data/protein1.pdb', ligand: 'examples/data/ligand1.sdf'
2. receptor: 'examples/data/protein2.pdb', ligand: 'examples/data/ligand2.sdf'"

**Expected Response**: Should handle multiple scoring requests either sequentially or in parallel.

### Prompt 17: Screening Multiple Targets
"Submit virtual screening for these targets against the same ligand library:
- Target 1: 'examples/data/protein1.pdb'
- Target 2: 'examples/data/protein2.pdb'
- Ligand directory: 'examples/data/ligands/'"

**Expected Response**: Should submit multiple screening jobs.

## Real-World Scenarios

### Scenario 1: Drug Discovery Pipeline
"I'm working on a drug discovery project targeting EGFR. Please:
1. Score my lead compound 'lead_compound.sdf' against 'egfr_structure.pdb'
2. If the affinity is good (<-8 kcal/mol), submit a flexible docking job
3. While it's running, analyze the molecular properties of the compound
4. Check the job status and get results when complete"

**Expected Flow**:
1. Execute scoring → check affinity
2. Submit flexible docking (if affinity passes threshold)
3. Execute molecular analysis
4. Monitor job and retrieve results

### Scenario 2: Virtual Screening Study
"I want to screen a compound library for new inhibitors:
1. Submit virtual screening of 'compound_library/' against 'target_protein.pdb'
2. Set parameters: top 100 hits, affinity cutoff -7.5 kcal/mol
3. Monitor the job progress
4. When complete, analyze the molecular properties of the top 10 hits"

**Expected Flow**:
1. Submit screening job
2. Monitor progress
3. Retrieve results
4. Analyze top hits

### Scenario 3: CNN Model Comparison
"Compare CNN model performance for protein-ligand scoring:
1. Submit CNN comparison job with receptor 'test_protein.pdb' and ligand 'test_ligand.sdf'
2. Test 3 iterations with different models
3. Monitor the benchmarking process
4. Analyze which model performs best"

**Expected Flow**:
1. Submit CNN comparison
2. Monitor progress
3. Get comparative results

### Scenario 4: Structure-Activity Relationship Study
"For SAR analysis, I need to:
1. Score multiple analogs against the same target
2. Compare their molecular properties
3. Identify the best compounds for follow-up docking"

**Expected Flow**:
1. Batch scoring of analogs
2. Molecular analysis
3. Ranking and selection

## Error Handling Tests

### Prompt 18: Invalid File Format
"Score a complex with receptor 'protein.xyz' (invalid format)"

**Expected Response**: Should return error about unsupported file format.

### Prompt 19: Missing Parameters
"Submit molecular docking without specifying receptor file"

**Expected Response**: Should return error about missing required parameters.

### Prompt 20: Invalid Job ID
"Get status of job 'nonexistent_job_id'"

**Expected Response**: Should return error that job not found.

### Prompt 21: Conflicting Parameters
"Submit docking with both autobox_ligand and explicit center/size coordinates"

**Expected Response**: Should handle parameter conflicts gracefully.

## Performance Tests

### Prompt 22: Large File Handling
"Score a complex with a large protein (>5000 atoms) and multi-conformer ligand"

**Expected Response**: Should handle large files appropriately, potentially recommending submit API.

### Prompt 23: Concurrent Job Submission
"Submit 3 docking jobs simultaneously with different ligands"

**Expected Response**: Should queue jobs appropriately with proper job management.

## Expected Success Criteria

For each test prompt:

1. **Response Time**:
   - Sync tools: < 1 minute for simple operations
   - Submit tools: Immediate job submission with job_id

2. **Data Formats**:
   - Proper JSON responses with status indicators
   - Clear error messages for failures
   - Structured results for successful operations

3. **Job Management**:
   - Unique job IDs for tracking
   - Accurate status updates
   - Complete logs and results

4. **Error Handling**:
   - Graceful failure with helpful error messages
   - No server crashes or hangs
   - Consistent error response format

5. **File Handling**:
   - Proper path resolution
   - File existence validation
   - Format validation for molecular files

## Notes for Testers

- Replace file paths with actual paths to test data
- Monitor system resources during long-running jobs
- Check log files for any unexpected errors
- Verify output files are created and valid
- Test with various molecular file formats (PDB, SDF, PDBQT)