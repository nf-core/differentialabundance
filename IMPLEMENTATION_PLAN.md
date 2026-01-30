# Implementation Plan for Issue #472

## Problem

Currently paramsheet params override CLI flags and `-params-file`, preventing command-line overrides.

## Solution

Use Nextflow's `session.cliParams` API (available since v25.10.0) for proper parameter priority.

## Steps

1. Bump minimum Nextflow version to 25.10.0
2. Update param merging logic in getParamsheetConfigurations()
3. Update CI and documentation
4. Add minimal test for CLI override

## Files to Modify

- nextflow.config
- subworkflows/local/utils_nfcore_differentialabundance_pipeline/main.nf
- .github/workflows/nf-test.yml
- README.md
- ro-crate-metadata.json
- tests/ (new test file)
