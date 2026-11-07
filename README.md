# 2ITY MultiGrid FPS

Generate receptor energy grids for Similarity Island footprint clustering
in DOCK6.

## Inputs

- `inputs/`: prepared 2ITY receptor, binding-site box, and parameter files.
- `software/`: GRID executable and native MultiGrid Python utilities.
- `selection/receptor_residue_map.tsv`: receptor residue mapping.
- `input_manifest.tsv`: original input paths and checksums.

The scripts require these prepared inputs.

## Grid settings

Grid spacing: 0.3 Å; all-atom model; 6–9 VDW potential;
distance-dependent dielectric with factor 4; energy cutoff: 9999 Å.

Short relative filenames avoid the legacy GRID filename-buffer limit.

## Residue components

The current script generates 323 components: residues 1–322 individually
and residue 323 alone in the remainder grid.

Custom important-residue selection is not yet supported by this script.

## Run

For a fresh generation:

```bash
python3 001.prepare_multigrid.py prepare
sbatch 002.generate_multigrid.sh
```

To use another partition:

```bash
sbatch --partition=long-28core 002.generate_multigrid.sh
```

Generation refuses to overwrite existing energy or bump grids.
Do not rerun preparation or generation for the completed grid set.

## Outputs

- `grids/`: component grids, receptor fragments, and GRID output.
- `grids/grid_checksums.tsv`: validated grid checksums.
- `multigrid_parameters.in`: Similarity Island grid configuration.
- `generation_inputs.json`: generation-input checksums.
- `logs/generation_result.json`: generation status and elapsed time.

Validation checks component atom counts, grid geometry, energy-file sizes,
and finite energy values. Cartesian footprint agreement and clustering
performance require separate experiments.

## Current result

All 323 grid pairs were generated and passed file validation.
Development helpers and the original generation scripts are preserved
under `archive/development/`.
# Mult-Grid-Generation-in-Dock6
# Mult-Grid-Generation-in-Dock6
