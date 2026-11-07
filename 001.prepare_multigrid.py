#!/usr/bin/env python3
"""Prepare or generate one energy grid per receptor residue."""

from pathlib import Path
import array
import csv
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "inputs"
SOFTWARE = ROOT / "software"
GRIDS = ROOT / "grids"
SELECTION = ROOT / "selection"
LOGS = ROOT / "logs"

def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()

def load_residues():
    # This mapping was obtained directly from the receptor atom records.
    with (SELECTION / "receptor_residue_map.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    positions = [int(row["footprint_position"]) for row in rows]
    ids = [int(row["mol2_residue_id"]) for row in rows]

    if positions != list(range(1, len(rows) + 1)) or positions != ids:
        raise ValueError("Expected matching sequential footprint and MOL2 IDs.")
    if len(rows) < 2:
        raise ValueError("At least two receptor residues are required.")
    return rows

def components(rows):
    result = []
    for row in rows[:-1]:
        resid = int(row["mol2_residue_id"])
        result.append((row, f"resid_{resid:03d}"))

    # A one-residue remainder preserves the final residue as its own component.
    result.append((rows[-1], "resid_remaining"))
    return result

def prepare():
    rows = load_residues()

    if list(GRIDS.glob("*.nrg")) or list(GRIDS.glob("*.bmp")):
        raise ValueError("Existing grids found; preparation would overwrite their setup.")

    # Obtain the native utilities from the installed DOCK distribution.
    dock = Path(os.environ.get(
        "DOCK_ROOT",
        "/gpfs/projects/rizzo/iamanor/DOCK6_Development/Similarity_Island/dock6_dev",
    ))
    names = ("grid", "multigrid_fp_gen.py", "mol2.py")
    missing = [name for name in names if not (SOFTWARE / name).is_file()]

    for name in missing:
        if not (dock / "bin" / name).is_file():
            raise FileNotFoundError(dock / "bin" / name)

    SOFTWARE.mkdir(exist_ok=True)
    for name in missing:
        shutil.copy2(dock / "bin" / name, SOFTWARE / name)

    for directory in (GRIDS, SELECTION, LOGS):
        directory.mkdir(exist_ok=True)

    required = [
        INPUTS / "receptor.mol2",
        INPUTS / "box.pdb",
        INPUTS / "vdw_AMBER_parm99.defn",
        INPUTS / "chem.defn",
        SOFTWARE / "grid",
        SOFTWARE / "multigrid_fp_gen.py",
        SOFTWARE / "mol2.py",
    ]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    # These settings match the original 2ITY grid-generation input.
    # The native wrapper requires the temporary filenames used below.
    grid_input = f"""allow_non_integral_charges yes
compute_grids yes
grid_spacing 0.3
output_molecule yes
contact_score no
chemical_score no
energy_score yes
energy_cutoff_distance 9999
atom_model a
attractive_exponent 6
repulsive_exponent 9
distance_dielectric yes
dielectric_factor 4
bump_filter yes
bump_overlap 0.75
receptor_file temp.mol2
box_file ../inputs/box.pdb
vdw_definition_file ../inputs/vdw_AMBER_parm99.defn
chemical_definition_file ../inputs/chem.defn
score_grid_prefix temp.rec
receptor_out_file temp.rec.grid.mol2
"""
    # Reject filenames that would overflow the legacy GRID buffers.
    filename_keys = {
        "receptor_file", "box_file", "vdw_definition_file",
        "chemical_definition_file", "score_grid_prefix", "receptor_out_file",
    }
    for line in grid_input.splitlines():
        fields = line.split()
        if fields and fields[0] in filename_keys:
            if len(fields) != 2 or len(fields[1].encode()) > 80:
                raise ValueError(f"Unsafe GRID filename: {line}")
    (GRIDS / "grid.in").write_text(grid_input)

    # Residue fragments can have fractional charges after peptide bonds are cut.
    # Allowing those charges preserves the supplied atomic partial charges.
    selected = [row["mol2_residue_id"] for row in rows[:-1]]
    (SELECTION / "selected_residue_ids.txt").write_text(
        "\n".join(selected) + "\n"
    )

    with (SELECTION / "grid_components.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow([
            "grid_index", "footprint_position", "mol2_residue_id",
            "residue_name", "grid_prefix",
        ])
        for index, (row, prefix) in enumerate(components(rows)):
            writer.writerow([
                index, row["footprint_position"], row["mol2_residue_id"],
                row["residue_name"], str(GRIDS / prefix),
            ])

    # Modify only a local working copy of the native wrapper.
    # Checked shell calls stop generation immediately if GRID or a rename fails.
    original = (SOFTWARE / "multigrid_fp_gen.py").read_text()
    pattern = r"(?m)^([ \t]*)gridpath[ \t]*=.*$"
    if len(re.findall(pattern, original)) != 1:
        raise ValueError("Unexpected native generator gridpath assignment.")
    if "os.system(" not in original:
        raise ValueError("Unexpected native generator command structure.")

    checked = re.sub(
        pattern,
        lambda match: (
            match.group(1) + "gridpath = " +
            repr(shlex.quote(str(SOFTWARE / "grid")))
        ),
        original,
    )
    checked = checked.replace("os.system(", "checked_system(")

    helper = '''
import subprocess

def checked_system(command):
    print("Executing:", command, flush=True)
    subprocess.run(command, shell=True, check=True)

'''
    # Insert after the future import, which must remain at the beginning.
    anchor = "from __future__ import print_function"
    if checked.count(anchor) != 1:
        raise ValueError("Expected native generator future import.")
    checked = checked.replace(anchor, anchor + "\n" + helper, 1)

    working = SOFTWARE / "multigrid_fp_gen_checked.py"
    working.write_text(checked)
    compile(checked, str(working), "exec")
    compile(
        (SOFTWARE / "mol2.py").read_text(),
        str(SOFTWARE / "mol2.py"),
        "exec",
    )

    # These parameters are for Similarity Island FPS, not primary rescoring.
    parameters = [
        "similarity_island_fps_type multigrid",
        "similarity_island_fps_compare_type euclidean",
        "similarity_island_fps_score_type vdw_es",
        "similarity_island_fps_vdw_rep_rad_scale 1.0",
        f"similarity_island_fps_number_of_grids {len(rows)}",
    ]
    for index, (_, prefix) in enumerate(components(rows)):
        parameters.append(
            f"similarity_island_fps_grid_prefix{index} {GRIDS / prefix}"
        )
    (ROOT / "multigrid_parameters.in").write_text(
        "\n".join(parameters) + "\n"
    )

    # Freeze the generation inputs so a later run detects accidental edits.
    tracked = required + [
        GRIDS / "grid.in",
        SELECTION / "selected_residue_ids.txt",
        SELECTION / "grid_components.tsv",
        working,
        Path(__file__).resolve(),
    ]
    manifest = {
        str(path.relative_to(ROOT)): digest(path) for path in tracked
    }
    (ROOT / "generation_inputs.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


    print(f"Prepared {len(rows)} residue components.")
    print(f"Individual grids: {len(rows) - 1}")
    print(f"Remainder contains only residue {rows[-1]['mol2_residue_id']}.")
    print("No grids generated.")
    print("Instructions: README.md")

def validate(rows):
    signature = None
    checksums = []

    for row, prefix in components(rows):
        # The native splitter must preserve the expected number of atoms.
        molecule = GRIDS / f"{prefix}.mol2"
        lines = molecule.read_text().splitlines()
        marker = lines.index("@<TRIPOS>MOLECULE")
        count = int(lines[marker + 2].split()[0])
        if count != int(row["atom_count"]):
            raise ValueError(f"Unexpected component atom count: {prefix}")

        energy = GRIDS / f"{prefix}.nrg"
        bump = GRIDS / f"{prefix}.bmp"

        with energy.open("rb") as stream:
            header_bytes = stream.read(16)
        with bump.open("rb") as stream:
            geometry_bytes = stream.read(32)

        if len(header_bytes) != 16 or len(geometry_bytes) != 32:
            raise ValueError(f"Incomplete grid header: {prefix}")

        size, model, attractive, repulsive = struct.unpack("=4i", header_bytes)
        bsize, spacing, ox, oy, oz, nx, ny, nz = struct.unpack(
            "=if3f3i", geometry_bytes
        )
        if (
            size <= 0 or size != bsize or size != nx * ny * nz
            or min(nx, ny, nz) < 4
            or (attractive, repulsive) != (6, 9)
            or not math.isclose(spacing, 0.3, abs_tol=1e-6)
            or not all(math.isfinite(v) for v in (spacing, ox, oy, oz))
            or energy.stat().st_size != 16 + 12 * size
        ):
            raise ValueError(f"Invalid grid dimensions or model: {prefix}")

        current = (model, attractive, repulsive, geometry_bytes)
        if signature is None:
            signature = current
        elif current != signature:
            raise ValueError(f"Grid geometry or energy model differs: {prefix}")

        # Check the three floating-point energy arrays without loading all grids.
        with energy.open("rb") as stream:
            stream.seek(16)
            while True:
                chunk = stream.read(4 * 262144)
                if not chunk:
                    break
                values = array.array("f")
                values.frombytes(chunk)
                if not all(math.isfinite(value) for value in values):
                    raise ValueError(f"Nonfinite energy values: {prefix}")

        for path in (energy, bump):
            checksums.append([
                path.name, path.stat().st_size, digest(path)
            ])

    with (GRIDS / "grid_checksums.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["file", "bytes", "sha256"])
        writer.writerows(checksums)

def generate():
    rows = load_residues()
    manifest = json.loads((ROOT / "generation_inputs.json").read_text())
    for name, expected in manifest.items():
        if digest(ROOT / name) != expected:
            raise ValueError(f"Generation input changed: {name}")

    # Fresh outputs prevent failed native calls from reusing older grid files.
    if list(GRIDS.glob("*.nrg")) or list(GRIDS.glob("*.bmp")):
        raise ValueError("Grid files already exist; refusing to overwrite them.")

    selected = (SELECTION / "selected_residue_ids.txt").read_text().split()
    command = [
        sys.executable, "-u",
        str(SOFTWARE / "multigrid_fp_gen_checked.py"),
        str(INPUTS / "receptor.mol2"),
        "resid", "grid.in", *selected,
    ]
    started = time.perf_counter()
    result = {"command": command, "status": "running"}
    result_file = LOGS / "generation_result.json"
    result_file.write_text(json.dumps(result, indent=2) + "\n")

    try:
        with (LOGS / "grid_generation.log").open("w") as stream:
            subprocess.run(
                command, cwd=GRIDS, stdout=stream,
                stderr=subprocess.STDOUT, check=True,
            )
        validate(rows)
        result["status"] = "generation_and_file_validation_passed"
        print(f"Generated and validated {len(rows)} grid pairs.")
    except Exception as error:
        result["status"] = "failed"
        result["error"] = str(error)
        raise
    finally:
        result["elapsed_seconds"] = time.perf_counter() - started
        result_file.write_text(json.dumps(result, indent=2) + "\n")

if __name__ == "__main__":
    try:
        if sys.argv[1:] == ["prepare"]:
            prepare()
        elif sys.argv[1:] == ["generate"]:
            generate()
        else:
            raise ValueError("Usage: 001.prepare_multigrid.py prepare|generate")
    except Exception as error:
        raise SystemExit(f"ERROR: {error}")
