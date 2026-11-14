#!/usr/bin/env python
"""Dock every row of an active/control CSV with AutoDock Vina.

The script is intentionally CLI-only so the exact receptor, box and Vina
settings are recorded. It supports parallel workers and writes progress after
each completed molecule, allowing interrupted runs to resume.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def dock_one(index: int, row: dict, args: argparse.Namespace, run_dir: Path) -> tuple[int, float | None, str]:
    ligand = run_dir / f"ligand_{index}.pdbqt"
    pose = run_dir / f"pose_{index}.pdbqt"
    try:
        subprocess.run(
            ["obabel", f"-:{row['smiles']}", "-O", str(ligand), "--gen3d", "-p", "7.4"],
            check=True, capture_output=True, text=True, timeout=args.prep_timeout,
        )
        cmd = [
            args.vina, "--receptor", str(args.receptor), "--ligand", str(ligand),
            "--center_x", str(args.center[0]), "--center_y", str(args.center[1]),
            "--center_z", str(args.center[2]), "--size_x", str(args.size[0]),
            "--size_y", str(args.size[1]), "--size_z", str(args.size[2]),
            "--exhaustiveness", str(args.exhaustiveness), "--cpu", "1", "--out", str(pose),
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=args.dock_timeout)
        match = re.search(r"VINA RESULT:\s*(-?\d+\.\d+)", pose.read_text())
        return index, float(match.group(1)) if match else None, "ok" if match else "no_score"
    except subprocess.TimeoutExpired:
        return index, None, "timeout"
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "failed").strip().splitlines()
        return index, None, detail[-1][:160] if detail else "failed"


def write_rows(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receptor", type=Path, required=True)
    parser.add_argument("--center", type=float, nargs=3, required=True)
    parser.add_argument("--size", type=float, nargs=3, required=True)
    parser.add_argument("--vina", default="vina")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--prep-timeout", type=int, default=180)
    parser.add_argument("--dock-timeout", type=int, default=900)
    parser.add_argument("--run-dir", type=Path, default=Path.home() / ".aipharma_alk5_benchmark")
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    for field in ("score", "dock_status"):
        if field not in fields:
            fields.append(field)
    args.run_dir.mkdir(parents=True, exist_ok=True)

    pending = [(index, row) for index, row in enumerate(rows) if not row.get("score")]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(dock_one, index, row, args, args.run_dir): index
            for index, row in pending
        }
        completed = 0
        for future in as_completed(futures):
            index, score, status = future.result()
            rows[index]["score"] = "" if score is None else score
            rows[index]["dock_status"] = status
            completed += 1
            print(f"[{completed}/{len(pending)}] row={index + 1} score={score} status={status}", flush=True)
            write_rows(args.output, rows, fields)

    successful = sum(bool(row.get("score")) for row in rows)
    print(f"completed {successful}/{len(rows)} dockings; wrote {args.output}")


if __name__ == "__main__":
    main()
