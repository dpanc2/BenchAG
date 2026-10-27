#!/usr/bin/env python3
"""Validate inputs, record provenance, and launch the imported expression workflow."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

def prepare(config, check_only=False):
    import yaml
    import pandas as pd
    os.environ.setdefault("BENCHMARK_ROOT", str(ROOT))
    os.environ.setdefault("DATA_ROOT", str(ROOT / "data"))
    cfg = yaml.safe_load(os.path.expandvars(Path(config).read_text()))
    if cfg.get("organism") not in ("HOMO_SAPIENS", "MUS_MUSCULUS"):
        raise ValueError("Explicit organism is required")
    if (cfg.get("sequence_length"), cfg.get("gene_mask_width")) not in ((1048576, 200001), (16384, 10001)):
        raise ValueError("Supported (sequence_length, gene_mask_width): (1048576,200001), (16384,10001)")
    if cfg.get("prediction_functions") != ["intervals"]:
        raise ValueError("The supported workflow uses intervals only; variants remain an archival ablation")
    if not cfg.get("splits") or not set(cfg["splits"]) <= {"test", "valid"}:
        raise ValueError("splits must contain test and/or valid")
    for key in ("assembly", "annotation_release", "truth_normalization", "split_provenance", "model_version"):
        if not cfg.get(key) or "REQUIRED" in str(cfg[key]):
            raise ValueError(f"Fill in {key}; metadata must describe the actual data")
    input_keys = ["ontology_json", "path_to_selected"]
    # Upstream loads all four interval tables for gene-ID mapping.
    input_keys += ["genes_forward_val", "genes_reverse_val", "genes_forward_test", "genes_reverse_test"]
    input_keys += ["true_" + split + "_file" for split in cfg["splits"]]
    if cfg.get("tissues"):
        input_keys.append("tissues")
    hashes = {}
    for key in input_keys:
        path = Path(cfg[key]).expanduser()
        if not path.is_absolute():
            path = ROOT / path
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"{key}: {path}")
        cfg[key] = str(path)
        hashes[key] = hashlib.sha256(path.read_bytes()).hexdigest()
    mapping = json.loads(Path(cfg["ontology_json"]).read_text())
    if not mapping or not all(isinstance(v, str) for v in mapping.values()):
        raise ValueError("ontology_json must map track names to single sample IDs")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Mapping creates duplicate sample IDs; resolve before scoring")
    split_ids = {}
    prefix = "ENSG" if cfg["organism"] == "HOMO_SAPIENS" else "ENSMUSG"
    for split, suffix in (("test", "test"), ("valid", "val")):
        ids = set()
        for direction, strand in (("forward", "+"), ("reverse", "-")):
            df = pd.read_csv(cfg[f"genes_{direction}_{suffix}"], sep="\t")
            required = {"gene_id", "gene_id_unversioned", "chrom", "TSS", "gene_strand"}
            if not required <= set(df.columns):
                raise ValueError(f"Missing interval columns: {required - set(df.columns)}")
            if df.gene_id_unversioned.duplicated().any() or ids & set(df.gene_id_unversioned):
                raise ValueError(f"Duplicate interval genes in {split}")
            if not df.gene_strand.eq(strand).all() or not df.gene_id_unversioned.str.startswith(prefix).all():
                raise ValueError(f"Organism/strand mismatch in {split} {direction}")
            tss = pd.to_numeric(df.TSS, errors="raise")
            if (tss < 0).any() or (tss % 1 != 0).any():
                raise ValueError("TSS must be a nonnegative integer coordinate")
            ids.update(df.gene_id_unversioned)
        split_ids[split] = ids
        if split in cfg["splits"]:
            truth = pd.read_csv(cfg[f"true_{split}_file"])
            if "gene_id" not in truth or truth.gene_id.duplicated().any():
                raise ValueError("Truth requires unique gene_id")
            if not set(truth.gene_id.str.split(".").str[0]) <= ids:
                raise ValueError(f"Truth genes missing from {split} interval tables")
            if len(set(truth.columns) & set(mapping.values())) < 2:
                raise ValueError("Fewer than two mapped truth tracks; confirm sample ID mapping")
    if split_ids["test"] & split_ids["valid"]:
        raise ValueError("Gene overlap between test and validation tables")
    # Different settings/inputs receive isolated checkpoints and predictions.
    hashes["runner_source"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for source in ("Alphagenome_prediction_benchmarking_json.py", "score_ct_specificity.py"):
        hashes[source] = hashlib.sha256((ROOT / "benchmarks/expression/scripts" / source).read_bytes()).hexdigest()
    try:
        hashes["alphagenome_client"] = importlib.metadata.version("alphagenome")
    except importlib.metadata.PackageNotFoundError:
        hashes["alphagenome_client"] = "not installed"
    identity = hashlib.sha256(json.dumps({"config": cfg, "inputs": hashes}, sort_keys=True).encode()).hexdigest()[:16]
    out = ROOT / "results" / "expression" / identity
    cfg["output_dir"] = str(out)
    cfg["checkpoint_path"] = str(out / "checkpoint_{split}_{prediction_functions_name}.csv")
    cfg["save_pred_path"] = "predictions_{split}_{prediction_functions_name}.csv"
    if check_only:
        print(f"Input checks passed; run ID {identity}. This does not verify model training-region holdouts.")
        return cfg
    if cfg.get("run_prediction", True) and not os.environ.get("ALPHAGENOME_API_KEY"):
        raise ValueError("Set ALPHAGENOME_API_KEY")
    out.mkdir(parents=True, exist_ok=True)
    safe = {k: v for k, v in cfg.items() if k != "api_key"}
    (out / "resolved_config.yaml").write_text(yaml.safe_dump(safe))
    versions = {}
    for name in ("alphagenome", "pandas", "numpy", "biopython", "PyYAML"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not installed"
    git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True)
    (out / "manifest.json").write_text(json.dumps({"git_commit": git.stdout.strip(), "input_sha256": hashes, "packages": versions, "run_id": identity}, indent=2))
    env = os.environ.copy()
    env["ALPHAGENOME_ONTOLOGY_JSON"] = cfg["ontology_json"]
    subprocess.run([sys.executable, str(ROOT / "benchmarks/expression/scripts/Alphagenome_prediction_benchmarking_json.py"), "--config", str(out / "resolved_config.yaml")], cwd=ROOT, env=env, check=True)
    return cfg

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--check", action="store_true", help="Validate external inputs without API calls")
    args = parser.parse_args()
    try:
        prepare(args.config, args.check)
    except (ValueError, FileNotFoundError) as exc:
        parser.exit(2, f"Input error: {exc}\n")
