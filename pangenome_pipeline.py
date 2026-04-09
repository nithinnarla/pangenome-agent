#!/usr/bin/env python3
"""
Pangenome Pipeline 4 — VCF Region Extraction + Human-Readable Parser
Lucas Borges dos Santos + Nithin Narla — April 2026

Wraps bash steps 1-10 in Python with checkpoint validation at every step.
Uses Docker containers (Mac-compatible). Singularity not required.

Usage:
  python pangenome_pipeline.py --gfa chr6.gfa --vcf chr6.vcf.gz \
      --region "MM26#0#Hg_chrom6_MM26:520000-570000"

  python pangenome_pipeline.py --gfa chr6.gfa --vcf chr6.vcf.gz \
      --start 520000 --end 570000

Options:
  --reset        Clear all checkpoints and run from scratch
  --skip-step N  Skip to step N if earlier steps keep erroring

Current inputs: GFA (.gfa / .gfa.gz) + VCF (.vcf / .vcf.gz)
Future inputs:  GFF and FASTA support will be added when available from Lucas.
"""

import os
import sys
import re
import gzip
import shutil
import subprocess
import argparse
from pathlib import Path
from datetime import datetime

# ── Terminal colors ──────────────────────────────────────────────────────────
GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"

def log(msg, color=RESET):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"{color}[{ts}] {msg}{RESET}", flush=True)

def ok(msg):   log(f"✅  {msg}", GREEN)
def err(msg):  log(f"❌  {msg}", RED)
def warn(msg): log(f"⚠️   {msg}", YELLOW)
def info(msg): log(f"ℹ️   {msg}", CYAN)

# ── Checkpoint system ────────────────────────────────────────────────────────
CHECKPOINT_FILE = "pipeline_checkpoint.txt"

def load_checkpoints():
    if not os.path.exists(CHECKPOINT_FILE):
        return set()
    with open(CHECKPOINT_FILE) as f:
        return set(line.strip() for line in f if line.strip())

def save_checkpoint(step_name):
    with open(CHECKPOINT_FILE, "a") as f:
        f.write(step_name + "\n")
    ok(f"Checkpoint saved — {step_name}")

def clear_checkpoints():
    if os.path.exists(CHECKPOINT_FILE):
        os.remove(CHECKPOINT_FILE)
    node_range = "node_range.txt"
    if os.path.exists(node_range):
        os.remove(node_range)
    info("Checkpoints cleared — pipeline will run from scratch.")

# ── Docker helpers ───────────────────────────────────────────────────────────
def docker_run(image, command, workdir, capture_output=False):
    """Run a command in a Docker container with /data mounted to workdir."""
    cmd = ["docker", "run", "--rm"]
    if image == "odgi":
        cmd += ["--platform", "linux/amd64"]
    cmd += ["-v", f"{workdir}:/data", "-w", "/data", image] + command.split()
    info(f"Docker: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=capture_output, timeout=900)
        return result
    except subprocess.TimeoutExpired:
        err("Docker command timed out after 15 minutes.")
        return None
    except FileNotFoundError:
        err("Docker not found — is Docker Desktop running?")
        return None

def check_docker():
    result = subprocess.run(["docker", "info"], capture_output=True, text=True)
    if result.returncode != 0:
        err("Docker is not running. Start Docker Desktop and retry.")
        sys.exit(1)
    ok("Docker is running.")

def check_images():
    result = subprocess.run(["docker", "images", "--format", "{{.Repository}}"],
                            capture_output=True, text=True)
    images = result.stdout.strip().split("\n")
    missing = [x for x in ["odgi", "vg"] if x not in images]
    if missing:
        err(f"Missing Docker images: {missing}")
        print("\nRun these to install:")
        print("  docker pull pangenome/odgi:1770740056")
        print("  docker pull quay.io/vgteam/vg:v1.73.0")
        print("  docker tag pangenome/odgi:1770740056 odgi")
        print("  docker tag quay.io/vgteam/vg:v1.73.0 vg\n")
        sys.exit(1)
    ok("Docker images 'odgi' and 'vg' confirmed.")

# ── Step helpers ─────────────────────────────────────────────────────────────
def already_done(step):
    return step in load_checkpoints()

def abort(step_num, msg, tip=None):
    err(msg)
    if tip:
        warn(tip)
    warn(f"To resume from next step: re-run with --skip-step {step_num + 1}")
    sys.exit(1)

# ── Auto-detect region from GFA ──────────────────────────────────────────────
def extract_reference_name(gfa_10_path):
    """
    Read P lines from GFA 1.0 and return the first reference name.
    After vg converts 1.1 -> 1.0, P lines appear:
      P  MM26#0#Hg_chrom6_MM26  ...
    We grab column 2 from the first P line.
    """
    try:
        with open(gfa_10_path, "rb") as f:
            for raw in f:
                try:
                    line = raw.decode("utf-8", errors="ignore")
                except Exception:
                    continue
                if line.startswith("P\t"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        return parts[1]
    except Exception as e:
        err(f"Could not read GFA for reference name: {e}")
    return None

def build_region(gfa_10_path, start, end):
    """Build region string like MM26#0#Hg_chrom6_MM26:520000-570000"""
    ref = extract_reference_name(gfa_10_path)
    if ref is None:
        abort(0, "Could not auto-detect reference name from GFA 1.0.",
              "Use --region manually instead.")
    region = f"{ref}:{start}-{end}"
    ok(f"Auto-detected region: {region}")
    return region

# ────────────────────────────────────────────────────────────────────────────
# STEP 1 — Resolve absolute paths
# ────────────────────────────────────────────────────────────────────────────
def step1(gfa_path, vcf_path):
    STEP = "step1_paths"
    if already_done(STEP):
        ok("Step 1 — skipping (already done)")
        return str(Path(gfa_path).resolve()), str(Path(vcf_path).resolve())

    info("── STEP 1: Resolving absolute paths ──────────────────────────")
    gfa_abs = str(Path(gfa_path).resolve())
    vcf_abs = str(Path(vcf_path).resolve())

    if not os.path.exists(gfa_abs):
        abort(1, f"GFA file not found: {gfa_abs}", "Check the --gfa path argument.")
    if not os.path.exists(vcf_abs):
        abort(1, f"VCF file not found: {vcf_abs}", "Check the --vcf path argument.")

    ok(f"GFA: {gfa_abs}")
    ok(f"VCF: {vcf_abs}")
    save_checkpoint(STEP)
    return gfa_abs, vcf_abs

# ────────────────────────────────────────────────────────────────────────────
# STEP 2 — Unzip .gz files if needed
# ────────────────────────────────────────────────────────────────────────────
def step2(gfa_path, vcf_path, workdir):
    STEP = "step2_unzip"
    gfa_out = gfa_path.replace(".gz", "") if gfa_path.endswith(".gz") else gfa_path

    if already_done(STEP):
        ok("Step 2 — skipping (already done)")
        return gfa_out, vcf_path

    info("── STEP 2: Checking compression ──────────────────────────────")

    if gfa_path.endswith(".gz"):
        gfa_out = os.path.join(workdir, Path(gfa_path).stem)
        info(f"Unzipping GFA → {gfa_out}")
        try:
            with gzip.open(gfa_path, "rb") as fi, open(gfa_out, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            ok(f"GFA unzipped: {gfa_out}")
        except Exception as e:
            abort(2, f"Failed to unzip GFA: {e}", "Verify the .gz file is not corrupted.")
    else:
        ok("GFA not compressed — no action needed.")

    if vcf_path.endswith(".gz"):
        ok("VCF is .gz — will be handled by Step 9 (no unzip needed).")
    else:
        ok("VCF not compressed — no action needed.")

    save_checkpoint(STEP)
    return gfa_out, vcf_path

# ────────────────────────────────────────────────────────────────────────────
# STEP 3 — Check GFA version; convert 1.1 → 1.0 if needed
# ────────────────────────────────────────────────────────────────────────────
def step3(gfa_path, workdir):
    STEP = "step3_version"
    stem = Path(gfa_path).stem
    gfa_10_path = os.path.join(workdir, stem + ".10.gfa")

    if already_done(STEP):
        ok("Step 3 — skipping (already done)")
        return gfa_10_path if os.path.exists(gfa_10_path) else gfa_path

    info("── STEP 3: Checking GFA version ──────────────────────────────")

    version = None
    try:
        with open(gfa_path) as f:
            for line in f:
                if line.startswith("H"):
                    if "VN:Z:1.1" in line:
                        version = "1.1"
                    elif "VN:Z:1.0" in line:
                        version = "1.0"
                    break
    except Exception as e:
        abort(3, f"Cannot read GFA: {e}")

    if version is None:
        warn("GFA version header not found — assuming 1.1 and converting.")
        version = "1.1"

    if version == "1.0":
        ok("GFA is already version 1.0 — no conversion needed.")
        dst = os.path.join(workdir, Path(gfa_path).name)
        if not os.path.exists(dst):
            shutil.copy2(gfa_path, dst)
        save_checkpoint(STEP)
        return dst

    info("GFA is version 1.1 — converting to 1.0 with vg convert...")
    gfa_filename = Path(gfa_path).name
    result = docker_run("vg", f"vg convert -gfW {gfa_filename}", workdir, capture_output=True)

    if result is None or result.returncode != 0:
        stderr = result.stderr.decode() if result and isinstance(result.stderr, bytes) else str(result.stderr if result else "N/A")
        abort(3, f"vg convert failed.\n  stderr: {stderr}",
              "Verify Docker image 'vg' is available and GFA file is valid.")

    stdout_data = result.stdout if isinstance(result.stdout, bytes) else result.stdout.encode()
    with open(gfa_10_path, "wb") as f:
        f.write(stdout_data)

    if os.path.getsize(gfa_10_path) == 0:
        abort(3, "Converted GFA is empty.", "The vg convert command may have failed silently.")

    ok(f"GFA converted to 1.0: {gfa_10_path}")
    save_checkpoint(STEP)
    return gfa_10_path

# ────────────────────────────────────────────────────────────────────────────
# STEP 4 — Convert GFA 1.0 → ODGI (.og)
# ────────────────────────────────────────────────────────────────────────────
def step4(gfa_10_path, workdir):
    STEP = "step4_build_og"
    og_name = Path(gfa_10_path).stem + ".og"
    og_path = os.path.join(workdir, og_name)

    if already_done(STEP):
        ok("Step 4 — skipping (already done)")
        return og_path

    info("── STEP 4: Building ODGI graph (.og) ─────────────────────────")
    gfa_name = Path(gfa_10_path).name

    result = docker_run("odgi", f"build -g {gfa_name} -o {og_name} -P", workdir)
    if result is None or result.returncode != 0:
        abort(4, "odgi build failed.", "Verify the GFA 1.0 file is valid.")

    if not os.path.exists(og_path):
        abort(4, f"Expected .og file not created: {og_path}")

    ok(f"ODGI graph built: {og_path}")
    save_checkpoint(STEP)
    return og_path

# ────────────────────────────────────────────────────────────────────────────
# STEP 5 — Sort the ODGI graph
# ────────────────────────────────────────────────────────────────────────────
def step5(og_path, workdir):
    STEP = "step5_sort"
    sorted_name = Path(og_path).stem + ".sorted.og"
    sorted_path = os.path.join(workdir, sorted_name)

    if already_done(STEP):
        ok("Step 5 — skipping (already done)")
        return sorted_path

    info("── STEP 5: Sorting ODGI graph ────────────────────────────────")
    og_name = Path(og_path).name

    result = docker_run("odgi", f"sort -O -i {og_name} -o {sorted_name} -P", workdir)
    if result is None or result.returncode != 0:
        abort(5, "odgi sort failed.")

    if not os.path.exists(sorted_path):
        abort(5, f"Sorted .og file not created: {sorted_path}")

    ok(f"Sorted graph: {sorted_path}")
    save_checkpoint(STEP)
    return sorted_path

# ────────────────────────────────────────────────────────────────────────────
# STEP 6 — Extract region from sorted graph
# ────────────────────────────────────────────────────────────────────────────
def step6(sorted_og_path, region, workdir):
    STEP = "step6_extract"
    sub_name = Path(sorted_og_path).stem + ".subgraph.og"
    sub_path = os.path.join(workdir, sub_name)

    if already_done(STEP):
        ok("Step 6 — skipping (already done)")
        return sub_path

    info(f"── STEP 6: Extracting region {region} ────────────────────────")
    info("This may take several minutes for large graphs.")
    sorted_name = Path(sorted_og_path).name

    result = docker_run("odgi",
                        f"extract -i {sorted_name} -r {region} -o - -P -c 0 -E",
                        workdir, capture_output=True)

    if result is None or result.returncode != 0:
        abort(6, "odgi extract failed.",
              "Check that the region reference name matches your GFA.")

    data = result.stdout if isinstance(result.stdout, bytes) else result.stdout.encode("latin-1")
    with open(sub_path, "wb") as f:
        f.write(data)

    if os.path.getsize(sub_path) == 0:
        abort(6, f"Subgraph is empty — region '{region}' may not exist in this graph.")

    ok(f"Subgraph extracted: {sub_path}")
    save_checkpoint(STEP)
    return sub_path

# ────────────────────────────────────────────────────────────────────────────
# STEP 7 — Convert subgraph .og → .gfa
# ────────────────────────────────────────────────────────────────────────────
def step7(sub_og_path, workdir):
    STEP = "step7_view"
    sub_gfa_name = Path(sub_og_path).name + ".gfa"
    sub_gfa_path = os.path.join(workdir, sub_gfa_name)

    if already_done(STEP):
        ok("Step 7 — skipping (already done)")
        return sub_gfa_path

    info("── STEP 7: Converting subgraph .og → .gfa ────────────────────")
    sub_name = Path(sub_og_path).name

    result = docker_run("odgi", f"view -i {sub_name} -g", workdir, capture_output=True)

    if result is None or result.returncode != 0:
        abort(7, "odgi view failed.")

    data = result.stdout if isinstance(result.stdout, bytes) else result.stdout.encode()
    with open(sub_gfa_path, "wb") as f:
        f.write(data)

    if os.path.getsize(sub_gfa_path) == 0:
        abort(7, "Subgraph GFA is empty after conversion.")

    ok(f"Subgraph GFA: {sub_gfa_path}")
    save_checkpoint(STEP)
    return sub_gfa_path

# ────────────────────────────────────────────────────────────────────────────
# STEP 8 — Extract min/max node IDs from subgraph
# ────────────────────────────────────────────────────────────────────────────
def step8(sub_gfa_path):
    STEP = "step8_nodes"
    RANGE_FILE = "node_range.txt"

    if already_done(STEP):
        ok("Step 8 — skipping (already done)")
        if os.path.exists(RANGE_FILE):
            with open(RANGE_FILE) as f:
                a, b = f.read().strip().split()
                return int(a), int(b)
        warn("Checkpoint found but node_range.txt missing — re-running step 8.")

    info("── STEP 8: Extracting node ID range from subgraph ────────────")

    nodes = []
    try:
        with open(sub_gfa_path, "rb") as f:
            for raw in f:
                line = raw.decode("utf-8", errors="ignore")
                parts = line.strip().split("\t")
                if parts and parts[0] == "S":
                    try:
                        nodes.append(int(parts[1]))
                    except (ValueError, IndexError):
                        pass
    except Exception as e:
        abort(8, f"Cannot read subgraph GFA: {e}")

    if not nodes:
        abort(8, "No segment lines (S) found in subgraph.")

    min_n, max_n = min(nodes), max(nodes)
    ok(f"Min node: {min_n}  |  Max node: {max_n}  |  Total nodes: {len(nodes)}")

    with open(RANGE_FILE, "w") as f:
        f.write(f"{min_n} {max_n}\n")

    save_checkpoint(STEP)
    return min_n, max_n

# ────────────────────────────────────────────────────────────────────────────
# STEP 9 — Filter VCF to node range
# ────────────────────────────────────────────────────────────────────────────
def step9(vcf_path, min_n, max_n, workdir):
    STEP = "step9_filter"
    subset_path = os.path.join(workdir, "subset.vcf")

    if already_done(STEP):
        ok("Step 9 — skipping (already done)")
        return subset_path

    info(f"── STEP 9: Filtering VCF (nodes {min_n}–{max_n}) ─────────────")

    headers = []
    kept = []
    skipped = 0

    open_fn = gzip.open if vcf_path.endswith(".gz") else open

    try:
        with open_fn(vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    headers.append(line)
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    skipped += 1
                    continue
                node_ids = re.split(r"[<>]", parts[2])
                matched = False
                for nid in node_ids:
                    if nid.strip().isdigit():
                        n = int(nid)
                        if min_n <= n <= max_n:
                            matched = True
                            break
                if matched:
                    kept.append(line)
                else:
                    skipped += 1
    except Exception as e:
        abort(9, f"Failed to read VCF: {e}", "Verify the VCF file is not corrupted.")

    with open(subset_path, "w") as f:
        f.writelines(headers)
        f.writelines(kept)

    ok(f"VCF filtered — kept: {len(kept)} variants, skipped: {skipped}")
    if len(kept) == 0:
        warn("No variants matched the node range.")

    save_checkpoint(STEP)
    return subset_path


# ────────────────────────────────────────────────────────────────────────────
# LUCAS STEP 10 HELPERS — VCF filtering + summary (Lucas Borges dos Santos)
# ────────────────────────────────────────────────────────────────────────────
def extract_structural_variants(vcf_filepath: str, min_length: int = 50) -> list:
    """
    Parses a VCF file and extracts variants where REF or any ALT exceeds min_length.
    Uses actual sample names from VCF header.
    """
    import os
    if not os.path.exists(vcf_filepath):
        raise FileNotFoundError(f"VCF file not found at: {vcf_filepath}")

    filtered_variants = []
    sample_names = []

    with open(vcf_filepath, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                if len(columns) > 9:
                    sample_names = columns[9:]
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom    = parts[0]
            pos      = parts[1]
            ref_seq  = parts[3]
            alt_seqs = parts[4].split(',')
            ref_len  = len(ref_seq)
            alt_lens = [len(alt) for alt in alt_seqs]
            if ref_len > min_length or any(a > min_length for a in alt_lens):
                genotypes = {}
                if len(parts) > 9:
                    for i, sample_data in enumerate(parts[9:]):
                        s_name = sample_names[i] if i < len(sample_names) else f"Sample_{i+1}"
                        gt_call = sample_data.split(':')[0]
                        genotypes[s_name] = gt_call
                filtered_variants.append({
                    "chromosome": chrom,
                    "position":   pos,
                    "ref_length": ref_len,
                    "alt_lengths": alt_lens,
                    "genotypes":  genotypes
                })
    return filtered_variants


def summarize_variants(variants_data: list) -> str:
    """
    Generates a human-readable summary of the extracted variants (Lucas style).
    """
    if not variants_data:
        return "No variants found exceeding the length threshold."
    lines = []
    lines.append("--- Structural Variants Summary ---")
    lines.append(f"Total large variants detected: {len(variants_data)}\n")
    for i, var in enumerate(variants_data, 1):
        lines.append(f"Variant {i} | Location: chr{var['chromosome']}:{var['position']}")
        lines.append(f"  - REF Length: {var['ref_length']} bp")
        alt_str = ", ".join([f"{l} bp" for l in var['alt_lengths']])
        lines.append(f"  - ALT Lengths: {alt_str}")
        if var['genotypes']:
            gt_str = ", ".join([f"{s}: {g}" for s, g in var['genotypes'].items()])
            lines.append(f"  - Genotypes: {gt_str}")
        else:
            lines.append("  - Genotypes: None available")
        lines.append("-" * 40)
    return "\n".join(lines)

# ────────────────────────────────────────────────────────────────────────────
# STEP 10 — Parse VCF → Human-readable report (variants > 50bp)
# ────────────────────────────────────────────────────────────────────────────
def step10(subset_vcf_path, workdir):
    STEP = "step10_parse"
    out_path = os.path.join(workdir, "variants_human_readable.txt")

    if already_done(STEP):
        ok("Step 10 — skipping (already done)")
        return out_path

    info("── STEP 10: Parsing variants > 50bp → human-readable report ──")

    variants = []
    skipped_small = 0
    malformed = 0

    try:
        with open(subset_vcf_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    malformed += 1
                    continue

                chrom  = parts[0]
                pos    = parts[1]
                vid    = parts[2]
                ref    = parts[3]
                alt    = parts[4]
                qual   = parts[5] if len(parts) > 5 else "."
                filt   = parts[6] if len(parts) > 6 else "."
                info_f = parts[7] if len(parts) > 7 else "."

                ref_len = len(ref)
                alt_len = len(alt) if alt not in (".", "*") else 0
                var_len = max(ref_len, alt_len)

                if var_len <= 50:
                    skipped_small += 1
                    continue

                if ref_len == alt_len:
                    vtype = "SNP" if ref_len == 1 else "MNP"
                elif alt_len > ref_len:
                    vtype = "INSERTION"
                elif ref_len > alt_len:
                    vtype = "DELETION"
                else:
                    vtype = "COMPLEX_SV"

                haps = re.findall(r"[><](\d+)", vid)

                info_dict = {}
                if info_f not in (".", ""):
                    for item in info_f.split(";"):
                        if "=" in item:
                            k, v = item.split("=", 1)
                            info_dict[k] = v
                        else:
                            info_dict[item] = True

                gt_data = {}
                if len(parts) > 8:
                    fmt_keys = parts[8].split(":")
                    for idx, sample in enumerate(parts[9:], 1):
                        vals = sample.split(":")
                        gt_data[f"Sample_{idx}"] = dict(zip(fmt_keys, vals))

                variants.append({
                    "chrom": chrom, "pos": pos, "id": vid,
                    "ref": ref, "alt": alt,
                    "qual": qual, "filter": filt,
                    "type": vtype,
                    "ref_len": ref_len, "alt_len": alt_len,
                    "var_len": var_len,
                    "haplotypes": haps,
                    "info": info_dict,
                    "genotypes": gt_data
                })

    except Exception as e:
        abort(10, f"Failed to parse VCF: {e}")

    # Generate Lucas-style summary first
    lucas_variants = extract_structural_variants(subset_vcf_path, min_length=50)
    lucas_summary  = summarize_variants(lucas_variants)

    with open(out_path, "w") as o:
        o.write(lucas_summary)
        o.write("\n\n")
        o.write("=" * 80 + "\n")
        o.write("  PANGENOME VARIANT REPORT — Human-Readable Output\n")
        o.write(f"  Generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        o.write(f"  Input VCF : {subset_vcf_path}\n")
        o.write("  Filter    : Variants > 50 bp only\n")
        o.write("=" * 80 + "\n\n")

        o.write("SUMMARY\n")
        o.write("-" * 40 + "\n")
        o.write(f"  Variants reported (>50bp)  : {len(variants)}\n")
        o.write(f"  Variants skipped (<=50bp)  : {skipped_small}\n")
        o.write(f"  Malformed lines skipped    : {malformed}\n\n")

        type_counts = {}
        for v in variants:
            type_counts[v["type"]] = type_counts.get(v["type"], 0) + 1
        if type_counts:
            o.write("  By variant type:\n")
            for vt, ct in sorted(type_counts.items()):
                o.write(f"    {vt:<20} {ct}\n")
        o.write("\n")

        if not variants:
            o.write("No variants > 50bp found in this region.\n")
        else:
            o.write("=" * 80 + "\n")
            o.write("  VARIANTS\n")
            o.write("=" * 80 + "\n\n")

            for i, v in enumerate(variants, 1):
                o.write(f"[{i}] {v['type']} — {v['var_len']} bp\n")
                o.write(f"    Chromosome  : {v['chrom']}\n")
                o.write(f"    Position    : {v['pos']}\n")
                o.write(f"    Variant ID  : {v['id']}\n")
                o.write(f"    REF length  : {v['ref_len']} bp\n")
                o.write(f"    ALT length  : {v['alt_len']} bp\n")
                o.write(f"    Quality     : {v['qual']}\n")
                o.write(f"    Filter      : {v['filter']}\n")

                ref_d = v["ref"] if len(v["ref"]) <= 60 else v["ref"][:57] + "..."
                alt_d = v["alt"] if len(v["alt"]) <= 60 else v["alt"][:57] + "..."
                o.write(f"    REF seq     : {ref_d}\n")
                o.write(f"    ALT seq     : {alt_d}\n")

                if v["haplotypes"]:
                    o.write(f"    Node path   : {' → '.join(v['haplotypes'])}\n")

                if v["info"]:
                    o.write("    INFO        :\n")
                    for k, val in v["info"].items():
                        o.write(f"      {k}: {val}\n")

                if v["genotypes"]:
                    o.write("    Genotypes   :\n")
                    for samp, gt in v["genotypes"].items():
                        o.write(f"      {samp}: GT={gt.get('GT','?')}\n")

                o.write("\n")

    ok(f"Report written: {out_path}")
    ok(f"Variants reported: {len(variants)}  |  Skipped (<=50bp): {skipped_small}")
    save_checkpoint(STEP)
    return out_path

# ────────────────────────────────────────────────────────────────────────────
# Final summary
# ────────────────────────────────────────────────────────────────────────────
def print_summary(output_path):
    done = load_checkpoints()
    steps = [
        ("step1_paths",    "Step 1  — Resolve paths"),
        ("step2_unzip",    "Step 2  — Decompress files"),
        ("step3_version",  "Step 3  — GFA version check/convert"),
        ("step4_build_og", "Step 4  — Build ODGI graph"),
        ("step5_sort",     "Step 5  — Sort graph"),
        ("step6_extract",  "Step 6  — Extract region"),
        ("step7_view",     "Step 7  — Subgraph → GFA"),
        ("step8_nodes",    "Step 8  — Get node range"),
        ("step9_filter",   "Step 9  — Filter VCF"),
        ("step10_parse",   "Step 10 — Parse → human-readable"),
    ]
    print()
    print(BOLD + "=" * 60 + RESET)
    print(BOLD + "  PIPELINE COMPLETE" + RESET)
    print(BOLD + "=" * 60 + RESET)
    for key, label in steps:
        mark = GREEN + "✅" + RESET if key in done else RED + "❌" + RESET
        print(f"  {mark}  {label}")
    print()
    print(BOLD + f"  Output → {output_path}" + RESET)
    print(BOLD + "=" * 60 + RESET)
    print()

# ────────────────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Pangenome Pipeline 4 — VCF Region Extractor + Human-Readable Parser"
    )
    parser.add_argument("--gfa",     required=True, help="Input .gfa or .gfa.gz")
    parser.add_argument("--vcf",     required=True, help="Input .vcf or .vcf.gz")
    parser.add_argument("--region",  required=False, default=None,
        help='Region string e.g. "MM26#0#Hg_chrom6_MM26:520000-570000" — auto-detected if not provided')
    parser.add_argument("--start",   type=int, default=None,
        help="Region start coordinate (used with --end when --region not given)")
    parser.add_argument("--end",     type=int, default=None,
        help="Region end coordinate (used with --start when --region not given)")
    parser.add_argument("--workdir", default="pipeline_output",
        help="Output directory for all intermediate + final files")
    parser.add_argument("--reset",   action="store_true",
        help="Clear checkpoints and run from scratch")
    parser.add_argument("--skip-step", type=int, dest="skip_step",
        help="Mark steps 1 through N-1 as done and start from step N")
    args = parser.parse_args()

    if args.region is None and (args.start is None or args.end is None):
        parser.error("Provide either --region or both --start and --end.")

    os.makedirs(args.workdir, exist_ok=True)

    print(BOLD + CYAN)
    print("=" * 60)
    print("  PANGENOME PIPELINE 4")
    print(f"  GFA    : {args.gfa}")
    print(f"  VCF    : {args.vcf}")
    if args.region:
        print(f"  Region : {args.region}")
    else:
        print(f"  Region : auto-detect (start={args.start}, end={args.end})")
    print(f"  Output : {args.workdir}")
    print("=" * 60)
    print(RESET, flush=True)

    if args.reset:
        clear_checkpoints()

    if args.skip_step:
        step_keys = [
            "step1_paths", "step2_unzip", "step3_version",
            "step4_build_og", "step5_sort", "step6_extract",
            "step7_view", "step8_nodes", "step9_filter"
        ]
        warn(f"Skipping to step {args.skip_step} — marking steps 1–{args.skip_step-1} complete.")
        for i in range(min(args.skip_step - 1, len(step_keys))):
            if step_keys[i] not in load_checkpoints():
                save_checkpoint(step_keys[i])

    check_docker()
    check_images()

    workdir = str(Path(args.workdir).resolve())

    gfa_abs, vcf_abs = step1(args.gfa, args.vcf)
    gfa_unc, vcf_path = step2(gfa_abs, vcf_abs, workdir)

    gfa_wdir = os.path.join(workdir, Path(gfa_unc).name)
    if not os.path.exists(gfa_wdir):
        shutil.copy2(gfa_unc, gfa_wdir)

    gfa_10 = step3(gfa_wdir, workdir)

    if args.region:
        region = args.region
    else:
        region = build_region(gfa_10, args.start, args.end)

    og_path   = step4(gfa_10, workdir)
    sorted_og = step5(og_path, workdir)
    sub_og    = step6(sorted_og, region, workdir)
    sub_gfa   = step7(sub_og, workdir)
    min_n, max_n = step8(sub_gfa)

    vcf_wdir = os.path.join(workdir, Path(vcf_path).name)
    if not os.path.exists(vcf_wdir):
        shutil.copy2(vcf_path, vcf_wdir)

    subset_vcf  = step9(vcf_wdir, min_n, max_n, workdir)
    output_path = step10(subset_vcf, workdir)

    print_summary(output_path)

if __name__ == "__main__":
    main()
