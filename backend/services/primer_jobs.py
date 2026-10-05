from __future__ import annotations

import datetime as dt
import json
import re
import shutil
import subprocess
import threading
import uuid
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple

from pydantic import BaseModel, Field, model_validator

from ..api_limits import get_api_limits
from ..config import JOB_DIR, ORTHOLOG_BIN, ROOT
from ..db import connect
from .job_files import write_json


PRIMER_LIMITS = get_api_limits().primer
PRIMER_CANDIDATE_POOL = 1000
PRIMER_VISIBLE_PAIR_LIMIT = 20
PRIMER_BOWTIE2_THREADS = 2
PRIMER_LOGIC_VERSION = 5
CLONING_PRIMER3_LOCK = threading.Lock()


def clean_dna(value: str, *, allow_empty: bool = False) -> str:
    parts = [line.strip() for line in str(value or "").splitlines() if not line.lstrip().startswith(">")]
    sequence = re.sub(r"\s+", "", "".join(parts)).upper()
    if not sequence and allow_empty:
        return ""
    if not sequence:
        raise ValueError("sequence is empty")
    if re.search(r"[^ACGT]", sequence):
        raise ValueError("sequence must contain only A, C, G, and T")
    return sequence


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class PrimerRequest(BaseModel):
    species_id: str = Field(min_length=1)
    source: Literal["gene", "sequence"] = "gene"
    template_type: Literal["gdna", "cdna", "cdna_cloning"] = "gdna"
    gene_id: Optional[str] = None
    transcript_id: Optional[str] = None
    sequence: Optional[str] = None
    seqid: Optional[str] = None
    start: Optional[int] = Field(default=None, ge=1)
    end: Optional[int] = Field(default=None, ge=1)
    flank: int = Field(default=0, ge=0, le=PRIMER_LIMITS.max_flank)
    genomic_region: Literal["gene_body", "tss_window"] = "gene_body"
    tss_upstream: int = Field(default=1000, ge=0, le=PRIMER_LIMITS.max_flank)
    tss_downstream: int = Field(default=300, ge=0, le=PRIMER_LIMITS.max_flank)
    junction_mode: Literal["no_junction", "allow_junction"] = "no_junction"
    stop_codon_mode: Literal["keep", "remove"] = "keep"
    forward_tail: str = Field(default="", max_length=200)
    reverse_tail: str = Field(default="", max_length=200)
    product_size_range: str = Field(default="100-300", pattern=r"^\d+-\d+(,\d+-\d+)*$")
    primer_min_size: int = Field(default=18, ge=15, le=40)
    primer_max_size: int = Field(default=24, ge=15, le=40)
    primer_min_tm: float = Field(default=57.0, ge=40.0, le=90.0)
    primer_max_tm: float = Field(default=63.0, ge=40.0, le=90.0)
    primer_min_gc: float = Field(default=55.0, ge=0.0, le=100.0)
    primer_max_gc: float = Field(default=65.0, ge=0.0, le=100.0)
    check_specificity: bool = True
    max_region_length: int = Field(default=PRIMER_LIMITS.max_region_length, ge=100, le=PRIMER_LIMITS.max_region_length)
    max_offtarget_product_size: int = Field(default=1000, ge=50, le=PRIMER_LIMITS.max_offtarget_product_size)

    @model_validator(mode="before")
    @classmethod
    def cloning_defaults(cls, values: Any) -> Any:
        if isinstance(values, dict) and values.get("template_type") == "cdna_cloning":
            return {"primer_min_size": 20, "primer_max_size": 40,
                    "primer_min_tm": 57.0, "primer_max_tm": 65.0,
                    "primer_min_gc": 30.0, "primer_max_gc": 70.0, **values}
        return values

    @model_validator(mode="after")
    def validate_target(self) -> "PrimerRequest":
        has_gene = bool(self.gene_id)
        has_region = self.seqid is not None or self.start is not None or self.end is not None
        has_sequence = bool(self.sequence and self.sequence.strip())
        if self.source == "sequence":
            if has_gene or has_region:
                raise ValueError("sequence source cannot include a gene or genomic region")
            if not has_sequence:
                raise ValueError("sequence is required for sequence source")
        elif has_sequence:
            raise ValueError("gene source cannot include a submitted sequence")
        if self.source == "gene" and has_gene and has_region:
            raise ValueError("provide either gene_id or seqid/start/end, not both")
        if self.source == "gene" and not has_gene and not (self.seqid and self.start and self.end):
            raise ValueError("provide gene_id or seqid/start/end")
        if self.start and self.end and self.start > self.end:
            raise ValueError("start must be <= end")
        if self.primer_min_size > self.primer_max_size:
            raise ValueError("primer_min_size must be <= primer_max_size")
        if self.template_type != "cdna_cloning" and self.primer_max_size > 35:
            raise ValueError("Non-cloning primer length must be <= 35 nt")
        if self.primer_min_tm > self.primer_max_tm:
            raise ValueError("primer_min_tm must be <= primer_max_tm")
        if self.primer_min_gc > self.primer_max_gc:
            raise ValueError("primer_min_gc must be <= primer_max_gc")
        if self.genomic_region == "tss_window" and (self.source != "gene" or not has_gene or self.template_type != "gdna"):
            raise ValueError("tss_window requires gene_id")
        if self.template_type != "gdna" and self.source != "gene":
            raise ValueError("cDNA template modes require a gene source")
        return self


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def job_root() -> Path:
    path = JOB_DIR / "primer"
    path.mkdir(parents=True, exist_ok=True)
    return path


def job_path(job_id: str) -> Path:
    if not re.fullmatch(r"[a-f0-9-]{36}", job_id):
        raise ValueError("invalid job_id")
    return job_root() / job_id


def read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def create_primer_job(request: PrimerRequest) -> Dict[str, Any]:
    target, _ = resolve_template(request)
    if request.check_specificity and request.template_type != "cdna_cloning":
        require_bowtie2(request.species_id)

    job_id = str(uuid.uuid4())
    path = job_path(job_id)
    path.mkdir(parents=True, exist_ok=False)
    created_at = utc_now()
    write_json(path / "request.json", request.model_dump())
    write_json(path / "target.json", target)
    write_json(
        path / "status.json",
        {
            "job_id": job_id,
            "status": "queued",
            "created_at": created_at,
            "updated_at": created_at,
            "species_id": request.species_id,
            "target": target,
            "result_count": 0,
        },
    )
    return {"job_id": job_id, "status": "queued", "created_at": created_at, "species_id": request.species_id, "target": target}


def execute_primer_job(job_id: str) -> None:
    path = job_path(job_id)
    request = PrimerRequest(**read_json(path / "request.json"))
    target = read_json(path / "target.json")
    status_path = path / "status.json"
    started_at = utc_now()
    write_json(
        status_path,
        {
            "job_id": job_id,
            "status": "running",
            "created_at": read_json(status_path)["created_at"],
            "started_at": started_at,
            "updated_at": started_at,
            "species_id": request.species_id,
            "target": target,
            "result_count": 0,
        },
    )
    try:
        target, sequence = resolve_template(request)
        diagnostics: Dict[str, Any] = {}
        if request.template_type == "cdna_cloning":
            candidates = full_length_cloning_pairs(sequence, request, diagnostics)
        else:
            primer3_output = run_primer3(path, request, target, sequence)
            candidates = parse_primer3(primer3_output)
        generated_count = len(candidates)
        if request.template_type == "cdna":
            candidates = filter_cdna_candidates(candidates, target, request.junction_mode)
        if request.template_type == "cdna":
            annotate_cdna_junctions(candidates, target)
        if request.template_type == "cdna_cloning":
            apply_cloning_tails(candidates, request.forward_tail, request.reverse_tail)
        if request.check_specificity and request.template_type != "cdna_cloning" and candidates:
            add_specificity(path, request, target, candidates)
        eligible_count = len(candidates)
        candidates = candidates[:3] if request.template_type == "cdna_cloning" else select_visible_candidates(
            candidates,
            enforce_spatial_diversity=request.template_type in {"gdna", "cdna"},
        )
        finished_at = utc_now()
        write_json(path / "results.json", {
            "target": target, "primers": candidates, "logic_version": PRIMER_LOGIC_VERSION,
            "design_engine": "fixed_end_enumeration+primer3-py" if request.template_type == "cdna_cloning" else "primer3_core",
            "diagnostics": diagnostics,
            "candidate_counts": {"generated": generated_count, "eligible": eligible_count, "displayed": len(candidates)},
            "empty_reason": (
                diagnostics.get("empty_reason", "No full-length CDS primer pairs meet the selected constraints.")
                if request.template_type == "cdna_cloning" else
                "No multi-exon amplicons meet the selected constraints. Single-exon products are excluded; a single-exon transcript requires Exon only mode."
                if request.template_type == "cdna" and request.junction_mode == "allow_junction" else
                "No primer pairs meet the selected size, Tm, GC and primer-placement constraints."
            ) if not candidates else None,
        })
        write_json(
            status_path,
            {
                "job_id": job_id,
                "status": "completed",
                "created_at": read_json(status_path)["created_at"],
                "started_at": started_at,
                "finished_at": finished_at,
                "updated_at": finished_at,
                "species_id": request.species_id,
                "target": target,
                "result_count": len(candidates),
            },
        )
    except Exception as exc:
        failed_at = utc_now()
        write_json(
            status_path,
            {
                "job_id": job_id,
                "status": "failed",
                "created_at": read_json(status_path)["created_at"],
                "started_at": started_at,
                "failed_at": failed_at,
                "updated_at": failed_at,
                "species_id": request.species_id,
                "target": target,
                "result_count": 0,
                "error": str(exc),
            },
        )


def get_primer_job(job_id: str) -> Dict[str, Any]:
    path = job_path(job_id)
    status_path = path / "status.json"
    if not status_path.exists():
        raise FileNotFoundError(job_id)
    status = read_json(status_path)
    request = PrimerRequest(**read_json(path / "request.json"))
    status["design"] = request.model_dump(exclude={"sequence"})
    status["design"]["check_specificity"] = request.check_specificity and request.template_type != "cdna_cloning"
    result_path = path / "results.json"
    status["results"] = read_json(result_path) if status["status"] == "completed" and result_path.exists() else None
    return status


def resolve_target(request: PrimerRequest) -> Dict[str, Any]:
    if request.gene_id:
        with connect() as conn:
            row = conn.execute(
                "select gene_id, seqid, start, end, strand from genes where species_id = ? and gene_id = ?",
                (request.species_id, request.gene_id),
            ).fetchone()
        if row is None:
            raise ValueError(f"gene not found: {request.species_id}/{request.gene_id}")
        seqid, start, end = row["seqid"], int(row["start"]), int(row["end"])
        strand = row["strand"]
        target_id = row["gene_id"]
    else:
        seqid, start, end = str(request.seqid), int(request.start or 0), int(request.end or 0)
        strand = "+"
        target_id = f"{seqid}:{start}-{end}"

    with connect() as conn:
        chrom = conn.execute(
            "select length from chromosomes where species_id = ? and seqid = ?",
            (request.species_id, seqid),
        ).fetchone()
    if chrom is None:
        raise ValueError(f"seqid not found: {request.species_id}/{seqid}")
    chrom_length = int(chrom["length"])
    if request.genomic_region == "tss_window" and request.gene_id:
        # The transcription start is the low coordinate on + and high coordinate on -.
        tss = start if strand != "-" else end
        if strand == "-":
            region_start = max(1, tss - request.tss_downstream)
            region_end = min(chrom_length, tss + request.tss_upstream)
        else:
            region_start = max(1, tss - request.tss_upstream)
            region_end = min(chrom_length, tss + request.tss_downstream)
    else:
        region_start = max(1, start - request.flank)
        region_end = min(chrom_length, end + request.flank)
    if region_end - region_start + 1 > request.max_region_length:
        raise ValueError(f"target region is too long; max {request.max_region_length} bp")
    return {
        "target_id": target_id,
        "seqid": seqid,
        "start": region_start,
        "end": region_end,
        "core_start": start,
        "core_end": end,
        "strand": strand,
        "genomic_region": request.genomic_region,
        "length": region_end - region_start + 1,
    }


def representative_transcript(rows: List[Any]) -> Any:
    for row in rows:
        transcript_id = str(row["transcript_id"])
        if transcript_id.endswith(".1") or transcript_id.endswith("-01"):
            return row
    return rows[0]


def transcript_template(request: PrimerRequest) -> Tuple[Dict[str, Any], str]:
    if not request.gene_id:
        raise ValueError("gene_id is required for transcript templates")
    with connect() as conn:
        transcript_rows = conn.execute(
            "select transcript_id, gene_id, seqid, start, end, strand from transcripts where species_id = ? and gene_id = ? order by transcript_id",
            (request.species_id, request.gene_id),
        ).fetchall()
        if not transcript_rows:
            raise ValueError(f"no transcripts found: {request.species_id}/{request.gene_id}")
        if request.transcript_id:
            transcript = next((row for row in transcript_rows if row["transcript_id"] == request.transcript_id), None)
            if transcript is None:
                raise ValueError("transcript does not belong to the selected gene")
        else:
            transcript = representative_transcript(transcript_rows)
        feature_type = "CDS" if request.template_type == "cdna_cloning" else "exon"
        features = conn.execute(
            "select start, end, strand from features where species_id = ? and transcript_id = ? and feature_type = ? order by start, end",
            (request.species_id, transcript["transcript_id"], feature_type),
        ).fetchall()
    if not features:
        raise ValueError(f"no {feature_type} features found for transcript {transcript['transcript_id']}")
    strand = str(transcript["strand"] or "+")
    ordered_features = list(reversed(features)) if strand == "-" else list(features)
    fragments: List[str] = []
    junctions: List[int] = []
    exon_segments: List[Dict[str, int]] = []
    offset = 0
    for index, feature in enumerate(ordered_features):
        fragment = fetch_species_sequence(request.species_id, str(transcript["seqid"]), int(feature["start"]), int(feature["end"]))
        if strand == "-":
            fragment = reverse_complement(fragment)
        fragments.append(fragment)
        exon_segments.append({
            "number": index + 1, "start": offset + 1, "end": offset + len(fragment),
            "genomic_start": int(feature["start"]), "genomic_end": int(feature["end"]),
        })
        offset += len(fragment)
        if index < len(ordered_features) - 1:
            junctions.append(offset)
    sequence = "".join(fragments)
    if request.template_type == "cdna_cloning" and request.stop_codon_mode == "remove":
        if len(sequence) < 6:
            raise ValueError("CDS template is too short to remove a stop codon")
        sequence = sequence[:-3]
        junctions = [junction for junction in junctions if junction < len(sequence)]
    template_label = "CDS" if request.template_type == "cdna_cloning" else "spliced cDNA"
    target = {
        "target_id": request.gene_id,
        "transcript_id": str(transcript["transcript_id"]),
        "seqid": str(transcript["seqid"]),
        "start": int(transcript["start"]),
        "end": int(transcript["end"]),
        "core_start": int(transcript["start"]),
        "core_end": int(transcript["end"]),
        "strand": strand,
        "genomic_region": request.template_type,
        "template_label": template_label,
        "junctions": junctions,
        "exon_segments": exon_segments if feature_type == "exon" else [],
        "length": len(sequence),
    }
    return target, sequence


def resolve_template(request: PrimerRequest) -> Tuple[Dict[str, Any], str]:
    if request.source == "sequence":
        sequence = clean_dna(request.sequence or "")
        return {
            "target_id": "submitted_sequence",
            "seqid": None,
            "start": 1,
            "end": len(sequence),
            "core_start": 1,
            "core_end": len(sequence),
            "strand": "+",
            "genomic_region": "submitted_sequence",
            "template_label": "submitted sequence",
            "junctions": [],
            "length": len(sequence),
        }, sequence
    if request.template_type in {"cdna", "cdna_cloning"}:
        return transcript_template(request)
    target = resolve_target(request)
    sequence = fetch_species_sequence(request.species_id, target["seqid"], target["start"], target["end"])
    target["template_label"] = "genomic DNA"
    target["junctions"] = []
    return target, sequence


def genome_paths(species_id: str) -> Tuple[Path, Path]:
    with connect() as conn:
        rows = conn.execute(
            "select asset_type, path from genome_assets where species_id = ? and asset_type in ('fasta', 'fai')",
            (species_id,),
        ).fetchall()
    values = {row["asset_type"]: Path(row["path"]) for row in rows}
    if "fasta" not in values:
        raise ValueError(f"genome FASTA asset not found: {species_id}")
    fasta = values["fasta"]
    fai = values.get("fai", Path(str(fasta) + ".fai"))
    if not fasta.exists():
        raise ValueError(f"genome FASTA does not exist: {fasta}")
    if not fai.exists():
        raise ValueError(f"genome FASTA index does not exist: {fai}")
    return fasta, fai


def fetch_species_sequence(species_id: str, seqid: str, start: int, end: int) -> str:
    fasta, fai = genome_paths(species_id)
    fai_row = None
    with fai.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == seqid:
                fai_row = fields
                break
    if fai_row is None:
        raise ValueError(f"seqid not found in FASTA index: {seqid}")
    length, offset, line_bases, line_width = map(int, fai_row[1:5])
    if end > length:
        raise ValueError(f"requested interval exceeds chromosome length: {seqid}:{start}-{end}")

    parts: List[str] = []
    with fasta.open("rb") as handle:
        pos = start
        while pos <= end:
            line_index = (pos - 1) // line_bases
            in_line = (pos - 1) % line_bases
            take = min(end - pos + 1, line_bases - in_line)
            byte_offset = offset + line_index * line_width + in_line
            handle.seek(byte_offset)
            parts.append(handle.read(take).decode("ascii"))
            pos += take
    return "".join(parts).upper()


def run_primer3(path: Path, request: PrimerRequest, target: Dict[str, Any], sequence: str) -> Dict[str, str]:
    optimal_size = min(max(20, request.primer_min_size), request.primer_max_size)
    optimal_tm = min(max(60.0, request.primer_min_tm), request.primer_max_tm)
    optimal_gc = min(max(50.0, request.primer_min_gc), request.primer_max_gc)
    product_size_range = f"{len(sequence)}-{len(sequence)}" if request.template_type == "cdna_cloning" else request.product_size_range.replace(",", " ")
    primer3_input = "\n".join(
        [
            f"SEQUENCE_ID={request.species_id}:{target['target_id']}",
            f"SEQUENCE_TEMPLATE={sequence}",
            "PRIMER_TASK=generic",
            "PRIMER_PICK_LEFT_PRIMER=1",
            "PRIMER_PICK_RIGHT_PRIMER=1",
            f"PRIMER_OPT_SIZE={optimal_size}",
            f"PRIMER_MIN_SIZE={request.primer_min_size}",
            f"PRIMER_MAX_SIZE={request.primer_max_size}",
            f"PRIMER_OPT_TM={optimal_tm:.1f}",
            f"PRIMER_MIN_TM={request.primer_min_tm:.1f}",
            f"PRIMER_MAX_TM={request.primer_max_tm:.1f}",
            f"PRIMER_OPT_GC_PERCENT={optimal_gc:.1f}",
            f"PRIMER_MIN_GC={request.primer_min_gc:.1f}",
            f"PRIMER_MAX_GC={request.primer_max_gc:.1f}",
            "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0",
            "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0",
            f"PRIMER_PRODUCT_SIZE_RANGE={product_size_range}",
            f"PRIMER_NUM_RETURN={PRIMER_CANDIDATE_POOL}",
            "PRIMER_EXPLAIN_FLAG=1",
            "=",
            "",
        ]
    )
    input_path = path / "primer3.input"
    output_path = path / "primer3.output"
    input_path.write_text(primer3_input, encoding="utf-8")
    cmd = [str(ORTHOLOG_BIN / "primer3_core"), str(input_path)]
    with output_path.open("w", encoding="utf-8") as out:
        subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.STDOUT, timeout=PRIMER_LIMITS.primer3_timeout_seconds)
    result: Dict[str, str] = {}
    for line in output_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line and line != "=":
            key, value = line.split("=", 1)
            result[key] = value
    if "PRIMER_ERROR" in result and result["PRIMER_ERROR"]:
        raise ValueError(result["PRIMER_ERROR"])
    return result


def parse_primer3(result: Dict[str, str]) -> List[Dict[str, Any]]:
    returned = int(result.get("PRIMER_PAIR_NUM_RETURNED", "0") or 0)
    primers: List[Dict[str, Any]] = []
    for idx in range(returned):
        left_pos, left_len = parse_position(result[f"PRIMER_LEFT_{idx}"])
        right_pos, right_len = parse_position(result[f"PRIMER_RIGHT_{idx}"])
        primers.append(
            {
                "rank": idx + 1,
                "left_sequence": result[f"PRIMER_LEFT_{idx}_SEQUENCE"],
                "right_sequence": result[f"PRIMER_RIGHT_{idx}_SEQUENCE"],
                "left_position": {"offset": left_pos, "length": left_len},
                "right_position": {"offset": right_pos, "length": right_len},
                "product_size": int(result[f"PRIMER_PAIR_{idx}_PRODUCT_SIZE"]),
                "penalty": float(result.get(f"PRIMER_PAIR_{idx}_PENALTY", 0.0)),
                "left_tm": float(result.get(f"PRIMER_LEFT_{idx}_TM", 0.0)),
                "right_tm": float(result.get(f"PRIMER_RIGHT_{idx}_TM", 0.0)),
                "left_gc_percent": float(result.get(f"PRIMER_LEFT_{idx}_GC_PERCENT", 0.0)),
                "right_gc_percent": float(result.get(f"PRIMER_RIGHT_{idx}_GC_PERCENT", 0.0)),
            }
        )
    return primers


def parse_position(value: str) -> Tuple[int, int]:
    offset, length = value.split(",", 1)
    return int(offset), int(length)


def full_length_cloning_pairs(sequence: str, request: PrimerRequest, diagnostics: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    # primer3-py exposes shared thermodynamic state; cloning jobs must not race it.
    with CLONING_PRIMER3_LOCK:
        return _full_length_cloning_pairs(sequence, request, diagnostics)


def _full_length_cloning_pairs(sequence: str, request: PrimerRequest, diagnostics: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """Enumerate fixed ends beyond Primer3 design's 36 nt limit; use its thermodynamics."""
    try:
        import primer3
    except ImportError as exc:
        raise ValueError("Full-length cloning requires primer3-py with thermodynamic parameters.") from exc
    diagnostics = diagnostics if diagnostics is not None else {}
    if len(sequence) < request.primer_min_size * 2:
        diagnostics["empty_reason"] = "CDS is too short for two primers at the selected minimum length."
        return []
    forward_tail = clean_dna(request.forward_tail, allow_empty=True)
    reverse_tail = reverse_complement(clean_dna(request.reverse_tail, allow_empty=True))
    if max(len(forward_tail), len(reverse_tail)) + request.primer_max_size > 60:
        raise ValueError("Cloning structure screening requires addition length + maximum primer length <= 60 nt.")
    ends = []
    reasons = []
    for name, tail in [("Forward", forward_tail), ("Reverse", reverse_tail)]:
        accepted, rejected = [], {}
        lengths = range(request.primer_min_size, min(request.primer_max_size, len(sequence)) + 1)
        for length in lengths:
            core = sequence[:length] if name == "Forward" else reverse_complement(sequence[-length:])
            gc = 100.0 * sum(base in "GC" for base in core) / length
            tm = primer3.calc_tm(core)
            failures = []
            if gc < request.primer_min_gc:
                failures.append("GC below minimum")
            if gc > request.primer_max_gc:
                failures.append("GC above maximum")
            if tm < request.primer_min_tm:
                failures.append("Tm below minimum")
            if tm > request.primer_max_tm:
                failures.append("Tm above maximum")
            if re.search(r"(A{6}|C{6}|G{6}|T{6})", core):
                failures.append("homopolymer longer than 5 nt")
            oligo = tail + core
            hairpin = primer3.calc_hairpin(oligo).tm
            homodimer = primer3.calc_homodimer(oligo).tm
            if hairpin > 47.0:
                failures.append("hairpin Tm above 47 C")
            if homodimer > 47.0:
                failures.append("homodimer Tm above 47 C")
            for reason in failures:
                rejected[reason] = rejected.get(reason, 0) + 1
            if not failures:
                accepted.append({"core": core, "oligo": oligo, "tm": tm, "gc": gc,
                                 "length": length, "hairpin": hairpin, "homodimer": homodimer})
        diagnostics[name.lower()] = {"considered": len(lengths), "accepted": len(accepted), "rejected": rejected}
        if not accepted:
            reasons.append(f"{name}: no eligible primer ({'; '.join(f'{key}: {count}' for key, count in rejected.items())}).")
        ends.append(accepted)
    pairs = []
    rejected_pairs = {"overlapping annealing regions": 0, "heterodimer Tm above 47 C": 0}
    optimum_tm = min(max(60.0, request.primer_min_tm), request.primer_max_tm)
    for left in ends[0]:
        for right in ends[1]:
            if left["length"] + right["length"] > len(sequence):
                rejected_pairs["overlapping annealing regions"] += 1
                continue
            heterodimer = primer3.calc_heterodimer(left["oligo"], right["oligo"]).tm
            if heterodimer > 47.0:
                rejected_pairs["heterodimer Tm above 47 C"] += 1
                continue
            delta = abs(left["tm"] - right["tm"])
            pairs.append({
                "left_sequence": left["core"], "right_sequence": right["core"],
                "left_position": {"offset": 0, "length": left["length"]},
                "right_position": {"offset": len(sequence) - 1, "length": right["length"]},
                "left_tm": left["tm"], "right_tm": right["tm"],
                "left_gc_percent": left["gc"], "right_gc_percent": right["gc"],
                "product_size": len(sequence), "full_length": True, "coverage": 100.0,
                "final_product_size": len(sequence) + len(forward_tail) + len(reverse_tail),
                "penalty": abs(left["tm"] - optimum_tm) + abs(right["tm"] - optimum_tm),
                "tm_difference": delta,
                "warnings": ["Primer Tm difference exceeds 5 C; review annealing conditions."] if delta > 5.0 else [],
                "structure": {"left_hairpin_tm": left["hairpin"], "right_hairpin_tm": right["hairpin"],
                              "left_homodimer_tm": left["homodimer"], "right_homodimer_tm": right["homodimer"],
                              "heterodimer_tm": heterodimer},
            })
    pairs.sort(key=lambda pair: (pair["tm_difference"] > 3.0, pair["tm_difference"], pair["penalty"]))
    for rank, pair in enumerate(pairs, 1):
        pair["rank"] = rank
    diagnostics["pair_rejections"] = rejected_pairs
    diagnostics["accepted_pairs"] = len(pairs)
    if not pairs:
        if not reasons:
            reasons = [f"Pair: {reason}: {count}." for reason, count in rejected_pairs.items() if count]
        diagnostics["empty_reason"] = "No full-length CDS primer pairs. " + " ".join(reasons)
    return pairs


def primer_spans_junction(offset: int, length: int, junctions: List[int]) -> bool:
    end_offset = offset + length - 1
    return any(offset < junction <= end_offset for junction in junctions)


def candidate_spans_junction(candidate: Dict[str, Any], junctions: List[int]) -> bool:
    left = candidate["left_position"]
    right = candidate["right_position"]
    left_spans = primer_spans_junction(int(left["offset"]), int(left["length"]), junctions)
    right_start = int(right["offset"]) - int(right["length"]) + 1
    right_spans = primer_spans_junction(right_start, int(right["length"]), junctions)
    return left_spans or right_spans


def candidate_product_spans_junction(candidate: Dict[str, Any], junctions: List[int]) -> bool:
    left_start = int(candidate["left_position"]["offset"])
    right_end = int(candidate["right_position"]["offset"])
    return any(left_start < junction <= right_end for junction in junctions)


def filter_cdna_candidates(candidates: List[Dict[str, Any]], target: Dict[str, Any], junction_mode: str) -> List[Dict[str, Any]]:
    """Separate single-exon and multi-exon products without primer overlap constraints."""
    junctions = [int(value) for value in target.get("junctions", [])]
    selected: List[Dict[str, Any]] = []
    for candidate in candidates:
        spans = candidate_product_spans_junction(candidate, junctions)
        if junction_mode == "no_junction" and spans:
            continue
        if junction_mode == "allow_junction" and not spans:
            continue
        selected.append(candidate)
    return selected


def exon_context(start_offset: int, end_offset: int, exon_segments: List[Dict[str, int]]) -> str:
    start = int(start_offset) + 1
    end = int(end_offset) + 1
    overlapping = [
        int(segment["number"])
        for segment in exon_segments
        if start <= int(segment["end"]) and end >= int(segment["start"])
    ]
    if not overlapping:
        return "inter-exon"
    if len(overlapping) == 1:
        return f"Exon {overlapping[0]}"
    return f"Exon {overlapping[0]}-{overlapping[-1]}"


def annotate_cdna_junctions(candidates: List[Dict[str, Any]], target: Dict[str, Any]) -> None:
    junctions = [int(value) for value in target.get("junctions", [])]
    exon_segments = target.get("exon_segments", [])
    for candidate in candidates:
        left = candidate["left_position"]
        right = candidate["right_position"]
        left_start = int(left["offset"])
        left_end = left_start + int(left["length"]) - 1
        right_end = int(right["offset"])
        right_start = right_end - int(right["length"]) + 1
        left_spans = primer_spans_junction(left_start, int(left["length"]), junctions)
        right_spans = primer_spans_junction(right_start, int(right["length"]), junctions)
        product_spans = any(left_start < junction <= right_end for junction in junctions)
        expected = expected_genome_pair(target, candidate)
        candidate["expected_genomic_product_size"] = expected["span"] if expected else None
        candidate["junction"] = {
            "left_context": exon_context(left_start, left_end, exon_segments),
            "right_context": exon_context(right_start, right_end, exon_segments),
            "left_spans_junction": left_spans,
            "right_spans_junction": right_spans,
            "pair_spans_junction": product_spans,
            "status": "junction-spanning primer" if left_spans or right_spans else "exon-only primer",
        }


def apply_cloning_tails(candidates: List[Dict[str, Any]], forward_tail: str, reverse_sense_tail: str) -> None:
    forward = clean_dna(forward_tail, allow_empty=True)
    reverse_sense = clean_dna(reverse_sense_tail, allow_empty=True)
    reverse_oligo_tail = reverse_complement(reverse_sense) if reverse_sense else ""
    for candidate in candidates:
        candidate["left_core_sequence"] = candidate["left_sequence"]
        candidate["right_core_sequence"] = candidate["right_sequence"]
        candidate["left_tail_sequence"] = forward
        candidate["right_tail_sequence"] = reverse_oligo_tail
        candidate["right_sense_tail_sequence"] = reverse_sense
        candidate["left_sequence"] = forward + candidate["left_core_sequence"]
        candidate["right_sequence"] = reverse_oligo_tail + candidate["right_core_sequence"]


def bowtie2_index_prefix(species_id: str) -> Optional[Path]:
    if not re.fullmatch(r"[A-Za-z0-9_-]+", species_id):
        return None
    with connect() as conn:
        has_assets = conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='search_assets'").fetchone()
        asset = conn.execute("SELECT index_prefix, status FROM search_assets WHERE species_id=? AND asset_type='bowtie2_genome'", (species_id,)).fetchone() if has_assets else None
    if asset is not None:
        if asset[1] != 'ready':
            return None
        prefix = Path(asset[0])
        if not prefix.is_absolute():
            prefix = ROOT / prefix
        from .bowtie_indexes import index_format
        try:
            index_format(prefix)
        except ValueError:
            return None
        return prefix
    # Compatibility for existing Cq/Da builds not yet registered in SQLite.
    prefix = ROOT / "data" / "species_assets" / species_id / "search" / "bowtie2" / f"{species_id}.genome"
    for extension in ("bt2", "bt2l"):
        suffixes = ("1", "2", "3", "4", "rev.1", "rev.2")
        if all(Path(f"{prefix}.{suffix}.{extension}").exists() for suffix in suffixes):
            return prefix
    return None


def bowtie2_binary() -> Optional[str]:
    configured = ORTHOLOG_BIN / "bowtie2"
    if configured.exists():
        return str(configured)
    return shutil.which("bowtie2")


def specificity_class(count: int) -> str:
    if count <= 0:
        return "none"
    if count == 1:
        return "unique"
    return "multi"


def write_fasta_record(handle: Any, identifier: str, sequence: str) -> None:
    handle.write(f">{identifier}\n{sequence}\n")


def require_bowtie2(species_id: str) -> Tuple[str, Path]:
    binary, prefix = bowtie2_binary(), bowtie2_index_prefix(species_id)
    if not binary:
        raise ValueError("Genome specificity requires Bowtie2; the executable is unavailable.")
    if not prefix:
        raise ValueError(f"Genome specificity requires a complete Bowtie2 index for {species_id}.")
    return binary, prefix


def expected_genome_binding(target: Dict[str, Any], offset: int, length: int, reverse: bool) -> Optional[Dict[str, Any]]:
    if not target.get("seqid"):
        return None
    if target.get("genomic_region") == "cdna":
        segment = next((s for s in target.get("exon_segments", [])
                        if int(s["start"]) - 1 <= offset and offset + length <= int(s["end"])), None)
        if segment is None:
            return None
        local_offset = offset - (int(segment["start"]) - 1)
        if target["strand"] == "-":
            end = int(segment["genomic_end"]) - local_offset
            start, strand = end - length + 1, "+" if reverse else "-"
        else:
            start = int(segment["genomic_start"]) + local_offset
            end, strand = start + length - 1, "-" if reverse else "+"
    else:
        start = int(target["start"]) + offset
        end, strand = start + length - 1, "-" if reverse else "+"
    return {"seqid": target["seqid"], "start": start, "end": end, "strand": strand}


def expected_genome_pair(target: Dict[str, Any], candidate: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    left, right = candidate["left_position"], candidate["right_position"]
    forward = expected_genome_binding(target, int(left["offset"]), int(left["length"]), False)
    reverse = expected_genome_binding(target, int(right["offset"]) - int(right["length"]) + 1, int(right["length"]), True)
    if forward is None or reverse is None:
        return None
    return {"forward": forward, "reverse": reverse,
            "span": max(forward["end"], reverse["end"]) - min(forward["start"], reverse["start"]) + 1}


def exact_sam_records(path: Path) -> Dict[str, List[Dict[str, Any]]]:
    records: Dict[str, List[Dict[str, Any]]] = {}
    # Missing/malformed output is a failed screen, never evidence for zero hits.
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError("Incomplete Bowtie2 SAM record.")
            flag = int(fields[1])
            if flag & 4:
                continue
            length = len(fields[9])
            if fields[5] not in (f"{length}M", f"{length}=") or "NM:i:0" not in fields[11:]:
                raise ValueError("Bowtie2 returned a non-exact or incomplete primer alignment.")
            record = {"seqid": fields[2], "start": int(fields[3]), "end": int(fields[3]) + length - 1,
                      "strand": "-" if flag & 16 else "+", "flag": flag,
                      "mate_seqid": fields[2] if fields[6] == "=" else fields[6], "mate_start": int(fields[7])}
            records.setdefault(fields[0], []).append(record)
    return records


def binding_key(record: Dict[str, Any]) -> Tuple[Any, ...]:
    return tuple(record[key] for key in ("seqid", "start", "end", "strand"))


def sam_pair_products(records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    products: Dict[Tuple[Any, ...], Dict[str, Any]] = {}
    for forward in records:
        if not forward["flag"] & 64 or not forward["flag"] & 2:
            continue
        for reverse in records:
            if not reverse["flag"] & 128 or not reverse["flag"] & 2:
                continue
            if (forward["mate_seqid"], forward["mate_start"]) != (reverse["seqid"], reverse["start"]):
                continue
            if (reverse["mate_seqid"], reverse["mate_start"]) != (forward["seqid"], forward["start"]):
                continue
            if forward["seqid"] != reverse["seqid"] or forward["strand"] == reverse["strand"]:
                continue
            key = binding_key(forward) + binding_key(reverse)
            products[key] = {"forward": {k: forward[k] for k in ("seqid", "start", "end", "strand")},
                             "reverse": {k: reverse[k] for k in ("seqid", "start", "end", "strand")}}
    return list(products.values())


def bowtie2_specificity(path: Path, request: PrimerRequest, target: Dict[str, Any],
                       candidates: List[Dict[str, Any]]) -> Dict[str, Any]:
    binary, prefix = require_bowtie2(request.species_id)
    unique_primers = list(dict.fromkeys(
        str(candidate.get(key) or candidate[key.replace("_core", "")]).upper()
        for candidate in candidates for key in ("left_core_sequence", "right_core_sequence")
    ))
    expected = [expected_genome_pair(target, candidate) for candidate in candidates]
    requested_product_max = max(int(bounds.split("-")[1]) for bounds in request.product_size_range.split(","))
    # Widen to include known targets when policy permits; expose longer targets as out of range.
    desired_span = max([request.max_offtarget_product_size, requested_product_max]
                       + [int(pair["span"]) for pair in expected if pair])
    screening_max = min(PRIMER_LIMITS.max_offtarget_product_size, desired_span)
    single_fasta, left_fasta, right_fasta = [path / name for name in
                                           ("specificity_primers.fa", "specificity_left.fa", "specificity_right.fa")]
    single_sam, pair_sam = path / "specificity_primers.sam", path / "specificity_pairs.sam"
    primer_ids = {f"primer{index}": sequence for index, sequence in enumerate(unique_primers)}
    with single_fasta.open("wt", encoding="utf-8") as handle:
        for identifier, sequence in primer_ids.items():
            write_fasta_record(handle, identifier, sequence)
    with left_fasta.open("wt", encoding="utf-8") as left_handle, right_fasta.open("wt", encoding="utf-8") as right_handle:
        for index, candidate in enumerate(candidates):
            write_fasta_record(left_handle, f"pair{index}", str(candidate.get("left_core_sequence") or candidate["left_sequence"]))
            write_fasta_record(right_handle, f"pair{index}", str(candidate.get("right_core_sequence") or candidate["right_sequence"]))
    common = [binary, "-x", str(prefix), "-f", "-k", "2", "--end-to-end", "-N", "0",
              "-L", str(max(4, min(20, min(map(len, unique_primers))))), "--score-min", "L,0,0",
              "--no-unal", "--quiet", "-p", str(PRIMER_BOWTIE2_THREADS)]
    timeout = PRIMER_LIMITS.specificity_timeout_seconds
    subprocess.run(common + ["-U", str(single_fasta), "-S", str(single_sam)], check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
    subprocess.run(common + ["-1", str(left_fasta), "-2", str(right_fasta), "-S", str(pair_sam),
                             "--no-mixed", "--no-discordant", "-I", "0", "-X", str(screening_max), "--fr"],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
    single_records, pair_records = exact_sam_records(single_sam), exact_sam_records(pair_sam)
    primer_counts = {sequence: min(2, len({binding_key(r) for r in single_records.get(identifier, [])}))
                     for identifier, sequence in primer_ids.items()}
    products = [sam_pair_products(pair_records.get(f"pair{index}", [])) for index in range(len(candidates))]
    return {"engine": "bowtie2_exact_k2", "primer_counts": primer_counts, "products": products,
            "expected_pairs": expected, "screening_max_product_size": screening_max}


def add_specificity(path: Path, request: PrimerRequest, target: Dict[str, Any], candidates: List[Dict[str, Any]]) -> None:
    if not candidates:
        return
    result = bowtie2_specificity(path, request, target, candidates)
    for index, primer in enumerate(candidates):
        left = str(primer.get("left_core_sequence") or primer["left_sequence"])
        right = str(primer.get("right_core_sequence") or primer["right_sequence"])
        left_count, right_count = result["primer_counts"].get(left, 0), result["primer_counts"].get(right, 0)
        products, expected = result["products"][index], result["expected_pairs"][index]
        count = min(2, len(products))
        on_target = sum(binding_key(p["forward"]) == binding_key(expected["forward"])
                        and binding_key(p["reverse"]) == binding_key(expected["reverse"])
                        for p in products) if expected else None
        beyond_limit = bool(expected and expected["span"] > result["screening_max_product_size"])
        if beyond_limit:
            status = "target_out_of_range"
        elif not target.get("seqid"):
            status = "genome_matches_only"
        elif expected is None:
            status = "junction_genome_screen"
        elif on_target:
            status = "target_confirmed" if count == 1 else "off_target_detected"
        else:
            status = "off_target_detected" if count else "target_not_detected"
        primer["specificity"] = {
            "engine": result["engine"], "scope": "exact_genome", "hit_count_cap": 2,
            "left_hit_count": left_count, "right_hit_count": right_count,
            "left_hits_capped": left_count >= 2, "right_hits_capped": right_count >= 2,
            "pair_amplicon_count": count, "pair_count_capped": count >= 2,
            "on_target_product_count": on_target,
            "off_target_product_count": len(products) - on_target if on_target is not None else None,
            "off_target_count_is_lower_bound": count >= 2,
            "left_match_class": specificity_class(left_count), "right_match_class": specificity_class(right_count),
            "pair_amplicon_class": specificity_class(count), "status": status,
            "screening_max_product_size": result["screening_max_product_size"],
            "expected_genomic_product_size": expected["span"] if expected else None,
            "products": products,
        }


def visible_candidate_key(candidate: Dict[str, Any]) -> Tuple[int, int, int, float, int]:
    specificity = candidate.get("specificity") or {}
    raw_pair_count = specificity.get("pair_amplicon_count")
    pair_count = int(raw_pair_count) if raw_pair_count is not None else 1
    left_count = int(specificity.get("left_hit_count", 1))
    right_count = int(specificity.get("right_hit_count", 1))
    status = specificity.get("status")
    if status == "junction_genome_screen":
        # Lack of a contiguous genomic product is desirable for a junction primer,
        # but says nothing about its specificity among spliced transcripts.
        return (0 if pair_count == 0 else 1, pair_count, left_count + right_count,
                float(candidate.get("penalty", 0.0)), int(candidate.get("product_size", 0)))
    if status and status not in {"target_confirmed", "genome_matches_only"}:
        return (2, pair_count, left_count + right_count,
                float(candidate.get("penalty", 0.0)), int(candidate.get("product_size", 0)))
    return (
        0 if pair_count == 1 else 1,
        0 if left_count == 1 and right_count == 1 else 1,
        pair_count,
        float(candidate.get("penalty", 0.0)),
        int(candidate.get("product_size", 0)),
    )


def primer_pair_geometry(candidate: Dict[str, Any]) -> Tuple[int, int, int]:
    left_start = int(candidate["left_position"]["offset"])
    right_end = int(candidate["right_position"]["offset"])
    right_start = right_end - int(candidate["right_position"]["length"]) + 1
    return left_start, right_start, right_end


def spatially_redundant_primer_pair(candidate: Dict[str, Any], existing: Dict[str, Any]) -> bool:
    """Reject gDNA pairs that reuse a primer site or represent the same product region."""
    left_start, right_start, right_end = primer_pair_geometry(candidate)
    existing_left, existing_right_start, existing_right_end = primer_pair_geometry(existing)

    # A repeated or nearly repeated primer-binding site does not add an independent assay option.
    if abs(left_start - existing_left) <= 30 or abs(right_start - existing_right_start) <= 30:
        return True

    overlap = max(0, min(right_end, existing_right_end) - max(left_start, existing_left) + 1)
    shorter_product = min(right_end - left_start + 1, existing_right_end - existing_left + 1)
    return overlap / max(shorter_product, 1) >= 0.5


def select_visible_candidates(
    candidates: List[Dict[str, Any]],
    enforce_spatial_diversity: bool = False,
) -> List[Dict[str, Any]]:
    selected: List[Dict[str, Any]] = []
    for candidate in sorted(candidates, key=visible_candidate_key):
        left_offset = int(candidate["left_position"]["offset"])
        right_offset = int(candidate["right_position"]["offset"])
        if any(
            abs(left_offset - int(existing["left_position"]["offset"])) <= 18
            and abs(right_offset - int(existing["right_position"]["offset"])) <= 18
            for existing in selected
        ):
            continue
        if enforce_spatial_diversity and any(spatially_redundant_primer_pair(candidate, existing) for existing in selected):
            continue
        selected.append(candidate)
        if len(selected) >= PRIMER_VISIBLE_PAIR_LIMIT:
            break
    for rank, candidate in enumerate(selected, start=1):
        candidate["rank"] = rank
    return selected
