"""`synaptik data`: download reference transcriptomes for off-target screening."""

from __future__ import annotations

import os
import urllib.request
from pathlib import Path

from ..main import Fmt, header

TRANSCRIPTOME_URLS: dict[str, tuple[str, str]] = {
    "human": (
        "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
    ),
    "mouse": (
        "https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
        "Mus_musculus.GRCm39.cdna.all.fa.gz",
    ),
}


def default_data_dir() -> Path:
    return Path.home() / ".synaptik" / "data"


def fetch_transcriptome(
    species: str = "human",
    out_path: str | Path | None = None,
    force: bool = False,
) -> Path:
    """Download the Ensembl cDNA transcriptome FASTA (.fa.gz) atomically."""
    key = species.strip().lower()
    if key not in TRANSCRIPTOME_URLS:
        raise ValueError(f"Unsupported species '{species}'. Choose from: {', '.join(sorted(TRANSCRIPTOME_URLS))}.")
    url, filename = TRANSCRIPTOME_URLS[key]
    dest = Path(out_path) if out_path else (default_data_dir() / filename)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        return dest

    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.tmp")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "SYNAPTIK-Research/0.9"})
        with urllib.request.urlopen(req, timeout=120) as resp, open(tmp, "wb") as out:
            while True:
                chunk = resp.read(1024 * 256)
                if not chunk:
                    break
                out.write(chunk)
        os.replace(tmp, dest)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
    return dest


def cmd_data(args) -> int:
    if args.data_action == "fetch-transcriptome":
        header(f"Fetching {args.species} reference transcriptome (Ensembl cDNA)")
        dest = fetch_transcriptome(species=args.species, out_path=args.out, force=args.force)
        size_mb = dest.stat().st_size / (1024 * 1024)
        print(Fmt.green(f"  Ready: {dest} ({size_mb:.1f} MB)"))
        print(Fmt.dim(f"  Use with: synaptik aso --sequence <SEQ> --transcriptome \"{dest}\""))
        return 0
    return 1


def register(sub, common) -> None:
    s = sub.add_parser(
        "data", parents=[common],
        help="manage local reference transcriptomes for ASO off-target screening",
    )
    dsub = s.add_subparsers(dest="data_action", required=True)
    ft = dsub.add_parser(
        "fetch-transcriptome", parents=[common],
        help="download Ensembl GRCh38 (human) or GRCm39 (mouse) cDNA FASTA",
    )
    ft.add_argument("--species", choices=sorted(TRANSCRIPTOME_URLS), default="human")
    ft.add_argument("--out", help="destination file path (default: ~/.synaptik/data/...)")
    ft.add_argument("--force", action="store_true", help="re-download even if file exists")
    ft.set_defaults(func=cmd_data)
