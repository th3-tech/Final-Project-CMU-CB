#!/usr/bin/env python3
"""Convert protein structures to FASTA; keep coordinate CIF files; remove redundant formats."""

import gzip
import re
import shutil
from pathlib import Path

PROTEINS_DIR = Path(__file__).resolve().parent / "Proteins"
FASTA_DIR = PROTEINS_DIR / "fasta"
CIF_DIR = PROTEINS_DIR / "cif"
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K",
}


def read_cif_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", errors="replace") as handle:
            return handle.read()
    return path.read_text(errors="replace")


def parse_semicolon_value(lines: list[str], start: int) -> tuple[str, int]:
    value_lines = []
    index = start + 1
    while index < len(lines):
        line = lines[index]
        if line.strip() == ";":
            break
        value_lines.append(line.strip())
        index += 1
    return "".join(value_lines), index + 1


def sequences_from_entity_poly(text: str) -> list[tuple[str, str]]:
    lines = text.splitlines()
    sequences: list[tuple[str, str]] = []

    index = 0
    while index < len(lines):
        line = lines[index]
        if line.strip() == "loop_":
            headers: list[str] = []
            index += 1
            while index < len(lines) and lines[index].startswith("_"):
                headers.append(lines[index].strip().split()[0])
                index += 1

            if "_entity_poly.pdbx_seq_one_letter_code" not in headers:
                continue

            seq_col = headers.index("_entity_poly.pdbx_seq_one_letter_code")
            entity_col = headers.index("_entity_poly.entity_id") if "_entity_poly.entity_id" in headers else None
            strand_col = (
                headers.index("_entity_poly.pdbx_strand_id")
                if "_entity_poly.pdbx_strand_id" in headers
                else None
            )

            while index < len(lines):
                row_line = lines[index]
                stripped = row_line.strip()
                if not stripped or stripped.startswith("#") or stripped.startswith("_") or stripped == "loop_":
                    break

                values: list[str] = []
                if stripped.startswith(";"):
                    sequence, index = parse_semicolon_value(lines, index)
                    values = [sequence]
                else:
                    parts = re.findall(r"'(?:\\'|[^'])*'|\"(?:\\\"|[^\"])*\"|[^\s]+", row_line)
                    values = [part.strip("'\"") for part in parts]
                    index += 1

                if len(values) <= seq_col:
                    continue

                sequence = re.sub(r"\s+", "", values[seq_col])
                if not sequence or not re.fullmatch(r"[A-Za-z]+", sequence):
                    continue

                label = "entity"
                if entity_col is not None and len(values) > entity_col:
                    label = values[entity_col]
                if strand_col is not None and len(values) > strand_col and values[strand_col] not in ("?", ""):
                    label = values[strand_col]

                sequences.append((label, sequence.upper()))
            continue

        index += 1

    return sequences


def sequences_from_poly_seq_scheme(text: str) -> list[tuple[str, str]]:
    lines = text.splitlines()
    chains: dict[str, list[tuple[int, str]]] = {}

    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue

        headers: list[str] = []
        index += 1
        while index < len(lines) and lines[index].startswith("_"):
            headers.append(lines[index].strip().split()[0])
            index += 1

        if "_pdbx_poly_seq_scheme.mon_id" not in headers:
            continue

        mon_col = headers.index("_pdbx_poly_seq_scheme.mon_id")
        seq_col = headers.index("_pdbx_poly_seq_scheme.seq_id")
        chain_col = (
            headers.index("_pdbx_poly_seq_scheme.pdb_strand_id")
            if "_pdbx_poly_seq_scheme.pdb_strand_id" in headers
            else None
        )

        while index < len(lines):
            row_line = lines[index].strip()
            if not row_line or row_line.startswith("#") or row_line.startswith("_") or row_line == "loop_":
                break

            parts = re.findall(r"'(?:\\'|[^'])*'|\"(?:\\\"|[^\"])*\"|[^\s]+", lines[index])
            values = [part.strip("'\"") for part in parts]
            index += 1

            if len(values) <= max(mon_col, seq_col):
                continue

            residue = values[mon_col].upper()
            if residue not in AA3TO1:
                continue

            chain = values[chain_col] if chain_col is not None and len(values) > chain_col else "A"
            chains.setdefault(chain, []).append((int(values[seq_col]), AA3TO1[residue]))

    sequences: list[tuple[str, str]] = []
    for chain, residues in sorted(chains.items()):
        residues.sort(key=lambda item: item[0])
        sequences.append((chain, "".join(residue for _, residue in residues)))
    return sequences


def sequences_from_pdb(text: str) -> list[tuple[str, str]]:
    chains: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line.startswith("SEQRES"):
            continue
        chain = line[11]
        residues = line.split()[4:]
        chains.setdefault(chain, []).extend(AA3TO1.get(residue.upper(), "X") for residue in residues)

    return [(chain, "".join(residues)) for chain, residues in sorted(chains.items())]


def choose_primary_chain(sequences: list[tuple[str, str]]) -> tuple[str, str]:
    if not sequences:
        raise ValueError("no sequences")

    for label, sequence in sequences:
        if label == "A":
            return label, sequence

    label, sequence = max(sequences, key=lambda item: len(item[1]))
    return label, sequence


def extract_sequences(path: Path) -> list[tuple[str, str]]:
    text = read_cif_text(path)

    sequences = sequences_from_poly_seq_scheme(text)
    if sequences:
        return sequences

    sequences = sequences_from_entity_poly(text)
    if sequences:
        return sequences

    if path.suffixes[-2:] == [".pdb", ".gz"]:
        with gzip.open(path, "rt", errors="replace") as handle:
            text = handle.read()
    elif path.suffix == ".pdb":
        text = path.read_text(errors="replace")
    else:
        return []

    return sequences_from_pdb(text)


def preferred_cif_path(base_id: str, cif_paths: list[Path]) -> Path:
    exact = [path for path in cif_paths if path.name == f"{base_id}.cif.gz"]
    if exact:
        return exact[0]

    assembly_one = [path for path in cif_paths if path.name == f"{base_id}-assembly1.cif.gz"]
    if assembly_one:
        return assembly_one[0]

    return sorted(cif_paths, key=lambda path: path.name)[0]


def write_fasta(path: Path, protein_id: str, label: str, sequence: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{protein_id}|chain={label}|length={len(sequence)}\n{sequence}\n")


def write_cif(source_cif: Path, destination: Path) -> None:
    if source_cif.suffix == ".gz":
        with gzip.open(source_cif, "rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(source_cif, destination)


def protein_base_id(filename: str) -> str:
    name = filename.replace(".cif.gz", "")
    if "-assembly" in name:
        name = name.split("-assembly", 1)[0]
    return name.lower()


def main() -> None:
    coordinate_cifs = [
        path
        for path in PROTEINS_DIR.rglob("*.cif.gz")
        if not path.name.endswith("-sf.cif.gz")
    ]

    grouped: dict[str, list[Path]] = {}
    for path in coordinate_cifs:
        grouped.setdefault(protein_base_id(path.name), []).append(path)

    staging_fasta = PROTEINS_DIR / ".staging-fasta"
    staging_cif = PROTEINS_DIR / ".staging-cif"
    for staging_dir in (staging_fasta, staging_cif):
        if staging_dir.exists():
            shutil.rmtree(staging_dir)
        staging_dir.mkdir()

    converted = 0
    skipped: list[str] = []

    for protein_id, cif_paths in sorted(grouped.items()):
        source_cif = preferred_cif_path(protein_id, cif_paths)
        sequences = extract_sequences(source_cif)
        if not sequences:
            skipped.append(protein_id)
            continue

        label, sequence = choose_primary_chain(sequences)
        write_cif(source_cif, staging_cif / f"{protein_id}.cif")
        write_fasta(staging_fasta / f"{protein_id}.fasta", protein_id, label, sequence)
        converted += 1

    FASTA_DIR.mkdir(exist_ok=True)
    CIF_DIR.mkdir(exist_ok=True)

    for child in FASTA_DIR.iterdir():
        child.unlink()
    for child in CIF_DIR.iterdir():
        child.unlink()

    for item in staging_fasta.iterdir():
        shutil.move(str(item), FASTA_DIR / item.name)
    for item in staging_cif.iterdir():
        shutil.move(str(item), CIF_DIR / item.name)

    staging_fasta.rmdir()
    staging_cif.rmdir()

    for child in PROTEINS_DIR.iterdir():
        if child.name in {".staging-fasta", ".staging-cif", "fasta", "cif"}:
            continue
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()

    chou_data = Path(__file__).resolve().parent / "Chou-Fasman" / "data"
    for test_file in ("polyA.fasta", "polyV.fasta", "1oai_test.fasta", "6jhz.fasta"):
        test_path = chou_data / test_file
        if test_path.exists():
            test_path.unlink()

    print(f"Converted {converted} proteins to FASTA and kept matching CIF files.")
    if skipped:
        print(f"Skipped {len(skipped)} files with no extractable sequence: {', '.join(skipped[:10])}")
        if len(skipped) > 10:
            print(f"... and {len(skipped) - 10} more")


if __name__ == "__main__":
    main()
