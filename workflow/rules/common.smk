import csv


def metadata_proteome_rows():
    rows = []
    path = config["metadata"]["species"]
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            protein_fasta = row.get("protein_fasta", "").strip()
            if protein_fasta and protein_fasta.upper() != "NA":
                rows.append(row)
    return rows


def expansion_candidate_rows():
    rows = []
    path = config["metadata"].get("expansion_candidates")
    if not path:
        return rows
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            protein_fasta = row.get("protein_fasta", "").strip()
            if (
                row.get("download_status", "").strip() == "downloaded"
                and row.get("recommended_use", "").strip() == "outgroup_or_sensitivity_analysis"
                and protein_fasta
                and protein_fasta.upper() != "NA"
            ):
                rows.append(row)
    return rows


PROTEOME_ROWS = metadata_proteome_rows()
PROTEOME_FASTA_BY_SPECIES = {
    row["species_id"]: row["protein_fasta"]
    for row in PROTEOME_ROWS
}
PROTEOME_SPECIES = sorted(PROTEOME_FASTA_BY_SPECIES)

EXPANSION_CANDIDATE_ROWS = expansion_candidate_rows()
EXPANSION_FASTA_BY_SPECIES = {
    row["species_id"]: row["protein_fasta"]
    for row in EXPANSION_CANDIDATE_ROWS
}
EXPANSION_SPECIES = sorted(EXPANSION_FASTA_BY_SPECIES)
ORTHOFINDER_EXPANSION_SPECIES = sorted(PROTEOME_SPECIES + EXPANSION_SPECIES)
