# -*- coding: utf-8 -*-
"""Script to convert FASTA files from DNA to RNA format.

@author: Laura C. Kühle
"""

if "snakemake" in locals():
    smk = snakemake

    def main() -> None:
        """Convert data via Snakemake."""
        print(f"Query: {smk.input[0]}")
        print(f"Output: {smk.output[0]}")
        print(f"Mode: {smk.params['mode']}\n")

        print("Convert FASTA...")
        convert(
            qry_path=smk.input[0], output_path=smk.output[0], mode=smk.params["mode"]
        )
        print("Conversion completed!\n")


def convert(qry_path: str, output_path: str, mode: str) -> None:
    with open(qry_path, "r", encoding="utf-8") as qry_file:
        with open(output_path, "w", encoding="utf-8") as output_file:
            for qry_line in qry_file.readlines():
                # Convert line (if not header)
                if not qry_line.startswith(">"):
                    match mode:
                        case "ref":
                            seq = convert_ref(seq=qry_line.rstrip("\n"))
                        case "qry":
                            seq = convert_qry(seq=qry_line.rstrip("\n"))
                        case _:
                            raise Exception(f"Mode {mode} not implemented.")

                    qry_line = seq + "\n"

                # Write updated line to file
                output_file.write(qry_line)


def convert_qry(seq: str):
    return "".join([convert_to_unmodified(base) for base in seq])


def convert_to_unmodified(base: str):
    if base in ["U", "A", "C", "G"]:
        return base
    return "N"


def convert_ref(seq: str):
    return "".join([convert_to_rna(base) for base in seq])


def convert_to_rna(base: str):
    if base == "T":
        return "U"
    return base


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
