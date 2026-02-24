# -*- coding: utf-8 -*-
"""Script to update reference FASTA based on additional sources.

@author: Laura C. Kühle
"""

import re
import polars as pl


MOD_LIST = pl.read_excel("assets/modification_dictionary.xlsx")


if "snakemake" in locals():
    smk = snakemake

    def main() -> None:
        """Update reference via Snakemake."""
        print(f"Reference: {smk.input['ref']}")
        print(f"Query: {smk.input['qry']}")
        print(f"Alignment: {smk.input['align']}")
        print(f"Output: {smk.output['ref']}")
        print(f"Meta: {smk.output['meta']}")

        print("Update reference FASTA...")
        update(
            ref_path=smk.input["ref"],
            qry_path=smk.input["qry"],
            align_path=smk.input["align"],
            output_path=smk.output["ref"],
            meta_path=smk.output["meta"],
        )
        print("Update completed!\n")


def update(ref_path: str, qry_path: str, align_path: str, output_path: str,
           meta_path: str) -> None:
    ref_seqs = get_sequences(input_path=ref_path)
    ref_header = list(ref_seqs.keys())[0]
    ref_seq = ref_seqs[ref_header]

    qry_seqs = get_sequences(input_path=qry_path)

    alignments = get_alignments(input_path=align_path)

    update_dict = {}
    for qry_id in alignments.keys():
        ref_seq, update_dict = update_ref(
            ref_pos=alignments[qry_id][0],
            cigar=alignments[qry_id][1],
            ref_seq=ref_seq,
            qry_seq=qry_seqs[qry_id],
            col_name=qry_id.split("_")[-1],
            update_dict=update_dict,
        )

    # Write updated line to file
    with open(output_path, "w", encoding="utf-8") as output_file:
        output_file.write(">" + ref_header + "\n")
        output_file.write(ref_seq)

    # Write update meta to file
    pos_list = update_dict.keys()
    pos_list = sorted(pos_list)
    meta = pl.from_dict({"pos": pos_list, "name": [update_dict[pos] for pos in pos_list]})
    meta.write_csv(meta_path, separator="\t")


def update_ref(ref_pos, cigar, ref_seq, qry_seq, col_name, update_dict):
    qry_pos = 0
    ref_pos -= 1
    new_seq = ref_seq[:ref_pos]
    for val in flatten_cigar(cigar):
        match val:
            case "I":
                qry_pos += 1
            case "D":
                new_seq += ref_seq[ref_pos]
                ref_pos += 1
            case "M":
                if qry_seq[qry_pos] in ["U", "A", "C", "G"]:
                    new_seq += ref_seq[ref_pos]
                else:
                    mod_name = qry_seq[qry_pos]
                    mod_base = (
                        MOD_LIST.filter(pl.col(col_name) == mod_name)
                        .get_column("unmodified_base")
                        .to_list()[0]
                    )
                    mod_id = (
                        MOD_LIST.filter(pl.col(col_name) == mod_name)
                        .get_column("short_name")
                        .to_list()[0]
                    )
                    # print(mod_name, mod_base, mod_id)
                    update_dict[ref_pos+1] = mod_id
                    new_seq += mod_base

                qry_pos += 1
                ref_pos += 1
            case _:
                print("Yeah, that is a problem: ", val)

    new_seq += ref_seq[ref_pos:]

    return new_seq, update_dict


def flatten_cigar(cigar):
    # Flatten blocks in CIGAR string
    entries = re.split(r"([MIDNSHP=X])", cigar)[:-1]
    return "".join(
        [
            int(entries[2 * idx]) * entries[2 * idx + 1]
            for idx in range(len(entries) // 2)
        ]
    )


def get_alignments(input_path: str) -> dict:
    align_dict = {}
    with open(input_path, "r", encoding="utf-8") as input_file:
        for line in input_file.readlines():
            if line.startswith("@"):
                continue

            entries = line.rstrip("\n").split("\t")
            align_dict[entries[0]] = (int(entries[3]), entries[5])

    return align_dict


def get_sequences(input_path: str) -> dict:
    header = ""
    seq_dict = {}
    with open(input_path, "r", encoding="utf-8") as input_file:
        line = input_file.readline()
        while line != "":
            if line.startswith(">"):
                header = line.lstrip(">").rstrip("\n")
            else:
                if not header in seq_dict.keys():
                    seq_dict[header] = ""
                seq_dict[header] += line.rstrip("\n")

            line = input_file.readline()

    return seq_dict


if __name__ == "__main__":
    if "snakemake" in locals():
        main()
    else:
        print("Not Defined.")
