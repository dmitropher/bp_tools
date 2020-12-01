#!/usr/bin/env python3
import json
import click


def ss_elements_from_param_set(param_set):
    """
    """


def bp_string_from_param_set(param_set, abego=False):
    """
    """
    first = True
    helix_type = "H" if abego else "HA"
    sheet_type = "E" if abego else "ED"
    loop_type = "L" if abego else "LD"
    tmpType = "LD"
    ss_elements = ss_elements_from_param_set(param_set)
    bp_string = ""
    for ssType, ssLength in ss_elements:
        if ssType.lower() == "h":
            tmpType = helix_type
        if ssType.lower() == "l":
            tmpType = loop_type
        if ssType.lower() == "e":
            tmpType = sheet_type
        if first:

            line = "{} {} {}".format(1, "A", tmpType) + "\n"
            bp_string.append(line)
            first = False
            for ii in range(1, ssLength):
                line = "{} {} {}".format(0, "x", tmpType) + "\n"
                bp_string.append(line)
        else:
            for ii in range(0, ssLength):
                line = "{} {} {}".format(0, "x", tmpType) + "\n"
                bp_string.append(line)


def build_files_from_param_set(param_set):
    """
    """
    first = True
    bp_string = bp_string_from_param_set(param_set)


@click.command()
@click.argument("entry")
@click.argument("json_path")
@click.option("-o", "--output-dir", default=".")
@click.option(
    "-r",
    "--rosetta-flags-file",
    "rosetta_flags_file",
    default="",
    help="Optional: include a rosetta flags file",
)
def main(
    entry,
    json_path,
    output_dir=".",
    append=False,
    abego=False,
    rosetta_flags_file="",
):
    """
    Reads a particular entry from a store and builds bp and cst files (copies others)
    """
    with open(json_path, "r") as f:
        param_list = json.read(f)

    param_set = param_list[entry]
    bp_string, cst_string = build_files_from_param_set(param_set, abego=abego)

    with open(f"{output_dir}/design.bp", "w") as f:
        f.write(bp_string)
    with open(f"{output_dir}/design.cst", "w") as f:
        f.write(cst_string)


if __name__ == "__main__":
    main()
