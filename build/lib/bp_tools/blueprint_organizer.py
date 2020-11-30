#!/usr/bin/env python3
from itertools import product

from argparse import ArgumentParser

import json

parser = ArgumentParser()


import click


@click.command()
@click.option("-o", "--output-dir", default=".")
@click.option("-f", "--fragment-file", default="")
@click.option(
    "-s",
    "--struct-params",
    "struct_params",
    default=[],
    show_default=True,
    multiple=True,
)
@click.option(
    "-w",
    "--write-frag-file",
    "write_frag_file",
    default=False,
    show_default=True,
)
@click.option(
    "-l",
    "--write-param-sampler-file",
    "write_param_sampler_file",
    default="",
    show_default=True,
)
@click.option(
    "-e",
    "--extra-files-dir",
    "extra_files_dir",
    default=".",
    show_default=True,
)
@click.option("-p", "--extra-pdb", "extra_pdb", default="", show_default=True)
@click.option(
    "-a/ ",
    "--append-mode/--prepend-mode",
    "append",
    default=False,
    show_default=True,
)
@click.option(
    "-r",
    "--rosetta-flags-file",
    "rosetta_flags_file",
    default="",
    help="Optional: include a rosetta flags file",
)
def main(
    output_dir=".",
    struct_params=[],
    fragment_file="",
    extra_files_dir=".",
    extra_pdb="",
    append=False,
    write_param_sampler_file="",
    rosetta_flags_file="",
):
    ""
    if struct_params and fragment_file:
        raise ValueError(
            "either a fragment file or struct params must be given, but not both"
        )
    sse_sampler_list = (
        build_from_params(struct_params)
        if struct_params
        else build_from_file(fragment_file)
    )

    # dump the sampler in case we have various options for param interpretation
    # For instance, inserting different fragments from a pdb into an existing param set
    if write_param_sampler_file:
        with open(write_param_sampler_file, "w") as f:
            json.dump(sse_sampler_list, f)

    fragerator = product(*(f.get_ss_elements_list() for f in sse_sampler_list))

    dump_list = []

    # run through all the permutations and save each collection of frags as a list
    for ss_elements in fragerator:
        print(ss_elements)
        elem_dict_list = [sse.to_dict() for sse in ss_elements]
        dump_list.append(elem_dict_list)

    with open("frag_params.json", "w") as f:
        json.dump(dump_list, f)


if __name__ == "__main__":
    main()
