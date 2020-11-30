#!/usr/bin/env python3
from itertools import product

import pyrosetta

from argparse import ArgumentParser
import os
import json

parser = ArgumentParser()


import click


def read_flag_file(filename):
    """
    Reads the flag file, ignoring comment lines


    """
    with open(filename, "r") as myfile:
        lines = myfile.read().splitlines()
    # filter the lines
    lines = [l for l in lines if l.startswith("-")]
    return " ".join(lines)


def run_pyrosetta_with_flags(flags_file_path, mute=False):
    if not flags_file_path:
        pyrosetta.init("-mute all " if mute else "", silent=mute)
        return
    flags = read_flag_file(flags_file_path)
    flags_str = " ".join(flags.replace("\n", " ").split())
    pyrosetta.init(
        f"-mute all {flags_str}" if mute else flags_str, silent=mute
    )


def safe_load_pdb(pdb, rosetta_flags_file=""):
    run_pyrosetta_with_flags(rosetta_flags_file, mute=False)
    try:
        return pyrosetta.pose_from_pdb(pdb)
    except RuntimeError as e:
        print(e)
        print(f"unable to load: {pdb}")
        return


class SecondaryStructElementSampler(object):
    """
    A container for describing the fragment type and sizes to use for b-print

    This container is used as an interface for sampling various lengths
    """

    def __init__(
        self, dssp_type, min_size, max_size, repeat_dist=0, repeat_dist_cst=0
    ):
        self.dssp_type = dssp_type
        self.min_size = min_size
        self.max_size = max_size
        self.repeat_dist = repeat_dist
        self.repeat_dist_cst = repeat_dist_cst

    def to_dict(self):
        return {
            "dssp_type": self.dssp_type,
            "min_size": self.min_size,
            "max_size": self.max_size,
            "repeat_dist": self.repeat_dist,
            "repeat_dist_cst": self.repeat_dist_cst,
        }

    def tuple_generator(self):
        return (
            (self.dssp_type, i, self.repeat_dist, self.repeat_dist_cst)
            for i in range(self.min_size, self.max_size + 1)
        )

    def from_dict(cls, dict):
        return cls(**dict)

    def __repr__(self):
        return f"""SecondaryStructElementSampler(**{self.to_dict()})"""

    def __str__(self):
        return str(
            "dssp_type"
            + self.dssp_type
            + "\n"
            + "min_size"
            + self.min_size
            + "\n"
            + "max_size"
            + self.max_size
            + "\n"
            + "repeat_dist"
            + self.repeat_dist
            + "\n"
            + "repeat_dist_cst"
            + self.repeat_dist_cst
            + "\n"
        )

    def get_ss_elements_list(self):
        return [
            SecondaryStructElement(
                self.dssp_type, i, self.repeat_dist, self.repeat_dist
            )
            for i in range(self.min_size, self.max_size + 1)
        ]


class SecondaryStructElement(object):
    """
    A container for describing the specific fragment type to use for b-print

    This container is a convenience interface for saving run params
    """

    def __init__(self, dssp_type, size, repeat_dist=0, repeat_dist_cst=0):
        self.dssp_type = dssp_type
        self.size = size
        self.repeat_dist = repeat_dist
        self.repeat_dist_cst = repeat_dist_cst

    def to_dict(self):
        return {
            "dssp_type": self.dssp_type,
            "size": self.size,
            "repeat_dist": self.repeat_dist,
            "repeat_dist_cst": self.repeat_dist_cst,
        }

    def from_dict(cls, dict):
        return cls(**dict)

    def __repr__(self):
        return f"""SecondaryStructElement(**{self.to_dict()})"""

    def __str__(self):
        return str(
            "dssp_type"
            + self.dssp_type
            + "\n"
            + "size"
            + self.size
            + "\n"
            + "repeat_dist"
            + self.repeat_dist
            + "\n"
            + "repeat_dist_cst"
            + self.repeat_dist_cst
            + "\n"
        )


def get_default_xml():
    """
    returns a hardcoded xml
    """
    return """<ROSETTASCRIPTS>
    <TASKOPERATIONS>
    </TASKOPERATIONS>
    <SCOREFXNS>
    <ScoreFunction name="sfn_centroid" weights="abinitio_remodel_cen.wts">
        <Reweight scoretype="sheet" weight="4"/>
        <Reweight scoretype="ss_pair" weight="1"/>
    </ScoreFunction>
    <ScoreFunction name="sfn_motif" weights="empty">
    <Reweight scoretype="cen_pair_motifs" weight="1"/>
    </ScoreFunction>
    <ScoreFunction name="sfn_motif_degree" weights="empty">
    <Reweight scoretype="cen_pair_motif_degree" weight="1"/>
    </ScoreFunction>
    <ScoreFunction name="VDW" weights="empty">
    <Reweight scoretype="vdw" weight="1"/>
    </ScoreFunction>
    </SCOREFXNS>
    <FILTERS>
    <worst9mer name="worst9mer_h" threshold="0.15" only_helices="true"/>
    <ScoreType name="VDW" scorefxn="VDW" threshold="100" confidence="1" />
    <ScoreType name="motif_score" scorefxn="sfn_motif" threshold="-1" confidence="1" />
    <ScoreType name="motif_degree_score" scorefxn="sfn_motif_degree" threshold="-0.3" confidence="1" />
    <SSDegree name="ss_degree_avg" report_avg="true" ignore_terminal_ss="2"/>
    <SSDegree name="ss_degree_worst" report_avg="false" ignore_terminal_ss="2"/>
    <RepeatParameter name="radius" param_type="radius" numb_repeats="4" min="0" max="99999999" confidence="1"/> #replace with correct number of repeats
    <RepeatParameter name="rise" param_type="rise" numb_repeats="4" min="0" max="99999999" confidence="1"/>
    <RepeatParameter name="omega" param_type="omega" numb_repeats="4" min="0" max="99999999" confidence="1"/>
    </FILTERS>
    <MOVERS>
    <RemodelMover name="remodel_mover" blueprint="design.blueprint"/>
    </MOVERS>
    <PROTOCOLS>
    <Add mover_name="remodel_mover"/>
    <Add filter_name="VDW"/>
    <Add filter_name="worst9mer_h"/>
    <Add filter_name="motif_score"/>
    <Add filter_name="motif_degree_score"/>
    <Add filter_name="ss_degree_worst"/>
    <Add filter_name="radius"/>
    <Add filter_name="rise"/>
    <Add filter_name="omega"/>
    </PROTOCOLS>
    <OUTPUT scorefxn="sfn_centroid"/>
    </ROSETTASCRIPTS>"""


def prepare_design_dir(
    dirname,
    ss_elements,
    extra_files_dir,
    extra_pose=None,
    append=False,
    abego=False,
):  # xml_str=""):
    """
    prepares a directory to run the rosetta_scripts executable

    ss_elements should be a list of tuples (dssp_type,size,lattice_space,cst_tolerance)

    setting lattice space to 0 will create no constraints for that element
    """

    # forgive this ugly str parse
    name = "_".join(
        ["".join((type, str(size))) for type, size, lat, cst in ss_elements]
    )
    repeat_size = sum(size for type, size, lat, cst in ss_elements)
    if dirname != "":
        path_name = "{}/{}".format(dirname, name)
    if not os.path.exists(path_name):
        os.makedirs(path_name)

    # TODO remove
    copy_necessary_files(path_name, extra_files_dir)
    blueprint_elements = [(type, size) for type, size, lat, cst in ss_elements]
    create_blueprint(
        path_name + "/design.blueprint",
        blueprint_elements,
        extra_pose=extra_pose,
        append=append,
        abego=abego,
    )

    offset = 1
    cst_string = ""
    for type, size, lat, cst in ss_elements:
        if lat:
            cst_string += create_atom_pair_cst_string(
                range(offset, offset + size),
                lat,
                repeat_size,
                cst,
                sheet_mode=(type.lower() == "e"),
            )
        offset += size
    with open(path_name + "/lattice_csts.cst", "w") as f:
        f.write(cst_string)
    add_flags(path_name, name, repeat_size)
    # fl = open(path_name + "/design.xml", "w")
    # if not xml_str:
    #     xml_str = get_default_xml()
    # print(xml_str, file=fl)
    # fl.close()


def copy_necessary_files(name, dir):
    os.symlink(f"{dir}/flags_cst", name + "/flags")
    os.symlink(
        f"{dir}/abinitio_remodel_cen_stage0a.wts",
        name + "/abinitio_remodel_cen_stage0a.wts",
    )
    os.symlink(
        f"{dir}/abinitio_remodel_cen_stage0b.wts",
        name + "/abinitio_remodel_cen_stage0b.wts",
    )
    os.symlink(
        f"{dir}/abinitio_remodel_cen_stage1.wts",
        name + "/abinitio_remodel_cen_stage1.wts",
    )
    os.symlink(
        f"{dir}/abinitio_remodel_cen_stage2.wts",
        name + "/abinitio_remodel_cen_stage2.wts",
    )
    os.symlink(
        f"{dir}/abinitio_remodel_cen.wts", name + "/abinitio_remodel_cen.wts"
    )
    os.symlink(f"{dir}/cmd", name + "/cmd")
    os.symlink(f"{dir}/start.pdb", name + "/start.pdb")


def add_flags(path_name, name, design_length):
    fl = open(path_name + "/motif_flags", "w")
    fl.write("-score:motif_residues ")
    # design_length = get_design_length(name)
    # print("design_length:{}",design_length)
    for ii in range(design_length, design_length * 2):
        fl.write("{},".format(ii))
    fl.write("{}".format(design_length * 2))
    fl.close


def get_design_length(name):
    name_wo_filename = os.path.splitext(name)[0]
    component_list = name_wo_filename.split("_")
    length = 0
    for ii in range(1, (len(component_list) - 1)):
        tmpItem = component_list[ii].replace("h", "")
        tmpItem2 = tmpItem.replace("l", "")
        if tmpItem2[0] != "t":
            length += int(tmpItem2)
    return length


def create_atom_pair_cst_string(
    resnums,
    lattice_space,
    repeat_size,
    cst_tolerance,
    sheet_mode=False,
    sheet_sd=0.2,
):
    """
    Creates atom pair csts between repeats. Sheet mode adds N-O H-O also
    """
    lines = ""
    cst_template = "AtomPair CA %s CA %s HARMONIC %s %s"

    # cst_template_NO = "AtomPair N  %s O  %s HARMONIC 2.8 0.2"
    # cst_template_HO = "AtomPair H  %s O  %s HARMONIC 1.8 0.2"
    cst_template_CB = "AtomPair CB %s CB %s HARMONIC 4.8 0.5"
    c = str("%.2f" % lattice_space).rjust(4)
    d = str("%.2f" % cst_tolerance).rjust(2)
    # ref_0 = resnums[0]
    for i in resnums:
        # if ss_list[i-1] == 'h':
        a = str(int(i)).rjust(2)
        b = str(int(repeat_size + i)).rjust(3)
        line = cst_template % (a, b, c, d)
        lines += line
        lines += "\n"
        if sheet_mode:  # and (i - ref_0 + 1) % 2:
            a = str(int(i)).rjust(2)
            b = str(int(repeat_size + i - 1)).rjust(3)
            line_cb = cst_template_CB % (a, b)
            # line_no = cst_template_NO % (a, b)
            # line_ho = cst_template_HO % (a, b)
            lines += line_cb + "\n"
            # lines += line_ho + "\n"
            # lines += line_no + "\n"
    return lines


def create_csts(name, ss_elements, lattice_space, cst_tolerance):
    ss_lengths = [int(element[1:]) for element in ss_elements]
    # for i, element in enumerate(ss_elements):
    #    for x in range(ss_lengths[i]):
    # 	  ss_list.append(ss_elements[0][0])
    rep_len = sum(ss_lengths)
    cst_template = "AtomPair CA %s CA %s HARMONIC %s %s"
    c = str("%.2f" % lattice_space).rjust(4)
    d = str("%.2f" % cst_tolerance).rjust(2)
    with open(name + "/lattice_csts.cst", "w") as fout:
        for rep_unit in [0]:
            for i in range(1, ss_lengths[0] + 1):
                # if ss_list[i-1] == 'h':
                a = str(int(i + rep_unit * rep_len)).rjust(2)
                b = str(int(rep_len + i + rep_unit * rep_len)).rjust(3)
                line = cst_template % (a, b, c, d)
                fout.write(line + "\n")
            for i in range(
                ss_lengths[0] + ss_lengths[1] + 1,
                ss_lengths[0] + ss_lengths[1] + ss_lengths[2] + 1,
            ):
                # if ss_list[i-1] == 'h':
                a = str(int(i + rep_unit * rep_len)).rjust(2)
                b = str(int(rep_len + i + rep_unit * rep_len)).rjust(3)
                line = cst_template % (a, b, c, d)
                fout.write(line + "\n")


def pose_to_blueprint(pose, offset=0, clip=2):
    """
    """
    # chain_a = pose.split_by_chain()[1]
    out = ""
    for resi in range(1 + offset, len(pose.residues) + 1 - clip):
        resn = pose.residue(resi).name1()
        out += f"{resi} {resn} ."
        out += "\n"
    out = out[:-1]
    return out


def create_blueprint(
    path, ss_elements, extra_pose=None, append=False, abego=False
):
    """
    Dumps a blueprint file at the given path

    ss_elements should be a list of tuples (dssp_type,length) to be bprint built
    """
    fl = open(path, "w")
    first = True
    extra_pose_a = None
    if extra_pose:

        extra_pose_a = extra_pose.split_by_chain()[1]
        if append:

            # This dumps all but the last two positions in the pose to the bp
            # The last res before the dump is special, as is the last for
            # bp insertions
            print(pose_to_blueprint(extra_pose_a), file=fl)
    helix_type = "H" if abego else "HA"
    sheet_type = "E" if abego else "ED"
    loop_type = "L" if abego else "LD"
    tmpType = "LD"
    for ssType, ssLength in ss_elements:
        if ssType.lower() == "h":
            tmpType = helix_type
        if ssType.lower() == "l":
            tmpType = loop_type
        if ssType.lower() == "e":
            tmpType = sheet_type
        if first:
            if extra_pose_a and append:
                pre_append_pos = len(extra_pose_a.residues) - 1
                print(
                    "{} {} {}".format(
                        pre_append_pos,
                        extra_pose_a.residue(pre_append_pos).name1(),
                        tmpType,
                    ),
                    file=fl,
                )

            else:
                print("{} {} {}".format(1, "A", tmpType), file=fl)
            first = False
            for ii in range(1, ssLength):
                print("{} {} {}".format(0, "x", tmpType), file=fl)
        else:
            for ii in range(0, ssLength):
                print("{} {} {}".format(0, "x", tmpType), file=fl)
    if extra_pose_a and append:
        last_pos = len(extra_pose_a.residues)
        print(
            "{} {} {}".format(
                last_pos, extra_pose_a.residue(last_pos).name1(), tmpType
            ),
            file=fl,
        )
    if extra_pose_a and not append:
        print(
            "{} {} {}".format(2, extra_pose_a.residue(2).name1(), tmpType),
            file=fl,
        )
        print(pose_to_blueprint(extra_pose_a, offset=2, clip=0), file=fl)

    fl.close()


def build_from_params(params,):
    ""
    sse_list = []
    for pstr in params:
        try:
            dssp_type, min_size, max_size, repeat_dist, repeat_dist_cst = pstr.split(
                " "
            )
        except ValueError:
            print("no constraints entered for: ", pstr)
            dssp_type, min_size, max_size, = pstr.split(" ")
            repeat_dist, repeat_dist_cst = [0] * 2
        sse_list.append(
            SecondaryStructElementSampler(
                dssp_type,
                int(min_size),
                int(max_size),
                float(repeat_dist),
                float(repeat_dist_cst),
            )
        )
    return sse_list


def build_from_file(file):
    ""
    return [SecondaryStructElement.from_dict(d) for d in json.load(file)]


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
    "-d/ ",
    "--disable-abego/--dssp-mode",
    "abego",
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
    abego=False,
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
    fragerator = product(*(f.get_ss_elements_list() for f in sse_sampler_list))

    dump_list = []

    for ss_elements in fragerator:
        print(ss_elements)
        elem_dict_list = [sse.to_dict() for sse in ss_elements]
        dump_list.append(elem_dict_list)

    with open("frag_params.json", "w") as f:
        json.dump(dump_list, f)

    # extra_pose = (
    #     safe_load_pdb(extra_pdb, rosetta_flags_file=rosetta_flags_file)
    #     if extra_pdb
    #     else None
    # )
    # prepare_design_dir(
    #     output_dir,
    #     elem_list,
    #     extra_files_dir,
    #     extra_pose=extra_pose,
    #     append=append,
    #     abego=abego,
    # )


if __name__ == "__main__":
    main()
