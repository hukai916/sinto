import sys
from contextlib import closing

import pysam
from sinto.constants import OUT_FORMAT_CONVERSION


def header_line_to_str(line):
    return "\t".join(f"{k}:{v}" for k, v in line.items())


def build_header(header, tag_vals):
    new_rg_lines = []
    rg = {}
    for val in tag_vals:
        rg["ID"] = f"{val}"
        rg["SM"] = f"{val}"
        rg["LB"] = f"{val}"
        new_rg_lines.append("@RG\t" + header_line_to_str(rg))
    return str(header) + "\n".join(new_rg_lines) + "\n"


def tagtorg(bam, tag, output, out_format="t"):
    """Add tags to reads from individual cells

    Copies BAM entries to a new file, adding a read tag to cells matching an input table

    Parameters
    ----------
    bam : str
        Path to SAM/BAM file, or "-" to read from stdin.
    output : str
        Name for output SAM/BAM file, or "-" to write to stdout.
    out_format : str
        One of "t" (SAM), "b" (BAM), or "u" (uncompressed BAM) ("t" default)
    """

    infile   = pysam.AlignmentFile(bam)
    infile2  = pysam.AlignmentFile(bam)
    tag_vals = set()
    for rec in infile:
        try:
            tag_val = rec.get_tag(tag)
        except KeyError:
            pass
        else:
            tag_vals.add(tag_val)
    assert len(tag_vals) > 0, "supplied --tag \"" + tag + "\" not found in supplied bam file!"

    new_header = build_header(infile.header, tag_vals)

    outfile = pysam.AlignmentFile(
        output, "w" + OUT_FORMAT_CONVERSION[out_format], text=new_header
    )

    with closing(outfile) as outfile:
        for rec in infile2:
            try:
                tag_val = rec.get_tag(tag)
            except KeyError:
                pass
            else:
                rec.set_tag(
                    "RG",
                    tag_val,
                    value_type="Z",
                    replace=True)
            finally:
                outfile.write(rec)
