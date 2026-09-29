#!/usr/bin/env python3


"""
==================================================
PolyHIC Module 3

Hi-C assisted haplotype validation

==================================================
"""


import argparse
import subprocess
import os
import pysam
import pandas as pd



##############################################
# Merge local haplotypes
##############################################


def merge_haplotype_fasta(
        fasta_files,
        output):


    with open(output,"w") as out:


        for f in fasta_files:


            with open(f) as inp:

                for line in inp:

                    out.write(line)



##############################################
# Hi-C mapping
##############################################


def hic_mapping(
        reference,
        hic1,
        hic2,
        prefix,
        threads):


    subprocess.run(
        f"bwa-mem2 index {reference}",
        shell=True
    )


    cmd=f"""

    bwa-mem2 mem \
    -t {threads} \
    {reference} \
    {hic1} \
    {hic2} \
    | samtools sort \
    -o {prefix}.bam

    """


    subprocess.run(
        cmd,
        shell=True
    )


    subprocess.run(
        f"samtools index {prefix}.bam",
        shell=True
    )



##############################################
# Calculate Hi-C contact
##############################################


def hic_contact(
        bam):


    sam=pysam.AlignmentFile(
        bam,
        "rb"
    )


    contacts={}



    for read in sam.fetch():


        if not read.is_read1:

            continue


        mate_chr=read.next_reference_id


        chr1=read.reference_id


        if chr1 < 0 or mate_chr <0:

            continue


        key=(
            chr1,
            mate_chr
        )


        contacts[key]=(
            contacts.get(key,0)+1
        )



    return contacts



##############################################
# CIS calculation
##############################################


def calculate_CIS(
        contacts,
        hap1,
        hap2):


    intra1=contacts.get(
        (hap1,hap1),
        0
    )


    intra2=contacts.get(
        (hap2,hap2),
        0
    )


    inter=contacts.get(
        (hap1,hap2),
        0
    )


    cis=(
        intra1+
        intra2
    )/(inter+1)


    return cis



##############################################
# main
##############################################


def main():


    parser=argparse.ArgumentParser()


    parser.add_argument(
        "--hap1",
        required=True
    )


    parser.add_argument(
        "--hap2",
        required=True
    )


    parser.add_argument(
        "--hic1",
        required=True
    )


    parser.add_argument(
        "--hic2",
        required=True
    )


    parser.add_argument(
        "-t",
        "--threads",
        default=48
    )


    parser.add_argument(
        "-o",
        "--output",
        default="PolyHIC_HAV"
    )


    args=parser.parse_args()



    os.makedirs(
        args.output,
        exist_ok=True
    )



    ############################

    # merge fasta

    ############################


    ref=(
        args.output+
        "/haplotypes.fa"
    )


    merge_haplotype_fasta(
        [
            args.hap1,
            args.hap2
        ],
        ref
    )



    ############################

    # HiC mapping

    ############################


    hic_mapping(
        ref,
        args.hic1,
        args.hic2,
        args.output+"/hic",
        args.threads
    )



    ############################

    # contact

    ############################


    contacts=hic_contact(
        args.output+
        "/hic.bam"
    )



    # assume:
    # hap1 chromosome id=0
    # hap2 chromosome id=1


    cis=calculate_CIS(
        contacts,
        0,
        1
    )


    print(
        "Contact Isolation Score:",
        cis
    )



    if cis>2:


        result="validated"


    else:


        result="uncertain"



    with open(
        args.output+
        "/validation.txt",
        "w"
    ) as f:


        f.write(
            result+
            "\n"
        )



if __name__=="__main__":

    main()

