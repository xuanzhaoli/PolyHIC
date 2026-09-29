#!/usr/bin/env python3

"""
=========================================================
PolyHIC Module 1

Adaptive depth-based haplotype resolution classification

Input:
    1. Initial contig assembly fasta
    2. Illumina paired-end reads
    3. Genome ploidy

Output:
    1. Window-based depth file
    2. Estimated haplotype depth
    3. Window classification result


Author:
    PolyHIC development team

=========================================================
"""


import argparse
import subprocess
import pysam
import pandas as pd
import numpy as np

from sklearn.mixture import GaussianMixture



#########################################################
# Step 1
# Mapping Illumina reads
#########################################################


def mapping_reads(
        assembly,
        read1,
        read2,
        prefix,
        threads):


    print("[INFO] Indexing assembly...")


    subprocess.run(
        f"bwa-mem2 index {assembly}",
        shell=True,
        check=True
    )


    print("[INFO] Mapping reads...")


    cmd=f"""
    bwa-mem2 mem \
    -t {threads} \
    {assembly} \
    {read1} \
    {read2} \
    | samtools sort \
    -@ {threads} \
    -o {prefix}.bam
    """


    subprocess.run(
        cmd,
        shell=True,
        check=True
    )


    subprocess.run(
        f"samtools index {prefix}.bam",
        shell=True,
        check=True
    )


    print("[INFO] Mapping finished")



#########################################################
# Step 2
# Window-based depth calculation
#########################################################


def calculate_window_depth(
        bam,
        assembly,
        window_size):


    print("[INFO] Calculating window depth...")


    fasta=pysam.FastaFile(
        assembly
    )


    sam=pysam.AlignmentFile(
        bam,
        "rb"
    )


    results=[]


    for contig in fasta.references:


        length=fasta.get_reference_length(
            contig
        )


        print(
            "[INFO]",
            contig,
            length
        )


        for start in range(
                0,
                length,
                window_size
        ):


            end=min(
                start+window_size,
                length
            )


            depth_sum=0
            base_number=0



            for pileup in sam.pileup(
                    contig,
                    start,
                    end,
                    truncate=True):


                depth_sum += pileup.nsegments

                base_number += 1



            if base_number > 0:

                mean_depth = (
                    depth_sum /
                    base_number
                )

            else:

                mean_depth = 0



            results.append(
                [
                    contig,
                    start+1,
                    end,
                    mean_depth
                ]
            )



    sam.close()

    fasta.close()



    return pd.DataFrame(
        results,
        columns=[
            "contig",
            "start",
            "end",
            "depth"
        ]
    )



#########################################################
# Step 3
# Estimate haplotype depth using GMM
#########################################################


def estimate_haplotype_depth(
        depths,
        ploidy):


    print(
        "[INFO] Estimating haplotype depth..."
    )


    X=np.array(
        depths
    ).reshape(
        -1,
        1
    )


    #
    # number of possible depth peaks
    #
    # e.g.
    # tetraploid:
    # 1x 2x 3x 4x
    #

    n_components=min(
        ploidy,
        5
    )


    model=GaussianMixture(
        n_components=n_components,
        random_state=123
    )


    model.fit(X)



    peaks=sorted(
        model.means_.flatten()
    )


    hap_depth=peaks[0]



    return hap_depth, peaks



#########################################################
# Step 4
# Normalize depth and classify
#########################################################


def classify_window(
        depth,
        hap_depth):


    normalized_depth = (
        depth /
        hap_depth
    )


    if normalized_depth < 1.5:


        classification=(
            "haplotype-resolved"
        )


    elif normalized_depth < 2.5:


        classification=(
            "haplotype-fused"
        )


    else:


        classification=(
            "multi-copy"
        )


    return (
        normalized_depth,
        classification
    )



#########################################################
# Main
#########################################################


def main():


    parser=argparse.ArgumentParser(
        description=
        "PolyHIC Module 1: "
        "Depth-based haplotype classification"
    )


    parser.add_argument(
        "-a",
        "--assembly",
        required=True,
        help="initial contig fasta"
    )


    parser.add_argument(
        "-1",
        "--read1",
        required=True,
        help="Illumina read1"
    )


    parser.add_argument(
        "-2",
        "--read2",
        required=True,
        help="Illumina read2"
    )


    parser.add_argument(
        "--ploidy",
        type=int,
        required=True,
        help="genome ploidy"
    )


    parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=1000,
        help=
        "window size for depth calculation "
        "(default: 1000 bp)"
    )


    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=24
    )


    parser.add_argument(
        "-o",
        "--output",
        default="PolyHIC"
    )


    args=parser.parse_args()



    prefix=args.output



    #################################################
    # 1. Mapping
    #################################################

    mapping_reads(
        args.assembly,
        args.read1,
        args.read2,
        prefix,
        args.threads
    )



    #################################################
    # 2. Window depth
    #################################################

    depth_df=calculate_window_depth(
        prefix+".bam",
        args.assembly,
        args.window
    )


    depth_file=(
        prefix+
        ".depth.window.tsv"
    )


    depth_df.to_csv(
        depth_file,
        sep="\t",
        index=False
    )


    print(
        "[INFO] Depth file:",
        depth_file
    )



    #################################################
    # 3. Estimate haplotype depth
    #################################################

    hap_depth, peaks = (
        estimate_haplotype_depth(
            depth_df["depth"],
            args.ploidy
        )
    )



    print(
        "\nEstimated haplotype depth:"
    )

    print(
        hap_depth
    )


    print(
        "\nDepth peaks:"
    )

    print(
        peaks
    )



    with open(
        prefix+".haplotype_depth.txt",
        "w"
    ) as f:


        f.write(
            "haplotype_depth\t{}\n"
            .format(
                hap_depth
            )
        )


        f.write(
            "peaks\t{}\n"
            .format(
                ",".join(
                    map(
                        str,
                        peaks
                    )
                )
            )
        )



    #################################################
    # 4. Classification
    #################################################

    results=[]


    for _,row in depth_df.iterrows():


        nd,cls=classify_window(
            row["depth"],
            hap_depth
        )


        results.append(
            [
                row["contig"],
                row["start"],
                row["end"],
                row["depth"],
                nd,
                cls
            ]
        )



    result_df=pd.DataFrame(
        results,
        columns=[
            "contig",
            "start",
            "end",
            "depth",
            "normalized_depth",
            "classification"
        ]
    )


    outfile=(
        prefix+
        ".classification.tsv"
    )


    result_df.to_csv(
        outfile,
        sep="\t",
        index=False
    )


    print(
        "[INFO] Classification result:",
        outfile
    )


    print(
        "[INFO] PolyHIC Module 1 finished"
    )



if __name__=="__main__":

    main()
