#!/usr/bin/env python3


import argparse
import pandas as pd
import subprocess
import os



########################################
# merge fused windows
########################################


def merge_fused_regions(
        classification,
        min_windows=3):


    df=pd.read_csv(
        classification,
        sep="\t"
    )


    fused=df[
        df.classification==
        "haplotype-fused"
    ]


    regions=[]


    for contig,group in fused.groupby(
        "contig"
    ):


        group=group.sort_values(
            "start"
        )


        start=None
        end=None
        count=0


        for _,r in group.iterrows():


            if start is None:

                start=r.start
                end=r.end
                count=1


            elif r.start <= end:


                end=r.end
                count+=1


            else:


                if count>=min_windows:

                    regions.append(
                        [
                            contig,
                            start,
                            end
                        ]
                    )


                start=r.start
                end=r.end
                count=1



        if count>=min_windows:

            regions.append(
                [
                    contig,
                    start,
                    end
                ]
            )


    return pd.DataFrame(
        regions,
        columns=[
            "contig",
            "start",
            "end"
        ]
    )



########################################
# extract HiFi reads
########################################


def extract_hifi(
        bam,
        region,
        output):


    cmd=f"""

    samtools view \
    -b {bam} \
    {region} \
    | samtools fastq \
    - > {output}

    """


    subprocess.run(
        cmd,
        shell=True
    )



########################################
# local assembly
########################################


def run_hifiasm(
        hifi,
        prefix,
        threads):


    cmd=f"""

    hifiasm \
    -o {prefix} \
    -t {threads} \
    {hifi}

    """


    subprocess.run(
        cmd,
        shell=True
    )



########################################

def main():


    parser=argparse.ArgumentParser()


    parser.add_argument(
        "-c",
        "--classification",
        required=True
    )


    parser.add_argument(
        "--hifi-bam",
        required=True
    )


    parser.add_argument(
        "--threads",
        default=48
    )


    parser.add_argument(
        "-o",
        "--output",
        default="PolyHIC_LHR"
    )


    args=parser.parse_args()



    os.makedirs(
        args.output,
        exist_ok=True
    )



    regions=merge_fused_regions(
        args.classification
    )


    regions.to_csv(
        args.output+
        "/fused_regions.bed",
        sep="\t",
        index=False,
        header=False
    )


    for i,r in regions.iterrows():


        name=f"{r.contig}_{r.start}_{r.end}"


        region=(
            f"{r.contig}:"
            f"{r.start}-"
            f"{r.end}"
        )


        hifi_out=(
            args.output+
            "/"+
            name+
            ".hifi.fastq"
        )


        extract_hifi(
            args.hifi_bam,
            region,
            hifi_out
        )


        run_hifiasm(
            hifi_out,
            args.output+
            "/"+
            name+
            ".local",
            args.threads
        )



if __name__=="__main__":

    main()

