#!/usr/bin/env python3


"""
PolyHIC Module 4

Polyploid Haplotype-aware
Chromosome Reconstruction

"""


import argparse
import pandas as pd
import numpy as np

from scipy.cluster.hierarchy import (
    linkage,
    fcluster
)

from sklearn.preprocessing import normalize



########################################
# Build similarity matrix
########################################


def build_similarity(
        contact_file):


    df=pd.read_csv(
        contact_file,
        sep="\t"
    )


    contigs=sorted(
        list(
            set(df.contig1)
            |
            set(df.contig2)
        )
    )


    matrix=np.zeros(
        (
            len(contigs),
            len(contigs)
        )
    )


    index={
        c:i
        for i,c in enumerate(contigs)
    }



    for _,row in df.iterrows():

        i=index[row.contig1]

        j=index[row.contig2]


        matrix[i,j]=row.contact

        matrix[j,i]=row.contact



    return (
        contigs,
        matrix
    )



########################################
# Hi-C clustering
########################################


def cluster_contigs(
        matrix,
        n_chr):


    #
    # normalize contact profile
    #

    profile=normalize(
        matrix,
        axis=1
    )


    #
    # distance
    #

    distance=1-np.corrcoef(
        profile
    )


    distance[np.isnan(distance)]=1



    condensed=[]


    for i in range(
        len(distance)
    ):

        for j in range(
            i+1,
            len(distance)
        ):

            condensed.append(
                distance[i,j]
            )



    Z=linkage(
        condensed,
        method="average"
    )


    labels=fcluster(
        Z,
        n_chr,
        criterion="maxclust"
    )


    return labels



########################################
# main
########################################


def main():


    parser=argparse.ArgumentParser()


    parser.add_argument(
        "-c",
        "--contact",
        required=True
    )


    parser.add_argument(
        "--chromosome-number",
        type=int,
        required=True
    )


    parser.add_argument(
        "-o",
        "--output",
        default="chromosome_groups.tsv"
    )


    args=parser.parse_args()



    contigs,matrix=build_similarity(
        args.contact
    )



    labels=cluster_contigs(
        matrix,
        args.chromosome_number
    )



    result=pd.DataFrame(
        {
            "contig":contigs,
            "chromosome_group":labels
        }
    )


    result.to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__=="__main__":

    main()
