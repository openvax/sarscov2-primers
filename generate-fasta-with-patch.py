#!/usr/local/bin/python3

"""

Generate FASTA file from primer TSV and also generate version where one primer
is replaced to reduce dropout, proposed in:

    A proposal of an alternative primer for the ARTIC Network’s multiplex
    PCR to improve coverage of SARS-CoV-2 genome sequencing
    -by-
    Kentaro Itokawa*, Tsuyoshi Sekizuka, Masanori Hashino,
    Rina Tanaka, Makoto Kuroda*
"""

import sys
import pandas as pd
import os.path


def write_fasta(df, base_filename):

    with open("%s.fa" % base_filename, "w") as f:
        for name, row in df.iterrows():
            f.write(">%s pool=%s length=%d\n%s\n" % (
                name,
                row.pool,
                row.length,
                row.seq))
    pools = sorted(df.pool.unique())
    for pool in pools:
        with open("%s-pool-%s.fa" % (base_filename, pool), "w") as f:
            df_pool = df[df.pool == pool]
            for name, row in df_pool.iterrows():
                f.write(">%s length=%d\n%s\n" % (
                    name,
                    row.length,
                    row.seq))

def replace_primer(df):
    df = df.copy()
    primer_name = "nCoV-2019_76_RIGHT"
    row = df.loc[primer_name]
    old_primer = "ACACCTGTGCCTGTTAAACCAT"
    assert row.seq == old_primer
    new_primer = "TCTCTGCCAAATTGTTGGAAAGGCA"
    df.loc[primer_name] = [row.pool, new_primer, len(new_primer)]
    return df

def main():
    if len(sys.argv) < 2:
        raise ValueError("No path given for original primers")
    path = sys.argv[1]
    if not os.path.exists(path):
        raise ValueError("Primers not found at '%s'" % path)
    df = pd.read_csv(path, sep="\t").set_index("name")[["pool", "seq", "length"]]
    write_fasta(df, "artic-ncov2019-primers")

    df_itokawa = replace_primer(df)
    write_fasta(df_itokawa, "artic-ncov2019-primers-with-itokawa-patch")

if __name__ == "__main__":
    main()