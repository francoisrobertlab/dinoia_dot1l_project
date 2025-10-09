import argparse
import sys
import re
from typing import TextIO


class Gene:
    def __init__(self, ensembl: str, mgi: str):
        self.ensembl = ensembl
        self.mgi = mgi

class GoGene:
    def __init__(self, mgi: str):
        self.mgi = mgi


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(
        description="Extract lines associated with GO genes.")
    parser.add_argument('count_matrix', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Count matrix")
    parser.add_argument('--go', type=argparse.FileType('r'),
                        help="GO genes")
    parser.add_argument('--genes', type=argparse.FileType('r'),
                        help="All genes in GTF or GFF3 format")

    args = parser.parse_args(argv)

    genes = parse_gtf(args.genes)
    genes_per_id = {gene.ensembl: gene for gene in genes}

    go_genes = parse_go(args.go)
    go_genes_id = set([gene.mgi for gene in go_genes])

    print(args.count_matrix.readline().rstrip('\r\n'))
    for line in args.count_matrix:
        line = line.rstrip('\r\n')
        columns = line.split('\t')
        if columns[0] in genes_per_id:
            gene = genes_per_id[columns[0]]
            if gene.mgi in go_genes_id:
                print(f"{line}")


def parse_go(go: TextIO = None) -> list[GoGene]:
    genes = []
    mgi_id_pattern = re.compile(r"MGI:([\w\.:]+)")
    for line in go:
        mgi_id_match = mgi_id_pattern.search(line)
        if mgi_id_match:
            genes.append(GoGene(mgi_id_match.group(1)))
    return genes


def parse_gtf(gtf: TextIO = None) -> list[Gene]:
    genes = []
    gene_id_pattern = re.compile(r"gene_id=([\w\.]+)")
    mgi_id_pattern = re.compile(r"mgi_id=([\w\.:]+)")
    for line in gtf:
        gene_id_match = gene_id_pattern.search(line)
        mgi_id_match = mgi_id_pattern.search(line)
        if gene_id_match and mgi_id_match:
            genes.append(Gene(gene_id_match.group(1), mgi_id_match.group(1)))
    return genes


if __name__ == '__main__':
    main()
