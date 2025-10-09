import argparse
import sys
import re
import os
import itertools
from typing import Literal


class ReadsPerGene:
    def __init__(self, gene: str, reads: int, reads1: int, reads2: int):
        self.gene = gene
        self.reads = reads
        self.reads1 = reads1
        self.reads2 = reads2


def file_path(string: str):
    if not string or os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)


SPECIAL_GENES = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(
        description="Merges multiple 'ReadsPerGene.out.tab' outputs from STAR into a single count file.")
    parser.add_argument('reads', nargs='+', type=file_path,
                        help="'ReadsPerGene.out.tab' outputs from STAR")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Tab delimited output file containing merged counts")
    parser.add_argument('--column', choices=["unstranded", "reads1", "reads2"], default="unstranded",
                        help="Column containing counts to keep  (default: %(default)s)")
    parser.add_argument('--name', default=r"(.*)\.ReadsPerGene\.out\.tab",
                        help="Regex to obtain sample name from filename  (default: %(default)s)")
    parser.add_argument('--special', action=argparse.BooleanOptionalAction, default=False,
                        help="Include special gene like 'N_unmapped'  (default: %(default)s)")

    args = parser.parse_args(argv)

    samples = []
    sample_files = {}
    for reads_file in args.reads:
        sample_name_search = re.search(args.name, reads_file)
        if not sample_name_search:
            print(f"regex {args.name} not found in filename {reads_file}, sample name assumed to be filename",
                  file=sys.stderr)
        sample = sample_name_search.group(1) if sample_name_search else reads_file
        samples.append(sample)
        sample_files[sample] = reads_file

    counts_per_sample = {}
    for sample, reads_file in sample_files.items():
        reads = parse_reads(reads_file)
        counts_per_sample[sample] = {read.gene: read for read in reads}

    genes = list(set(list(itertools.chain.from_iterable(counts.keys() for counts in counts_per_sample.values()))))
    genes.sort()
    # Force 'special genes' to be at the top.
    [genes.remove(special_gene) for special_gene in SPECIAL_GENES]
    if args.special:
        genes = SPECIAL_GENES + genes

    args.output.write("Gene\t")
    args.output.write("\t".join(samples))
    args.output.write("\n")
    for gene in genes:
        args.output.write(f"{gene}")
        for sample in samples:
            sample_counts = counts_per_sample[sample]
            count = get_count(sample_counts[gene], args.column) if gene in sample_counts else 0
            args.output.write(f"\t{count}")
        args.output.write("\n")


def get_count(reads: ReadsPerGene, choice: Literal["unstranded", "reads1", "reads2"]) -> int:
    if choice == "reads1":
        return reads.reads1
    elif choice == "reads2":
        return reads.reads2
    else:
        return reads.reads


def parse_reads(reads: str) -> list[ReadsPerGene]:
    genes = []
    with open(reads, 'r') as reads_in:
        for line in reads_in:
            columns = line.rstrip('\r\n').split('\t')
            gene = columns[0]
            reads = int(columns[1])
            reads1 = int(columns[2])
            reads2 = int(columns[3])
            genes.append(ReadsPerGene(gene, reads, reads1, reads2))
    return genes


if __name__ == '__main__':
    main()
