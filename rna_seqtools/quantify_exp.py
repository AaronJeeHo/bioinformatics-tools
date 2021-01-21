"""
Author: Aaron Ho
Python Version: 3.8
"""


def RPK(read_count, gene_length):
    """
    Computes RPK (Reads Per Kilobase) for a given gene
    :param int read_count: number of reads mapped to current gene
    :param int gene_length: length of the current gene
    :return float: RPK value
    """
    return read_count / (gene_length / 1000)


def RPKM(read_count, gene_length, total_reads):
    """
    Computes RPKM (Reads Per Kilobase per Million reads) for a given gene
    :param int read_count: number of reads mapped to current gene
    :param int gene_length: length of the current gene
    :param int total_reads: total reads sequenced
    :return float: RPKM value
    """
    return read_count / ((gene_length / 1000) * (total_reads / 1e6))


def TPM(gene_rpk, scale_factor):
    """
    Computes TPM (Transcripts Per Million) for a given gene
    :param float gene_rpk: RPK value of current gene
    :param int scale_factor: total RPK values divided by one million
    :return float: TPM value
    """

    return gene_rpk / scale_factor


def get_expression(gene_lengths, read_counts):
    """
    Quantifies normalized gene expressions in TPM
    :param list gene_lengths: int list of gene lengths per gene
                                eg. [gene1-length, gene2-length, gene3-length]
    :param list read_counts: int list of read counts per gene
                                eg. [gene1-count, gene2-count, gene3-count]
    :return list: float list of TPM's per gene
                                eg. [gene1-TPM, gene2-TPM, gene3-TPM]
    """

    rpk_list = [read_counts[gene] / (gene_lengths[gene] / 1000)
                for gene in range(len(gene_lengths))]

    scale_factor = sum(rpk_list) / 1e6

    return [TPM(rpk_list[gene], scale_factor)
            for gene in range(len(gene_lengths))]


def main():
    pass


if __name__ == '__main__':
    main()