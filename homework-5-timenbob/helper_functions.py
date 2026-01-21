from collections import namedtuple
from typing import Dict, List, Tuple
import pandas as pd
from Bio import pairwise2
import numpy as np
from scipy.stats import hypergeom

GffEntry = namedtuple(
    "GffEntry",
    [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ],
)


GeneDict = Dict[str, GffEntry]


def read_gff(fname: str) -> Dict[str, GffEntry]:
    gene_dict = {}

    with open(fname) as f:
        for line in f:
            if line.startswith("#"):  # Comments start with '#' character
                continue

            parts = line.split("\t")
            parts = [p.strip() for p in parts]

            # Convert start and stop to ints
            start_idx = GffEntry._fields.index("start")
            parts[start_idx] = int(parts[start_idx]) - 1  # GFFs count from 1..
            stop_idx = GffEntry._fields.index("end")
            parts[stop_idx] = int(parts[stop_idx]) - 1  # GFFs count from 1..

            # Split the attributes
            attr_index = GffEntry._fields.index("attribute")
            attributes = {}
            for attr in parts[attr_index].split(";"):
                attr = attr.strip()
                k, v = attr.split("=")
                attributes[k] = v
            parts[attr_index] = attributes

            entry = GffEntry(*parts)

            gene_dict[entry.attribute["gene_name"]] = entry

    return gene_dict


def split_read(read: str) -> Tuple[str, str]:
    """Split a given read into its barcode and DNA sequence. The reads are
    already in DNA format, so no additional work will have to be done. This
    function needs only to take the read, and split it into the cell barcode,
    the primer, and the DNA sequence. The primer is not important, so we discard
    that.

    The first 12 bases correspond to the cell barcode.
    The next 24 bases corresond to the oligo-dT primer. (discard this)
    The reamining bases corresond to the actual DNA of interest.

    Parameters
    ----------
    read: str

    Returns
    -------
    str: cell_barcode
    str: mRNA sequence

    """
    koda=read[:12]
    seq=read[36:]
    return koda,seq
    #raise NotImplementedError()


def map_read_to_gene(read: str, ref_seq: str, genes: GeneDict) -> Tuple[str, float]:
    """
    Align a DNA sequence (read) against a reference sequence, and map it to the best matching gene.

    This function takes a DNA sequence, aligns it to a given reference sequence,
    and determines the best matching gene from a provided gene dictionary.
    It uses local alignment and computes the Hamming distance to evaluate the similarity of the read to gene sequences.
    The function returns the name of the best matching gene and the similarity score.
    You can use 'pairwise2.align.localxs(ref_seq, read, -1, -1)' to perform local alignment.

    Parameters
    ----------
    read: str
        The DNA sequence to be aligned. This sequence should not include the cell barcode or the oligo-dT primer. It represents the mRNA fragment obtained from sequencing.
    ref_seq: str
        The complete reference sequence (e.g., a viral genome) against which the read will be aligned. This sequence acts as a basis for comparison.
    genes: GeneDict
        A dictionary where keys are gene names and values are objects or tuples containing gene start and end positions in the reference sequence. This dictionary is used to identify specific genes within the reference sequence.

    Returns
    -------
    Tuple[str, float]
        - gene: str
            The name of the gene to which the read maps best. If the read aligns best to a non-gene region, return `None`.
        - similarity: float
            The similarity score between the read and the best matching gene sequence, calculated as the Hamming distance. This is a measure of how closely the read matches a gene, with higher values indicating better matches.

    Notes
    -----
    - The function performs local alignment of the read against the reference sequence and each gene segment.
    - If the read aligns better to a region outside of any gene, the function should return `None` for the gene name.
    - The function should handle cases where no alignment is found.
    """

    local_alignment = pairwise2.align.localxs(ref_seq, read, -1, -1)# vrne seq1, seq2, score, start, stop...vrne urejenoBio.Align.PairwiseAligner
    if not local_alignment:
        return None, 1.0
    
    best_alignment = local_alignment[0]
    aligned_ref_seq = best_alignment.seqA
    #print(aligned_ref_seq)
    aligned_read_seq = best_alignment.seqB  
    #print(aligned_read_seq)
    score= best_alignment.score
    start= best_alignment.start
    end= best_alignment.end

    matches = 0
    total_compared = 0
    dolzina=0
    
    for ref_base, read_base in zip(aligned_ref_seq, aligned_read_seq):
        if read_base != '-':
            dolzina+=1
        if ref_base == '-' or read_base == '-':
            continue
        total_compared += 1
        if ref_base == read_base:
            matches += 1
    
    similarity = matches / dolzina

    
    best_gene = None
    for gene_name, gene_info in genes.items():
        gene_start = gene_info.start
        gene_end = gene_info.end 
        
        if not (start < gene_start or end > gene_end):
            best_gene = gene_name
            break
    
    return best_gene, similarity





def generate_count_matrix(
    reads: List[str], ref_seq: str, genes: GeneDict, similarity_threshold: float
) -> pd.DataFrame:
    """

    Parameters
    ----------
    reads: List[str]
        The list of all reads that will be aligned.
    ref_seq: str
        The reference sequence that the read should be aligned against.
    genes: GeneDict
    similarity_threshold: float

    Returns
    -------
    count_table: pd.DataFrame
        The count table should be an N x G matrix where N is the number of
        unique cell barcodes in the reads and G is the number of genes in
        `genes`. The dataframe columns should be to a list of strings
        corrsponding to genes and the dataframe index should be a list of
        strings corresponding to cell barcodes. Each cell in the matrix should
        indicate the number of times a read mapped to a gene in that particular
        cell.

    """
    count_dict = {}
    for read in reads:
        cell_barcode, sequence = split_read(read)#raskosam ga 

        gene, similarity = map_read_to_gene(sequence, ref_seq, genes)#dobimo kermu genu i

        if similarity >= similarity_threshold and gene is not None:
            if cell_barcode not in count_dict:
                count_dict[cell_barcode] = {}
            if gene not in count_dict[cell_barcode]:
                count_dict[cell_barcode][gene] = 0
            count_dict[cell_barcode][gene] += 1

    count_table = pd.DataFrame.from_dict(count_dict, orient='index').fillna(0).astype(int)
    count_table = count_table.reindex(columns=genes.keys(), fill_value=0)   
    return count_table



def filter_matrix(
    count_matrix: pd.DataFrame,
    min_counts_per_cell: float,
    min_counts_per_gene: float,
) -> pd.DataFrame:
    """Filter a matrix by cell counts and gene counts.
    The cell count is the total number of molecules sequenced for a particular
    cell. The gene count is the total number of molecules sequenced that
    correspond to a particular gene.
    Filtering statistics should be computed on
    the original matrix. E.g. if you filter out the genes first, the filtered
    gene molecules should still count towards the cell counts.

    Parameters
    ----------
    count_matrix: pd.DataFrame
    min_counts_per_cell: float
    min_counts_per_gene: float

    Returns
    -------
    filtered_count_matrix: pd.DataFrame

    """
    cell_counts = count_matrix.sum(axis=1)
    gene_counts = count_matrix.sum(axis=0)

    # decide which cells and genes to keep
    keep_cells = cell_counts >= min_counts_per_cell
    keep_genes = gene_counts >= min_counts_per_gene

    # filter once, using masks
    filtered_count_matrix = count_matrix.loc[keep_cells, keep_genes]

    return filtered_count_matrix
    #raise NotImplementedError()


def normalize_expressions(expression_data: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize expressions by applying natural log-transformation with pseudo count 1,
    and scaling expressions of each sample to sum up to 10000.

    Parameters
    ----------
    expression_data: pd.DataFrame
        Expression matrix with cells as rows and genes as columns.

    Returns
    -------
    normalized_data: pd.DataFrame
        Normalized expression matrix with cells as rows and genes as columns.
        Matrix should have the same shape as the input matrix.
        Matrix should have the same index and column labels as the input matrix.
        Order of rows and columns should remain the same.
        Values in the matrix should be positive or zero.
    """
    log_data = np.log(expression_data + 1)
    vrstice=log_data.sum(axis=1)
    normirano=log_data.div(vrstice, axis=0)*10000
    return normirano
    
    #raise NotImplementedError()


def hypergeometric_pval(N: int, n: int, K: int, k: int) -> float:
    """
    Calculate the p-value using the following hypergeometric distribution.

    Parameters
    ----------
    N: int
        Total number of genes in the study (gene expression matrix)
    n: int
        Number of genes in your proposed gene set (e.g. from differential expression)
    K: int
        Number of genes in an annotated gene set (e.g. GO gene set)
    k: int
        Number of genes in both annotated and proposed geneset

    Returns
    -------
    p_value: float
        p-value from hypergeometric distribution of finding such or
        more extreme match at random
    """
    return hypergeom.sf(k - 1, N, K, n)
