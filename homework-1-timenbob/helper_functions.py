from typing import Tuple, Generator, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def codons(seq: str) -> Generator[str, None, None]:
    """Walk along the string, three nucleotides at a time. Cut off excess."""
    for i in range(0, len(seq) - 2, 3):
        yield seq[i:i + 3]


def extract_gt_orfs(record, start_codons, stop_codons, validate_cds=True, verbose=False):
    """Extract the ground truth ORFs as indicated by the NCBI annotator in the
    gene coding regions (CDS regins) of the genome.

    Parameters
    ----------
    record: SeqRecord
    start_codons: List[str]
    stop_codons: List[str]
    validate_cds: bool
        Filter out NCBI provided ORFs that do not fit our ORF criteria.
    verbose: bool

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    cds_regions = [f for f in record.features if f.type == "CDS"]

    orfs = []
    for region in cds_regions:
        loc = region.location
        #tu sem naredil spremembo , ker ne dela !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        seq = record.seq[int(loc.start):int(loc.end)]          # FIX: int(...) namesto .position
        #tu sem naredil spremembo , ker ne dela !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if region.strand == -1:
            seq = seq.reverse_complement()
            
        if not validate_cds:
            orfs.append((region.strand, int(loc.start), int(loc.end)))  # FIX: int(...) namesto .position
            continue

        try:
            assert str(seq[:3]) in start_codons, "Start codon not found!"   # FIX: str(...) za robustnost
            assert str(seq[-3:]) in stop_codons, "Stop codon not found!"    # FIX: str(...) za robustnost
            # Make sure there are no stop codons in the middle of the sequence
            for codon in codons(str(seq[3:-3])):                             # FIX: str(...) za generator kodonov
                assert (
                    codon not in stop_codons
                ), f"Stop codon {codon} found in the middle of the sequence!"

            # The CDS looks fine, add it to the ORFs
            orfs.append((region.strand, int(loc.start), int(loc.end)))       # FIX: int(...) namesto .position

        except AssertionError as ex:
            if verbose:
                print(
                    "Skipped CDS at region [%d - %d] on strand %d"
                    % (int(loc.start), int(loc.end), region.strand)          # FIX: int(...) namesto .position
                )
                print("\t", str(ex))
                
    # Some ORFs in paramecium have lenghts not divisible by 3. Remove these
    orfs = [orf for orf in orfs if (orf[2] - orf[1]) % 3 == 0]

    return orfs


def find_orfs(sequence, start_codons, stop_codons):
    """Find possible ORF candidates in a single reading frame.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int]]
        tuples of form (start_loc, stop_loc)

    """
    odg = []
    seq = str(sequence)
    i=0
    while i < len(seq)-2:
        if seq[i:i+3] in start_codons:
            zacetek = i
            konec = None
            for j in range(i + 3, len(seq) - 2, 3):
                if seq[j:j+3] in stop_codons:
                    konec = j + 3
                    odg.append((zacetek, konec))
                    i=konec
                    break
            if konec == None: i+=3
        else:
            i+=3
    return odg


def find_all_orfs(sequence, start_codons, stop_codons):
    """Find ALL the possible ORF candidates in the sequence using all six
    reading frames.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    odg = []
    seq = str(sequence)
    dolzina = len(seq)
    for i in range(3):
        orfs = find_orfs(seq[i:], start_codons, stop_codons)
        for (zacetek, konec) in orfs:
            odg.append((1, zacetek + i, konec + i))

    seq2=str(sequence.reverse_complement())
    for i in range(3):
        orfs = find_orfs(seq2[i:], start_codons, stop_codons)
        for (zacetek, konec) in orfs:
            odg.append((-1, dolzina - konec - i ,dolzina-zacetek-i))
    return odg


def translate_to_protein(seq):
    odg = ""
    # to tabelco mi je chatGPT naredil!!
    codon_table = {
        # --- Fenilalanin (F), Levcin (L) ---
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        # --- Izolevcin (I), Metionin (M/start) ---
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        # --- Valin (V) ---
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        # --- Serin (S), Prolin (P), Treonin (T), Alanin (A) ---
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        # --- Tirozin (Y), Histidin (H), Glutamin (Q), Asparagin (N), Lizin (K) ---
        "TAT": "Y", "TAC": "Y",
        "CAT": "H", "CAC": "H",
        "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N",
        "AAA": "K", "AAG": "K",
        # --- Aspartat (D), Glutamat (E) ---
        "GAT": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        # --- Cistein (C), Triptofan (W), Arginin (R), Serin (S), Glicin (G) ---
        "TGT": "C", "TGC": "C", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        # --- Stop kodoni ---
        "TAA": "", "TAG": "", "TGA": ""
    }
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in codon_table:
            aa = codon_table[codon]
            odg += aa

    return odg
    



def find_all_orfs_nested(sequence, start_codons, stop_codons):
    """Bonus problem: Find ALL the possible ORF candidates in the sequence using
    the updated definition of ORFs.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    odg = []
    seq = str(sequence)
    dolzina = len(seq)
    

    for i in range(3):
        orfs = seq[i:]
        for j in range(0,dolzina-i-2,3):
            if orfs[j:j+3] in start_codons:
                zac=i+j
                for k in range(j+3, dolzina-i-2,3):
                    if orfs[k:k+3] in stop_codons:
                        konec=i+k+3
                        odg.append((1,zac,konec))

    seq2=str(sequence.reverse_complement())
    for i in range(3):
        orfs = seq2[i:]
        for j in range(0,dolzina-i-2,3):
            if orfs[j:j+3] in start_codons:
                for k in range(j+3, dolzina-i-2,3):
                    if orfs[k:k+3] in stop_codons:
                        zac=dolzina-i-k-3
                        konec=dolzina-i-j
                        odg.append((-1,zac,konec))

            
    return odg
