import numpy as np
def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    ok_pari = [(a, b) for a, b in zip(reference_sequence, distant_sequence) if a != '-' and b != '-']
    
    #observed diff to je haming
    pb=0
    dolzina=0
    for a, b in ok_pari:
        if a != b : pb += 1
        dolzina += 1

    pt = pb / dolzina
    
    djc = -3/4 * np.log(1 - (4/3) * pt)*dolzina

    return djc


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    
    def transition(n1, n2):
        '''Spisu llm'''
        return (n1 == "A" and n2 == "G") or (n1 == "G" and n2 == "A") or (n1 == "C" and n2 == "T") or (n1 == "T" and n2 == "C")
    
    def transversion(n1, n2):
        '''spisu llm'''
        return (n1 == "A" and n2 == "C") or (n1 == "C" and n2 == "A") or (n1 == "A" and n2 == "T") or (n1 == "T" and n2 == "A") or (n1 == "G" and n2 == "C") or (n1 == "C" and n2 == "G") or (n1 == "G" and n2 == "T") or (n1 == "T" and n2 == "G")

    
    SUMTransitions = 0
    SUMTransversions = 0 
    dolzina = 0

    ok_pari = [(a, b) for a, b in zip(reference_sequence, distant_sequence) if a != '-' and b != '-']

    for a, b in ok_pari:
        if a != b : 
            if transition(a,b):
                SUMTransitions+=1
            elif transversion(a,b): 
                SUMTransversions+=1
        dolzina += 1
    

    P = SUMTransitions / dolzina
    Q = SUMTransversions / dolzina
    odg = (-1/2 * np.log((1 - 2 * P - Q))- 1/4 * np.log(1 - 2 * Q) ) * dolzina
    return odg
