def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    If multiple alignments achieve the same optimal score, you may return any
    one of them.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> def scoring_function(x, y):
    ...     return -5 if x == "-" or y == "-" else [-1, 5][x == y]
    ...
    >>> global_alignment("abracadabra", "dabarakadara", scoring_function)
    ('-ab-racadabra', 'dabarakada-ra', 29.0)

    Other alignments are not possible.

    """

    #inicalization
    M={(0,0):0}
    M.update({(i+1,0):scoring_function(c,"-")*(i+1) for i,c in enumerate(seq1)})
    M.update({(0,j+1):scoring_function("-",c)*(j+1) for j,c in enumerate(seq2)})

    P={(i,0):(i-1,0) for i in range(1,len(seq1)+1)}
    P.update({(0,j):(0,j-1) for j in range(1,len(seq2)+1)})

    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            M[i,j],P[i,j]=max(
                (M[i-1,j] + scoring_function(seq1[i-1],"-"),(i-1,j)),
                (M[i,j-1] + scoring_function("-",seq2[j-1]),(i,j-1)),
                (M[i-1,j-1] + scoring_function(seq1[i-1],seq2[j-1]),(i-1,j-1)))

    walk=[(len(seq1),len(seq2))]
    while P[walk[-1]] != (0,0):
        walk.append(P[walk[-1]])

    walk.reverse()

    si,ti=zip(*walk)
    out1,prev1=[],0
    for i in si:
        out1.append("-" if (i==prev1) else seq1[i-1])
        prev1=i
    out2,prev2=[],0
    for i in ti:
        out2.append("-" if (i==prev2) else seq2[i-1])
        prev2=i

    return("".join(out1),"".join(out2),M[(len(seq1),len(seq2))])
    

def local_alignment(seq1, seq2, scoring_function):
    M = {(0, 0): 0.0}
    for i in range(1, len(seq1) + 1):
        M[(i, 0)] = 0.0
    for j in range(1, len(seq2) + 1):
        M[(0, j)] = 0.0

    P = {}
    najbolsi_m = 0.0
    najbolsi_m_kje = (0, 0)


    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            up = M[(i - 1, j)] + scoring_function(seq1[i - 1], "-")
            left = M[(i, j - 1)] + scoring_function("-", seq2[j - 1])
            diag = M[(i - 1, j - 1)] + scoring_function(seq1[i - 1], seq2[j - 1])
            best = max(0.0, up, left, diag)


            M[(i, j)] = best
            if best == 0.0:
                P[(i, j)] = None
            elif best == diag:
                P[(i, j)] = (i - 1, j - 1)
            elif best == up:
                P[(i, j)] = (i - 1, j)
            else:
                P[(i, j)] = (i, j - 1)

            if best > najbolsi_m:
                najbolsi_m = best
                najbolsi_m_kje = (i, j)


    i, j = najbolsi_m_kje
    a1, a2 = [], []
    while (i, j) in P and P[(i, j)] is not None and M[(i, j)] > 0:
        pi, pj = P[(i, j)]
        if pi == i - 1 and pj == j - 1:
            a1.append(seq1[i - 1])
            a2.append(seq2[j - 1])
        elif pi == i - 1 and pj == j:
            a1.append(seq1[i - 1])
            a2.append("-")
        else:
            a1.append("-")
            a2.append(seq2[j - 1])
        i, j = pi, pj

    a1.reverse()
    a2.reverse()

    return "".join(a1), "".join(a2), najbolsi_m
