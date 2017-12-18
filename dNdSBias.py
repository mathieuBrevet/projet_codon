import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt


def nullspace(A, atol=1e-13, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = np.linalg.svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns


def equilibrium_frequencies(matrix):
    null_vec = nullspace(matrix.T)
    frequencies = null_vec.T[0]
    frequencies /= np.sum(frequencies)
    assert np.sum(frequencies) - 1 < 1e-13, "The frequencies does not sum to 1"
    assert np.sum(np.dot(frequencies, matrix)) < 1e-13, "The frequencies are not the nullspace of the transposed matrix"
    assert np.sum(np.dot(frequencies, expm(matrix)) - frequencies) < 1e-10, "The frequencies are not stationary"
    return frequencies


def mutation_between_codon(_codon_origin, _codon_target):
    assert len(_codon_origin) == len(_codon_target)
    return [(nuc_origin, nuc_target) for nuc_origin, nuc_target in zip(_codon_origin, _codon_target) if
            nuc_origin != nuc_target]


nucleotides = "ACGT"
codontable = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}
assert len(codontable.keys()) == 64, "There is not 64 codons in the codon table"

codons = [c for c, aa in codontable.items() if aa != 'X']
assert len(codons) == 61, "There is not 3 stop codons in the codon table"
assert not [n for n in "".join(codons) if (n not in nucleotides)], "There is a codon with an unrecognized nucleotide"

amino_acids_set = set(codontable.values())
amino_acids_set.remove('X')
amino_acids = "".join(amino_acids_set)
assert len(amino_acids) == 20, "There is not 20 amino-acids in the codon table"


aa_table = {}
for codon in codons:
    aa = codontable[codon]
    if aa not in aa_table:
        aa_table[aa] = []
    aa_table[aa].append(codons.index(codon))
assert len(aa_table) == 20, "There is not 20 amino-acids in the aa table"


codon_preference = True
if codon_preference:
    codon_preference = np.random.uniform(1, 2)

sigma_2 = np.random.uniform(1, 4)
assert 0 <= sigma_2 <= 4
sigma = np.sqrt(sigma_2)
fitness_aa = np.random.normal(loc=0, scale=sigma, size=len(amino_acids))
fitness_codons = np.zeros(len(codons))

for aa, codon_list in aa_table.items():
    prefered = np.random.choice(codon_list)
    fitness_codon = fitness_aa[amino_acids.find(aa)]
    fitness_codons[prefered] = fitness_codon + codon_preference
    for codon in codon_list:
        if codon != prefered:
            fitness_codons[codon] = fitness_codon - codon_preference

points = 50

at_percent = np.zeros(points)
at_obs_percent = np.zeros(points)
gamma = np.zeros(points)
gamma_n = np.zeros(points)
gamma_s = np.zeros(points)
mut_bias_range = np.logspace(-1, 1, points)
for mut_bias_index, mut_bias in enumerate(mut_bias_range):
    mutation_matrix = np.array([[0, 1, 1, mut_bias],
                                [mut_bias, 0, 1, mut_bias],
                                [mut_bias, 1, 0, mut_bias],
                                [mut_bias, 1, 1, 0]])
    mutation_matrix -= np.diagflat(np.sum(mutation_matrix, axis=1))
    assert np.sum(mutation_matrix) < 1e-10, "The mutation matrix don't have the rows summing to 0"

    nuc_frequencies = equilibrium_frequencies(mutation_matrix)
    at_percent[mut_bias_index] = nuc_frequencies[nucleotides.find("A")] + nuc_frequencies[nucleotides.find("T")]

    codon_matrix = np.zeros((len(codons), len(codons)))
    for c_origin_index, codon_origin in enumerate(codons):
        for c_target_index, codon_target in enumerate(codons):
            mutations = mutation_between_codon(codon_origin, codon_target)
            if len(mutations) == 1:
                nuc_origin, nuc_target = mutations[0]
                n_origin_index = nucleotides.find(nuc_origin)
                n_target_index = nucleotides.find(nuc_target)
                sel_coef = fitness_codons[c_target_index] - fitness_codons[c_origin_index]
                if abs(sel_coef) < 1e-10:
                    p_fix = 1
                else:
                    p_fix = sel_coef / (1. - np.exp(-sel_coef))
                codon_matrix[c_origin_index][c_target_index] = mutation_matrix[n_origin_index][n_target_index] * p_fix

    codon_matrix -= np.diagflat(np.sum(codon_matrix, axis=1))
    assert np.sum(codon_matrix) < 1e-10, "The codon matrix don't have the rows summing to 0"

    codon_frequencies = equilibrium_frequencies(codon_matrix)

    for codon_index, freq in enumerate(codon_frequencies):
        at_obs_percent[mut_bias_index] += sum([freq/3 for nuc in codons[codon_index] if nuc == "A" or nuc == "T"])

    dn, dn0, ds, ds0 = 0, 0, 0, 0
    for c_origin_index, codon_origin in enumerate(codons):
        for c_target_index, codon_target in enumerate(codons):
            mutations = mutation_between_codon(codon_origin, codon_target)
            if len(mutations) == 1:
                nuc_origin, nuc_target = mutations[0]
                n_origin_index = nucleotides.find(nuc_origin)
                n_target_index = nucleotides.find(nuc_target)
                d_partiel = codon_frequencies[c_origin_index] * codon_matrix[c_origin_index][c_target_index]
                d0_partiel = codon_frequencies[c_origin_index] * mutation_matrix[n_origin_index][n_target_index]
                if codontable[codon_origin] != codontable[codon_target]:
                    dn += d_partiel
                    dn0 += d0_partiel
                else:
                    ds += d_partiel
                    ds0 += d0_partiel
    d = dn + ds
    d0 = dn0 + ds0

    assert d != 0
    assert d0 != 0
    gamma[mut_bias_index] = d / d0
    gamma_n[mut_bias_index] = dn / dn0
    gamma_s[mut_bias_index] = ds / ds0


plt.subplot(211)
plt.plot(mut_bias_range, gamma, label='$\gamma$')
plt.plot(mut_bias_range, gamma_n, label='$\gamma_N$')
plt.plot(mut_bias_range, gamma_s, label='$\gamma_S$')
plt.plot(mut_bias_range, gamma_n / gamma_s, label='$\gamma_N / \gamma_S$')
plt.ylim((0, max(1, max(gamma_n / gamma_s))))
plt.xscale('log')
plt.ylabel('$\gamma$')
plt.title('Impact of mutational bias on $\gamma$ and %AT')
plt.legend()

plt.subplot(212)
plt.plot(mut_bias_range, at_percent, label='%AT pred')
plt.plot(mut_bias_range, at_obs_percent, label='%AT obs')
plt.plot(mut_bias_range, 1 - at_percent, label='%GC pred')
plt.plot(mut_bias_range, 1 - at_obs_percent, label='%GC obs')
plt.xscale('log')
plt.xlabel('$\lambda$')
plt.ylabel('% nucleotides')
plt.legend()
plt.show()
