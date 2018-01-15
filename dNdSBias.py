import numpy as np
import matplotlib.pyplot as plt

RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
GREEN = "#8FB03E"

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

codons = list(codontable.keys())
assert len(codons) == 64, "There is not 3 stop codons in the codon table"
assert not [n for n in "".join(codons) if (n not in nucleotides)], "There is a codon with an unrecognized nucleotide"

nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
z_ww = nbr_weak
z_ws = 2 * nbr_weak
z_sw = 2 * (3 - nbr_weak)
z_ss = 3 - nbr_weak

amino_acids_set = set(codontable.values())
amino_acids_set.remove('X')
amino_acids = "".join(amino_acids_set) + "X"
assert len(amino_acids) == 21, "There is not 21 amino-acids in the codon table"

aa_table = {}
for codon in codons:
    aa = codontable[codon]
    if aa not in aa_table:
        aa_table[aa] = []
    aa_table[aa].append(codons.index(codon))
assert len(aa_table) == 21, "There is not 21 amino-acids in the aa table"

fitness_codons = nbr_weak*nbr_weak

points = 50
mut_bias_obs = np.zeros(points)
gamma_ss = np.zeros(points)
gamma_sw = np.zeros(points)
gamma_ws = np.zeros(points)
gamma_ww = np.zeros(points)

mut_bias_range = np.logspace(-1, 1, points)
for mut_bias_id, mut_bias in enumerate(mut_bias_range):
    mutation_matrix = np.array([[0, 1, 1, mut_bias],
                                [mut_bias, 0, 1, mut_bias],
                                [mut_bias, 1, 0, mut_bias],
                                [mut_bias, 1, 1, 0]])

    codon_frequencies = np.power(mut_bias, nbr_weak)
    codon_frequencies *= np.exp(fitness_codons)
    codon_frequencies /= np.sum(codon_frequencies)

    mut_bias_obs[mut_bias_id] = sum(codon_frequencies * z_ws) / sum(codon_frequencies * z_sw)

    dss, dss0, dsw, dsw0 = 0, 0, 0, 0
    dww, dww0, dws, dws0 = 0, 0, 0, 0
    for c_origin_id, codon_origin in enumerate(codons):
        for c_target_id, codon_target in enumerate(codons):
            mutations = [(n_from, n_to) for n_from, n_to in zip(codon_origin, codon_target) if n_from != n_to]
            if len(mutations) == 1:
                nuc_origin, nuc_target = mutations[0]
                n_origin_id = nucleotides.find(nuc_origin)
                n_target_id = nucleotides.find(nuc_target)
                sel_coef = fitness_codons[c_target_id] - fitness_codons[c_origin_id]
                if abs(sel_coef) < 1e-10:
                    p_fix = 1
                else:
                    p_fix = sel_coef / (1. - np.exp(-sel_coef))
                d_partiel = codon_frequencies[c_origin_id] * mutation_matrix[n_origin_id][n_target_id] * p_fix
                d0_partiel = codon_frequencies[c_origin_id] * mutation_matrix[n_origin_id][n_target_id]
                if nuc_origin == "A" or nuc_origin == "T":
                    if nuc_target == "A" or nuc_target == "T":
                        dww += d_partiel
                        dww0 += d0_partiel
                    else:
                        dws += d_partiel
                        dws0 += d0_partiel
                else:
                    if nuc_target == "A" or nuc_target == "T":
                        dsw += d_partiel
                        dsw0 += d0_partiel
                    else:
                        dss += d_partiel
                        dss0 += d0_partiel

    gamma_ss[mut_bias_id] = dss / dss0
    gamma_sw[mut_bias_id] = dsw / dsw0
    gamma_ws[mut_bias_id] = dws / dws0
    gamma_ww[mut_bias_id] = dww / dww0


plt.subplot(211)
plt.plot(mut_bias_range, gamma_ws/gamma_sw, '-',
         label='$\gamma_{w \\rightarrow s}/\gamma_{s \\rightarrow w}$')
plt.plot(mut_bias_range, mut_bias_range/mut_bias_obs, ':',
         label='$\lambda / \lambda_{obs}$')
plt.plot(mut_bias_range, np.ones(len(mut_bias_range)), color="black")
plt.xscale('log')
plt.ylim((0, max(mut_bias_range/mut_bias_obs)*1.5))

plt.ylabel('$\gamma$')
plt.legend()

plt.subplot(212)
plt.plot(mut_bias_range, mut_bias_range, alpha=0.75, label='$\lambda$')
plt.plot(mut_bias_range, mut_bias_obs, '--', label='$\lambda_{obs}$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\lambda$')
plt.ylabel('$\lambda_{obs}$')
plt.legend()
plt.show()
