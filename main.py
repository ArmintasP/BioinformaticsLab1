import os
import math
from Bio import SeqIO, Data

codon_table = Data.CodonTable.unambiguous_dna_by_name["Standard"]


def frequency_per_1000(frequencies):
    relative_frequency = {}
    count = 0

    for k, v in frequencies.items():
        count += v

    # ratio = 1000 / count
    for k in frequencies:
        relative_frequency[k] = frequencies[k] * 1000 / count

    return relative_frequency


def get_codon_frequency(orfs):
    frequency = {codon: 0 for codon in codon_table.forward_table}
    return get_frequency(orfs, 3, frequency)


def get_dicodon_frequency(orfs):
    frequency = {}
    for codon_first in codon_table.forward_table:
        for codon_second in codon_table.forward_table:
            frequency[codon_first + codon_second] = 0

    return get_frequency(orfs, 6, frequency)


def get_frequency(strings, substring_length, frequency):
    for string in strings:
        for i in range(0, len(string), substring_length):
            substring = string[i: i + substring_length]
            if substring in frequency:
                frequency[substring] += 1
    return frequency


def get_open_reading_frames_both_strands(sequence, table, minimum_protein_length):
    orfs = []
    for seq in [sequence, sequence.reverse_complement()]:
        orfs += get_open_reading_frames(seq, table, minimum_protein_length)
    return orfs


def get_open_reading_frames(sequence, table, minimum_protein_length):
    orfs = []
    for frame in [0, 1, 2]:
        length = 3 * ((len(sequence) - frame) // 3)
        sequence = sequence[frame: frame + length]

        get_proteins(orfs, sequence, table, minimum_protein_length)
    return orfs


def get_proteins(proteins, sequence, table, minimum_protein_length):
    protein = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i: i + 3]
        if protein == "" and codon in table.start_codons:
            protein += codon
        elif protein != "":
            if codon in table.stop_codons:
                if len(protein) > minimum_protein_length:
                    proteins.append(protein)
                protein = ""
            else:
                protein += codon


def get_sequence_characteristics(sequence, sequence_name):
    filtered_orfs = get_open_reading_frames_both_strands(sequence, codon_table, minimum_protein_length=100)

    codon_frequency = get_codon_frequency(filtered_orfs)
    dicodon_frequency = get_dicodon_frequency(filtered_orfs)

    return sequence_name, frequency_per_1000(codon_frequency), frequency_per_1000(dicodon_frequency)


def get_euclidean_distance(frequencies1: dict, frequencies2: dict):
    total = 0
    for key in frequencies1:
        total += (frequencies1[key] - frequencies2[key]) ** 2

    return round(math.sqrt(total), 2)


def get_distance_matrices(characteristics):
    length = len(characteristics)

    dicodon_matrix = [[0 for _ in range(length + 1)] for _ in range(length)]
    codon_matrix = [[0 for _ in range(length + 1)] for _ in range(length)]

    for i in range(length):
        name, codon_frequency, dicodon_frequency = characteristics[i]
        codon_matrix[i][0] = name
        dicodon_matrix[i][0] = name

        for j in range(i + 1, length):
            name, codon_frequency2, dicodon_frequency2 = characteristics[j]

            distance = get_euclidean_distance(codon_frequency,
                                              codon_frequency2)

            codon_matrix[i][j + 1] = distance
            codon_matrix[j][i + 1] = distance

            distance = get_euclidean_distance(dicodon_frequency,
                                              dicodon_frequency2)

            dicodon_matrix[i][j + 1] = distance
            dicodon_matrix[j][i + 1] = distance

    return codon_matrix, dicodon_matrix


def write_to_file(matrix, file_name):
    with open("results/" + file_name, "w") as result_file:
        print(len(matrix), file=result_file)
        for row in matrix:
            print(" ".join(map(str, row)), file=result_file)


def main():
    sequence_characteristics = []

    for file in os.listdir("resources"):
        record = SeqIO.read("resources/" + file, "fasta")
        sequence_characteristics.append(get_sequence_characteristics(record.seq, record.name))
        print("Reading file: ", file, "...")

    codon_matrix, dicodon_matrix = get_distance_matrices(sequence_characteristics)
    write_to_file(codon_matrix, "codon-distance-matrix.txt")
    write_to_file(dicodon_matrix, "dicodon-distance-matrix.txt")
    print("Success! Distance matrices are saved to results/ directory.");


if __name__ == "__main__":
    main()
