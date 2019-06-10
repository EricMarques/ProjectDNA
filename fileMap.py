from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable
import codonTable as ct


gene = Seq("GGAGTACAATGATTCAGAACCCTGGGTGAGTGGAAGCCCCGGCAGATCGACAACCCAGAT" +
           "TACAAGGGCACTTGGATCCACCCAGAAATTGACAACCCCGAGTATTCTCCCGATCCCAGT" +
           "ATCTATGCCTATGATAGACTTCCTGCCACCCAAGAAGATAAAGGATCCTGATGCTTCAAA" +
           "ACCGGAAGACTGGGATGAGCGGGCCAAGATCGATGATCCCACAGACTCCAAGCCTGAGGA" +
           "CTGGGACAAGCCCGAGCATATCCCTGACCCTGAUGCTAAGAAGCCCGAGGACTGGGATGA" +
           "AGAGATGGACGGAGAGTGGGAACCCCCAACTTTGGCGTGCTGGGCCTGGACCTCTGGCAG" +
           "GTCAAGTCTGGCACCATCTTTGACAACTTCCTCATCACCAACGATGAGGCATACGCTGAG" +
           "GAGTTTGGCAACGAGACGTGGGGCGTAACAAAGGCAGCAGAGAAACAAATGAAGGACAAA" +
           "CAGGACGAGGAGCAGAGGCTTAAGGAGGAGGAAGAAGACAAGAAACGCAAAGAGGAGGAG" +
           "GAGGCAGAGGACAAGGAGGATGAA", generic_dna)

# print('DNA  = ' + gene)
# print('RNAm = ' + gene.transcribe())

# g = gene.translate(table='Standard')
# print('Prot = ' + g, '\n')
'''
for k, v in ct.codonTable.items():
    print(k, ' - ', v)
    if str(ct.startCodon) in gene:
        print('\n!!!!!!!!!!!!!1')
    print('NOPS')
'''
# standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
# print(standard_table)
print(ct.start_Codon_D)
print(ct.start_Codon_L)
start_c_l = str(ct.start_Codon_L).strip('[]').strip("'")
print(start_c_l)

condition = start_c_l in gene
print(condition)
print(gene.find(start_c_l) + 1)
# print(gene)
