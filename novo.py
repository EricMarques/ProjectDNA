import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

seq = "GCGCGCGTCCGTCCGTACTGCAGAGCCGCTGCCGGAGGGTCGTTTTAAAGGGCCCGCGCGTTGCCGCCCCC" \
	  "TCGGCCCGCCATGCTGCTATCCGTGCCGCTGCTGCTCGGCCTCCTCGGCCTGGCCGTCGCCGAGCCTGCC" \
	  "GTCTACTTCAAGGAGCAGTTTCTGGACGGAGGTAACGCCTGGTCCCGCCTCGAGGCCGCCCCGACGACGC" \
	  "GGCCGGCCCCCGATCCTGGATCTGCGTTGTCGCCCGTAATTACCGTTTAGAGGTCCAACACGGTGGCCTC" \
	  "CCGGGACTAGAGCCGCGGGCGATTTCTCTTCTGCGTCCCTGGGGAGCGCGGAGGGCGTAGCGGCCTCCCG" \
	  "CGGCGGGAGTTAGGGTTAGCCCGAGGATCTCTGAAGGCACCCGACGTGTCAAACTAGAGGTTGGAATGGG" \
	  "GAGTGTCGGGGATCTCCTTTCCTGTCCCCAGCAGCTTGTGGCTCTCGGCAGATGTTTGGTGTGGGGGGGG" \
	  "ATTAGCACAGCCGCTCTGACCTACCCCTCTAATCCCCCACTTAGACGGGTGGACTTCCCGCTGGATCGAA" \
	  "TCCAAACACAAGTCAGATTTTGGCAAATTCGTTCTCAGTTCCGGCAAGTTCTACGGTGACGAGGAGAAAG" \
	  "ATAAAGGTAAGAGCCTAGGAGTGGGTGCTCAGATCCGGGAGGACTTCCTGGCAGAAGTCCTTGTCTGTAC" \
	  "ACACACAGCCGGGACAGTCCCCTTGGAGGAGGACAGGTGGAGGAAGTGGGGGAGTCTTCTCTATTCTCTA" \
	  "AGTCGAGGGTCCTCGCGAGTCAAGGCCCAACGGTGACCTCACTACCGTCCCGTCTCAGGTTTGCAGACAA" \
	  "GCCAGGATGCACGCTTTTATGCTCTGTCGGCCAGTTTCGAGCCTTTCAGCAACAAAGGCCAGACGCTGGT" \
	  "GGTGCAGTTCACGGTGAAACATGAGCAGAACATCGACTGTGGGGGCGGCTATGTGAAGCTGTTTCCTAAT" \
	  "AGTTTGGACCAGACAGACATGCACGGAGACTCAGAATACAACATCATGTTTGGTGAGGGCCTGCTTCCTG" \
	  "GTGCTGATCTCTGTCCCATTAGTTAGAGGGAGACCCAGACCCCATTGACTTTCTTAATAATGATTTTTTT" \
	  "TGGAAGGGGAGCTAAAAGAATAAGTCCCAGCAACAATTTATTGCATTATGATCGCAGATCTAGGCTGTTA" \
	  "ATTTAATTTGCGTGTTTGTATATAGTTATTTCCCAATCTTACTAATGAGGATTTTGAGTTCTAGAGCACT" \
	  "GATTTTTTTTTTTTCTCCTTTAAACTTAAGGCTCCACCCACAGCCCATTCAGGACAGAATCAGGGTCTGA" \
	  "GTTTCTCTTCTCAGCCTTGACAGACCCGAGTTGAAGAACCAGGTCTTCCTTTTATAAAGAGGGGTGAGAG" \
	  "CCTCGAGATGATGGGTAGTCTCTGACTCTTAACTGGATCTGCTTCACACCTAGGTCCCGACATCTGTGGC" \
	  "CCTGGCACCAAGAAGGTTCATGTCATCTTCAACTACAAGGGCAAGAACGTGCTGATCAACAAGGACATCC" \
	  "GTTGCAAGGTGTGCCTGGGGGTGGTGGCAAATGGCTGTCATGGGGAGATTCAGAGGTCAGCCTCATTGGG" \
	  "GGGTGGCCCCCGCTCACCTTCTTCCTTCTTCAGGATGATGAGTTTACACACCTGTACACACTGATTGTGC" \
	  "GGCCAGACAACACCTATGAGGTGAAGATTGACAACAGCCAGGTGGAGTCCGGCTCCTTGGAAGACGATTG" \
	  "GGACTTCCTGCCACCCAAGAAGATAAAGGATCCTGATGCTTCAAAACCGGAAGACTGGGATGAGCGGGCC" \
	  "AAGATCGATGATCCCACAGACTCCAAGCCTGAGGTTGGTGTTTGGGCAGGGGCTCTGCTCTCCACATTGG" \
	  "AGGGTGTGGAAGACATCTGGGCCAACTCTGATCTCTTCATCTACCCCCCAGGACTGGGACAAGCCCGAGC" \
	  "ATATCCCTGACCCTGATGCTAAGAAGCCCGAGGACTGGGATGAAGAGATGGACGGAGAGTGGGAACCCCC" \
	  "AGTGATTCAGAACCCTGAGTACAAGGTGAGTTTGGGGCTCTGAGCAGGGCTGGGGCTCACAGTGGGGAGT" \
	  "GCACCAACCTTACTCACCCTTCGGTTTCCTTCTCCCTTCTGCAGGGTGAGTGGAAGCCCCGGCAGATCGA" \
	  "CAACCCAGATTACAAGGGCACTTGGATCCACCCAGAAATTGACAACCCCGAGTATTCTCCCGATCCCAGT" \
	  "ATCTATGCCTATGATAACTTTGGCGTGCTGGGCCTGGACCTCTGGCAGGTGAGACTTGGAGGAAAAAGGA" \
	  "GGATCCCTGGGGTACCTCAAGTGCATAAGATCACCCAAGAGGAAAGGGACAGGGTAGGCACCCCAGGTGA" \
	  "GTCTGACTCAAAAATGGTACTTCTTGTAAACAGTACTTCCTGGTCTGTCCCTGTGAAGTCCTCACAGCAA" \
	  "CCCCTTTAAGGTTATACTTGCTGTGCACCAAGTACTTCCCCAAGTACTTTTATGCAAATCAACTTCTTTA" \
	  "CCCCCAAAGACCTAGAAGGTGGTCAGGTAACCCAGTTAGTTAGCTGGGGCTGGGCACAGTGGCTCACCCT" \
	  "TACAATCACGGTACTTTGGGAGGCTGAGACAGAGGATTGCTTGAGGCCAGGAGTTACACAACTCAACCTA" \
	  "GCTTGGCAACACAGCGAGGAGACCCTATCTCTACAAAAAAAATTTTTTTTTTTGAGACAGAGTTTCACTC" \
	  "TTGTTGCTGAGGCTGGAGTGCAATGGCACGATCTCAGCTCACTGCGCCCTCCGTCTCCTGGTTTCAAGCG" \
	  "ATTCTCCTGCCTCAGCCTCCGGAGTAGCTGGGATTACAGGCATGTGCTACTATGGATGCCAGGCTAATTT" \
	  "TTTTTTTTTTTTTTTTTTTTGAGACCGTGCCTTGCTCTGTCGCCCAGGCTGGAGTGCAGTGGTGTGATCT" \
	  "CTGCTCACTGCAAGCTCCGCACGACCCCCCAGGTTCACTCCATTCTTCTGCCTCAGGGTCCCGAGTAACT" \
	  "GGGACTACAGGCACCCCCCACCATGCCTGGCTAATTTTTTTGTATTTTTTTTTTTAGTACAGACATGGTT" \
	  "TCACCGTGTTAGCCAGGATGGTCTCCATCTCCTGACCTCATGAACCACCCACCTTGGCCTCCCAAAGTGC" \
	  "TGGGATTACAGGCGTGAGCCACCTCACCCAGCCTTTTTGTAGAGACAGGGCTTCATGTTGCCCAGGTTGG" \
	  "TCTCGAACTCCTGGCCTCAGGTCATCTGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAAGGGTTAGC" \
	  "CACCATGCCTAGCCTCTACAAAAACTTTAAAAATTGGCGAGATGTCATGCATACCTGTAGTCCCAACTAC" \
	  "CAAGGAAGAAGGATGATCACTTGAGCCTGGGGCATCGAGGCTGCAGTGAGCCATGATTATGTCACTGCAC" \
	  "TCCAGCCTCGGTGACAGAGTGAGACCCTCTCAAAAAAAGTTGGGACTTGGCCGGACACAGTGGCTCACAC" \
	  "CTGTAATCCCAGCACTTTGGGAGGCCAAGGCGGGTGGATCACAAGGTCAGGAGATGGAGACCATCCTGGC" \
	  "TAACATGGTGAATGAAACCCCATCTCTAGTAAAAATACAAAAAATTTGCCAGGTGTGGTGGTGGGCGCCT" \
	  "GTAGTCCCAGCTACTCGGGAGGCTGAGGCAAAAGGATGACGTGAACCCGGGAGGCGGAGCTTGCAGTGAG" \
	  "CTGAGATCATGCCATTGCACTCCAGCCTGGGTGATAGCGAGACTCTGTCCCAAAAAAAAAAAAAAATGCT" \
	  "GGGACTGAATTTTTGTCTGTTTTGGTCACTGAAATACCTTCTGTGCCCAAGACAGTTCTGGCATGTAGTA" \
	  "GGTACCTGAAAAATACCTGAATAAGAGAGTGAGAAACAAGAAACAGGTGCAGAGAACTGAAGTCAGTGGC" \
	  "CCAAGGTCATGGGGGTAGGAAACCACAAAGCTGGGGTTTGAACCTGGGCAGTACAGCACCTGAGTCTCTC" \
	  "CATCTTTTTTTTTTTTTTTTTTTAAGACAGAGTCTTGCTCTGTCACCCAGGTTGGAGTGCAGTGGCTTGA" \
	  "TCTCGGCTCACTGCAGCCTCTGCCTTCCAGGTTCAAGTGATTCTCATGCCTCATCCTCTCGAGCAGCTGG" \
	  "AATTACAGGCATGCGCCACGACGCTGGGCTTTTTTTTTTTTGAGATGGAATTTCACTCTTGTTGCCCAGG" \
	  "CTGGAGTGCAATGATGCAATCTCGGCGGCTCACCACAACCTCTGCATCCCAGATTCAAGCGATTCTCCTG" \
	  "CCTCGGCCTCCTGAGTAGCTGGGATTACAGGGATGCGCCATCACAGACCCCGGGCTAATTTTTTTTAGTA" \
	  "GAGACAGAGTTTCACTATGTTGCCCAGGTTGGTCTCGAACTCCTGGCCTCAAGTGATCCGTTCGCCATGA" \
	  "CCTCCCAAAGTGCTGGGATTACAGGCATGAGCCCGTCCCGTCCCTGGCTGTCTCTCCATCTTTCCATCTT" \
	  "TTTTTTTTTTTTTTTTTTTTTTGGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCACGA" \
	  "TCTTGGCTCACTGCAAGCTCCGCCTCCTGGGTTCACATCATTCTCCTGTCTCAGCCTCCCAAATAGCTGG" \
	  "GACTACAGGCACTTGCCACCACGCCTGGCTGATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGT" \
	  "TAGCCAGGGTGGTCTCGATCTCCTGACCTCGTGATCCGCCCACCTTGGCCTCTGGGCGAGGATTACAGGC" \
	  "GTGATCCACCTCACCTGGCCTCTCCATCTTTTTAACTGCAGTGTCAGCGGTGTTCCTTGTCTTCTCTGCA" \
	  "GATGCAGGCAGCAGAATATAGTGGTTATAGGAACACAGGTGGAAACCCTGTCCAAAGCAAGGGCTATCGG" \
	  "GTATCACCTCTGACCATCCTTCCCATTCATCCTCCAGGTCAAGTCTGGCACCATCTTTGACAACTTCCTC" \
	  "ATCACCAACGATGAGGCATACGCTGAGGAGTTTGGCAACGAGACGTGGGGCGTAACAAAGGTGAGGCCTG" \
	  "GTCCTGGTCCTGATGTCGGGGGCGGGCAGGGCTGGCAGGGGGCAAGGCCCTGAGGTGTGTGCTCTGCCTG" \
	  "CAGGCAGCAGAGAAACAAATGAAGGACAAACAGGACGAGGAGCAGAGGCTTAAGGAGGAGGAAGAAGACA" \
	  "AGAAACGCAAAGAGGAGGAGGAGGCAGAGGACAAGGAGGATGATGAGGACAAAGATGAGGATGAGGAGGA" \
	  "TGAGGAGGACAAGGAGGAAGATGAGGAGGAAGATGTCCCCGGCCAGGCCAAGGACGAGCTGTAGAGAGGC" \
	  "CTGCCTCCAGGGCTGGACTGAGGCCTGAGCGCTCCTGCCGCAGAGCTGGCCGCGCCAAATAATGTCTCTG" \
	  "TGAGACTCGAGAACTTTCATTTTTTTCCAGGCTGGTTCGGATTTGGGGTGGATTTTGGTTTTGTTCCCCT" \
	  "CCTCCACTCTCCCCCACCCCCTCCCCGCCCTTTTTTTTTTTTTTTTTTAAACTGGTATTTTATCTTTGAT" \
	  "TCTCCTTCAGCCCTCACCCCTGGTTCTCATCTTTCTTGATCAACATCTTTTCTTGCCTCTGTCCCCTTCT" \
	  "CTCATCTCTTAGCTCCCCTCCAACCTGGGGGGCAGTGGTGTGGAGAAGCCACAGGCCTGAGATTTCATCT" \
	  "GCTCTCCTTCCTGGAGCCCAGAGGAGGGCAGCAGAAGGGGGTGGTGTCTCCAACCCCCCAGCACTGAGGA" \
	  "AGAACGGGGCTCTTCTCATTTCACCCCTCCCTTTCTCCCCTGCCCCCAGGACTGGGCCACTTCTGGGTGG" \
	  "GGCAGTGGGTCCCAGATTGGCTCACACTGAGAATGTAAGAACTACAAACAAAATTTCTATTAAATTAAAT" \
	  "TTTGTGTCTCC"

'''
#1
seg = "GCG"

# segmento = re.compile(r"{}".format(seg))
segmento = re.compile(r"(?={})".format(seg))

for padrao in segmento.finditer(seq):
       print(padrao.group(), padrao.start())
'''


'''
seq2 = "gcggcgtccgtccgtactgcagagccgctgccggagggtcgttttaaagggcccgcgcgt" \
"tgccgccccctcggcccgccatgctgctatccgtgccgctgctgctcggcctcctcggcc" \
"tggccgtcgccgagcctgccgtctacttcaaggagcagtttctggacggaggtaacgcct" \
"ggtcccgcctcgaggccgccccgacgacgcggccggcccccgatcctggatctgcgttgt" \
"cgcccgtaattaccgtttagaggtccaacacggtggcctcccgggactagagccgcgggc" \
"gatttctcttctgcgtccctggggagcgcggagggcgtagcggcctcccgcggcgggagt" \
"tagggttagcccgaggatctctgaaggcacccgacgtgtcaaactagaggttggaatggg" \
"gagtgtcggggatctcctttcctgtccccagcagcttgtggctctcggcagatgtttggt" \
"gtggggggggattagcacagccgctctgacctacccctctaatcccccacttagacgggt" \
"ggacttcccgctggatcgaatccaaacacaagtcagattttggcaaattcgttctcagtt" \
"ccggcaagttctacggtgacgaggagaaagataaaggtaagagcctaggagtgggtgctc" \
"agatccgggaggacttcctggcagaagtccttgtctgtacacacacagccgggacagtcc" \
"ccttggaggaggacaggtggaggaagtgggggagtcttctctattctctaagtcgagggt" \
"cctcgcgagtcaaggcccaacggtgacctcactaccgtcccgtctcaggtttgcagacaa" \
"gccaggatgcacgcttttatgctctgtcggccagtttcgagcctttcagcaacaaaggcc" \
"agacgctggtggtgcagttcacggtgaaacatgagcagaacatcgactgtgggggcggct" \
"atgtgaagctgtttcctaatagtttggaccagacagacatgcacggagactcagaataca" \
"acatcatgtttggtgagggcctgcttcctggtgctgatctctgtcccattagttagaggg" \
"agacccagaccccattgactttcttaataatgattttttttggaaggggagctaaaagaa" \
"taagtcccagcaacaatttattgcattatgatcgcagatctaggctgttaatttaatttg" \
"cgtgtttgtatatagttatttcccaatcttactaatgaggattttgagttctagagcact" \
"gattttttttttttctcctttaaacttaaggctccacccacagcccattcaggacagaat" \
"cagggtctgagtttctcttctcagccttgacagacccgagttgaagaaccaggtcttcct" \
"tttataaagaggggtgagagcctcgagatgatgggtagtctctgactcttaactggatct" \
"gcttcacacctaggtcccgacatctgtggccctggcaccaagaaggttcatgtcatcttc" \
"aactacaagggcaagaacgtgctgatcaacaaggacatccgttgcaaggtgtgcctgggg" \
"gtggtggcaaatggctgtcatggggagattcagaggtcagcctcattggggggtggcccc" \
"cgctcaccttcttccttcttcaggatgatgagtttacacacctgtacacactgattgtgc" \
"ggccagacaacacctatgaggtgaagattgacaacagccaggtggagtccggctccttgg" \
"aagacgattgggacttcctgccacccaagaagataaaggatcctgatgcttcaaaaccgg" \
"aagactgggatgagcgggccaagatcgatgatcccacagactccaagcctgaggttggtg" \
"tttgggcaggggctctgctctccacattggagggtgtggaagacatctgggccaactctg" \
"atctcttcatctaccccccaggactgggacaagcccgagcatatccctgaccctgatgct" \
"aagaagcccgaggactgggatgaagagatggacggagagtgggaacccccagtgattcag" \
"aaccctgagtacaaggtgagtttggggctctgagcagggctggggctcacagtggggagt" \
"gcaccaaccttactcacccttcggtttccttctcccttctgcagggtgagtggaagcccc" \
"ggcagatcgacaacccagattacaagggcacttggatccacccagaaattgacaaccccg" \
"agtattctcccgatcccagtatctatgcctatgataactttggcgtgctgggcctggacc" \
"tctggcaggtgagacttggaggaaaaaggaggatccctggggtacctcaagtgcataaga" \
"tcacccaagaggaaagggacagggtaggcaccccaggtgagtctgactcaaaaatggtac" \
"ttcttgtaaacagtacttcctggtctgtccctgtgaagtcctcacagcaacccctttaag" \
"gttatacttgctgtgcaccaagtacttccccaagtacttttatgcaaatcaacttcttta" \
"cccccaaagacctagaaggtggtcaggtaacccagttagttagctggggctgggcacagt" \
"ggctcacccttacaatcacggtactttgggaggctgagacagaggattgcttgaggccag" \
"gagttacacaactcaacctagcttggcaacacagcgaggagaccctatctctacaaaaaa" \
"aatttttttttttgagacagagtttcactcttgttgctgaggctggagtgcaatggcacg" \
"atctcagctcactgcgccctccgtctcctggtttcaagcgattctcctgcctcagcctcc" \
"ggagtagctgggattacaggcatgtgctactatggatgccaggctaattttttttttttt" \
"ttttttttttgagaccgtgccttgctctgtcgcccaggctggagtgcagtggtgtgatct" \
"ctgctcactgcaagctccgcacgaccccccaggttcactccattcttctgcctcagggtc" \
"ccgagtaactgggactacaggcaccccccaccatgcctggctaatttttttgtatttttt" \
"tttttagtacagacatggtttcaccgtgttagccaggatggtctccatctcctgacctca" \
"tgaaccacccaccttggcctcccaaagtgctgggattacaggcgtgagccacctcaccca" \
"gcctttttgtagagacagggcttcatgttgcccaggttggtctcgaactcctggcctcag" \
"gtcatctgcccgcctcggcctcccaaagtgctgggattacaagggttagccaccatgcct" \
"agcctctacaaaaactttaaaaattggcgagatgtcatgcatacctgtagtcccaactac" \
"caaggaagaaggatgatcacttgagcctggggcatcgaggctgcagtgagccatgattat" \
"gtcactgcactccagcctcggtgacagagtgagaccctctcaaaaaaagttgggacttgg" \
"ccggacacagtggctcacacctgtaatcccagcactttgggaggccaaggcgggtggatc" \
"acaaggtcaggagatggagaccatcctggctaacatggtgaatgaaaccccatctctagt" \
"aaaaatacaaaaaatttgccaggtgtggtggtgggcgcctgtagtcccagctactcggga" \
"ggctgaggcaaaaggatgacgtgaacccgggaggcggagcttgcagtgagctgagatcat" \
"gccattgcactccagcctgggtgatagcgagactctgtcccaaaaaaaaaaaaaaatgct" \
"gggactgaatttttgtctgttttggtcactgaaataccttctgtgcccaagacagttctg" \
"gcatgtagtaggtacctgaaaaatacctgaataagagagtgagaaacaagaaacaggtgc" \
"agagaactgaagtcagtggcccaaggtcatgggggtaggaaaccacaaagctggggtttg" \
"aacctgggcagtacagcacctgagtctctccatctttttttttttttttttttaagacag" \
"agtcttgctctgtcacccaggttggagtgcagtggcttgatctcggctcactgcagcctc" \
"tgccttccaggttcaagtgattctcatgcctcatcctctcgagcagctggaattacaggc" \
"atgcgccacgacgctgggcttttttttttttgagatggaatttcactcttgttgcccagg" \
"ctggagtgcaatgatgcaatctcggcggctcaccacaacctctgcatcccagattcaagc" \
"gattctcctgcctcggcctcctgagtagctgggattacagggatgcgccatcacagaccc" \
"cgggctaattttttttagtagagacagagtttcactatgttgcccaggttggtctcgaac" \
"tcctggcctcaagtgatccgttcgccatgacctcccaaagtgctgggattacaggcatga" \
"gcccgtcccgtccctggctgtctctccatctttccatctttttttttttttttttttttt" \
"ttggagatggagtctcactctgtcacccaggctggagtgcagtggcacgatcttggctca" \
"ctgcaagctccgcctcctgggttcacatcattctcctgtctcagcctcccaaatagctgg" \
"gactacaggcacttgccaccacgcctggctgattttttgtatttttagtagagacggggt" \
"ttcaccgtgttagccagggtggtctcgatctcctgacctcgtgatccgcccaccttggcc" \
"tctgggcgaggattacaggcgtgatccacctcacctggcctctccatctttttaactgca" \
"gtgtcagcggtgttccttgtcttctctgcagatgcaggcagcagaatatagtggttatag" \
"gaacacaggtggaaaccctgtccaaagcaagggctatcgggtatcacctctgaccatcct" \
"tcccattcatcctccaggtcaagtctggcaccatctttgacaacttcctcatcaccaacg" \
"atgaggcatacgctgaggagtttggcaacgagacgtggggcgtaacaaaggtgaggcctg" \
"gtcctggtcctgatgtcgggggcgggcagggctggcagggggcaaggccctgaggtgtgt" \
"gctctgcctgcaggcagcagagaaacaaatgaaggacaaacaggacgaggagcagaggct" \
"taaggaggaggaagaagacaagaaacgcaaagaggaggaggaggcagaggacaaggagga" \
"tgatgaggacaaagatgaggatgaggaggatgaggaggacaaggaggaagatgaggagga" \
"agatgtccccggccaggccaaggacgagctgtagagaggcctgcctccagggctggactg" \
"aggcctgagcgctcctgccgcagagctggccgcgccaaataatgtctctgtgagactcga" \
"gaactttcatttttttccaggctggttcggatttggggtggattttggttttgttcccct" \
"cctccactctcccccaccccctccccgcccttttttttttttttttttaaactggtattt" \
"tatctttgattctccttcagccctcacccctggttctcatctttcttgatcaacatcttt" \
"tcttgcctctgtccccttctctcatctcttagctcccctccaacctggggggcagtggtg" \
"tggagaagccacaggcctgagatttcatctgctctccttcctggagcccagaggagggca" \
"gcagaagggggtggtgtctccaaccccccagcactgaggaagaacggggctcttctcatt" \
"tcacccctccctttctcccctgcccccaggactgggccacttctgggtggggcagtgggt" \
"cccagattggctcacactgagaatgtaagaactacaaacaaaatttctattaaattaaat" \
"tttgtgtctcc"

seq3 = "CAACTTTGGCGTGCTGGGCCTGGACCTCTGGCAGTGAGCGGGCCA" \
	   "AGATCGATGATCCCACAGACTCCAAGCCTGAGGACTGGGACAAGCCCGAGCATATCCCTG" \
	   "GTCAAGTCTGACAACTTCCTCATCACCAACGATGAGGCATAGCACCATCTTTGCGCTGAG" \
	   "GAGTTTGGCAACGAGACGTGGGGCGTAACAAAGGCAGCAGAGAAACAAATGAAGGACAAA" \
	   "CAGGAGAGCAAAGCGGAGCAGAAGACAAGAAACAGGAGGAGGAGGCTTAAAGGAGGAG" \
	   "GAGAGGGGATGAGACGCG"

# query = "mavmaprtlvlllsgalaltqtwagshsmryfftsvsrpgrgeprfiavgyvddtqfvrfdsdaasqrmeprapwieqegpeywdgetrkvkahsqthrvdlgtlrgyynqseagshtvqrmygcdvgsdwrflrgyhqyaydgkdyialkedlrswtaadmaaqttkhkweaahvaeqlraylegtcvewlrrylengketlqrtdapkthmthhavsdheatlrcwalsfypaeitltwqrdgedqtqdtelvetrpagdgtfqkwaavvvpsgqeqrytchvqheglpkpltlrwepssqptipivgiiaglvlfgavitgavvaavmwrrkssdrkggsysqaassdsaqgsdvsltackv"


# subject = "mwtlvswvaltaglvagtrcpdgqfcpvaccldpggasysccrplldkwpttlsrhlggpcqvdahcsaghsciftvsgtssccpfpeavacgdghhccprgfhcsadgrscfqrsgnnsvgaiqcpdsqfecpdfstccvmvdgswgccpmpqasccedrvhccphgafcdlvhtrcitptgthplakklpaqrtnravalsssvmcpdarsrcpdgstccelpsgkygccpmpnatccsdhlhccpqdtvcdliqskclskenattdlltklpahtvgdvkcdmevscpdgytccrlqsgawgccpftqavccedhihccpagftcdtqkgtceqgphqvpwmekapahlslpdpqalkrdvpcdnvsscpssdtccqltsgewgccpipeavccsdhqhccpqgytcvaegqcqrgseivaglekmparraslshprdigcdqhtscpvgqtccpslggswaccqlphavccedrqhccpagytcnvkarscekevvsaqpatflarsphvgvkdvecgeghfchdnqtccrdnrqgwaccpyrqgvccadrrhccpagfrcaargtkclrreaprwdaplrdpalrqll"

'''
# print("Fasta: {}".format(len(str.strip(seq))))
# print("SEQ3: {}".format(len(str.strip(seq3))))


# 3
def long_mottif(data):
	substr = ''
	if len(data) > 1 and len(data[0]) > 0:
		for i in range(len(data[0])):
			for j in range(len(data[0])-i + 1):
				if j > len(substr) and all(data[0][i:i+j] in x for x in data):
					substr = data[0][i:i+j]
	return substr


path = "Samples/CALR/Proteins/"
dirs = os.listdir(path)

seq3 = open(path + "ProteinsAll.fasta", "r")

content = seq3.read()

lines = []

for data in content.split('>'):
	line = ''.join(data.split('\n')[1:]).replace('\n', '')
	if line:
		lines.append(line)

tr = Seq(long_mottif(lines), IUPAC.unambiguous_dna)
print(lines)
print("Seq Comum: {} \nTamanho: {} caracteres" .format(long_mottif(lines), len(long_mottif(lines))))
# print("Prot: {}\nTamanho: {} caracteres" .format(tr.translate(), len(tr.translate())))

seq3.close()
'''
#2
for seq_record in SeqIO.parse(seq3, 'fasta'):
	seqDna = seq_record.seq

se = str(seqDna)

if len(seqDna) % 3 == 0:
	print("Sequência Original: {}" .format(seqDna))
	print("Sequência Translate: {}" .format(seqDna.translate()))
	print("Len: {}".format(len(seqDna.translate())))
elif (len(seqDna) + 1) % 3 == 0:
	seq3N = Seq('N') + seqDna
	print("Sequência 1 N: {}" .format(seq3N))
	print("Sequencia 1 N:    {}" .format(seq3N[3:]))
	print("Sequência Translate: {}" .format(seq3N.translate()))
	print("Len: {}".format(len(seq3N.translate())))
elif (len(seqDna) + 2) % 3 == 0:
	seq3N = (2*Seq('N')) + seqDna
	print("Sequência 2 N: {}".format(seq3N))
	print("Sequencia 2 N:    {}".format(seq3N[3:]))
	print("Sequência Translate: {}".format(seq3N.translate()))
	print("Len: {}".format(len(seq3N.translate())))
'''
'''

GGAGTACAATGATTCAGAACCCTGGGTGAAGGAGTGCCCCGGCAGATCGACAACCCAGAT
TACAAGGGCACTTGGATCCACCCAGAAATTGACAACCCCGAGTCCCGATCCCAGTTCAAA
ATCTATGCCTATGATAGACTTTATTCCCTGCCACCCAAGAAGATAAAGGATCCTGATGCT
ACCGGAAGACTGGGAACCCTGATGCTAAGAAGCCCGAGGACTGGGATGAAGAGATGGACG
GAGAGTGGGAACCCCCAACTTTGGCGTGCTGGGCCTGGACCTCTGGCAGTGAGCGGGCCA
AGATCGATGATCCCACAGACTCCAAGCCTGAGGACTGGGACAAGCCCGAGCATATCCCTG
GTCAAGTCTGACAACTTCCTCATCACCAACGATGAGGCATAGCACCATCTTTGCGCTGAG
GAGTTTGGCAACGAGACGTGGGGCGTAACAAAGGCAGCAGAGAAACAAATGAAGGACAAA
CAGGAGAGCAAAGCGAGGAGCAGAAGACAAGAAACAGGAGGAGGAGGCTTAAAGGAGGAG
GAGAGGGGATGAGACAAAGGCAG
'''
