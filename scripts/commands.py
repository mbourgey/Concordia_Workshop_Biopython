from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
my_seq
my_seq.alphabet
from Bio.Alphabet import IUPAC
my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
my_seq
my_seq.alphabet
for index, letter in enumerate(my_seq):
   print("%i %s" % (index, letter))

len(my_seq)
my_seq.count("A")
my_seq.count("GT")
Seq("AAAA").count("AA")
my_seq[2:8]
p_seq = Seq("EVRNAK", IUPAC.protein)
d_seq = Seq('TACACT', IUPAC.unambiguous_dna)
d_seq + my_seq
p_seq + my_seq
my_seq.complement()

my_seq.reverse_complement()

p_seq = Seq("EVRNAK", IUPAC.protein)
p_seq.reverse_complement()

r_seq=my_seq.transcribe()
r_seq

r_seq.back_transcribe()

p_seq = r_seq.translate()
p_seq
p_seq = r_seq.translate(table="Vertebrate Mitochondrial")
p_seq
r_seq.translate(to_stop=True)
r_seq.translate(table=2, to_stop=True)
from Bio.SeqRecord import SeqRecord
simple_seq_r = SeqRecord(my_seq)
simple_seq_r

simple_seq_r.id = "THX1138"
simple_seq_r.name = "THX 1138 4EB"
simple_seq_r.description = "Made up sequence I wish I could write a paper about"
simple_seq_r

simple_seq_r.annotations["evidence"] = "None. I just made it up."
simple_seq_r.annotations

import random
simple_seq_r.letter_annotations["phred_quality"] =  random.sample(xrange(1, 50),len(simple_seq_r))
simple_seq_r.letter_annotations

simple_seq_r.format('fasta')

simple_seq_r.format('fastq')
from Bio import SeqIO
for seq_record in SeqIO.parse("data/NC_000913.fna","fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    print(len(seq_record.features))

for seq_record in SeqIO.parse("data/NC_000913.gbk","genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    print(len(seq_record.features))

identifiers = [seq_record.id for seq_record in SeqIO.parse("data/patato_pep.fasta","fasta")]
identifiers

records = list(SeqIO.parse("data/patato_pep.fasta","fasta"))
records[3]
records_dict = SeqIO.to_dict(SeqIO.parse("data/patato_pep.fasta","fasta"))
records_dict.keys()

records_dict['PGSC0003DMP400020381']
records_dict = SeqIO.index("data/patato_pep.fasta","fasta")
list(records_dict.keys())

records_dict['PGSC0003DMP400020381']
patato_pep = SeqIO.index_db("patato_pep.idx", "data/patato_pep.fasta","fasta")
patato_pep.keys()
patato_pep['PGSC0003DMP400040011']
import os
SeqIO.write(simple_seq_r, "testOut.fa",  "fasta")
os.system("cat testOut.fa")

for seq_record in SeqIO.parse("data/patato_pep.fasta","fasta") : 
  SeqIO.write(seq_record, "testOut.fa",  "fasta")

os.system("cat testOut.fa")

from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastp", "nr", seq_record.seq)
