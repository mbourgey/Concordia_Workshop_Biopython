**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

# Concordia Python workshop - module 3 - Intorduction to Biopython   

by Mathieu Bourgey, _Ph.D_  

This Workshop is an adaptation of some interesting point of the [general Biopython tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)

## Learning objectives
During this wiorkshop you will learn:  

 - Manipulate Sequence Objects
 - Annotate Sequence Objects
 - Read/write Fasta file
 - Blasting Sequences
 - Searching over Pubmed database

 

## The Seq Object
Biological sequences are arguably the central object in Bioinformatics. Biopython sequences are essentially strings of letters like **AGTACACTGGT** as seen in common biological file formats.

The Biopython `Seq` object is defined in the `Bio.Seq` module

The `Seq` object is different from traditional python strings:  

 1. It supports most of the string methods but it also comes with its specifc set of methods
  * translate() _- Turns a nucleotide sequence into a protein sequence._
  * reverse_complement() _- Returns the reverse complement sequence._
  * complement() _- Returns the complement sequence._
  * transcribe() _-Returns the RNA sequence from a DNA sequence._
  * back_transcribe() _- Returns the DNA sequence from an RNA sequence._
  * ungap() _- Return a copy of the sequence without the gap character(s)._
 2. It has an important attribute, the alphabet, which is an object describing the type of the sequence and how the characters should be interpreted. Biopython alphabet are define in the `Bio.Alphabet` module

 
[The detail API of the `Seq` object](http://biopython.org/DIST/docs/api/Bio.Seq.Seq-class.html)


[The detail API of the `Alphabet` object](http://biopython.org/DIST/docs/api/Bio.Alphabet-module.html)
 
### Sequence Manipulation
We’ll use the [IUPAC](http://www.chem.qmw.ac.uk/iupac/) alphabets here to deal with some of our favorite objects: DNA, RNA and Proteins.

We can create an ambiguous sequence with the default generic alphabet like this:

```{.python}
from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
my_seq
```

> Seq('AGTACACTGGT', Alphabet())

```{.python}
my_seq.alphabet
```

> Alphabet()

In the above example, we haven't specified an alphabet so we end up with
a default generic alphabet. Biopython doesn't know if this is a
nucleotide sequence or a protein rich in alanines, glycines, cysteines
and threonines. If *you* know, you should supply this information

```{.python}
from Bio.Alphabet import IUPAC
my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
my_seq
```

> Seq('AGTACACTGGT', IUPACUnambiguousDNA())

```{.python}
my_seq.alphabet
```

> IUPACUnambiguousDNA()

####  Seq as python strings

In many ways, we can deal with Seq objects as if they were normal Python strings, for example getting the length, or iterating over the elements

```{.python}
for index, letter in enumerate(my_seq):
   print("%i %s" % (index, letter))

```

> 0 A   
> 1 G   
> 2 T   
> 3 A   
> 4 C   
> 5 A   
> 6 C   
> 7 T   
> 8 G   
> 9 G   
> 10 T   

```{.python}
len(my_seq)
```

> 11

The Seq object has a `.count()` method, just like a string. Note that this means that like a Python string, this gives a non-overlapping count

```{.python}
my_seq.count("A")
```

> 3

```{.python}
my_seq.count("GT")
```

> 2

Note that this means that like a Python string, this gives a non-overlapping count:

```{.python}
Seq("AAAA").count("AA")
```

> 2

####  Slicing a Seq object

let’s get a slice of the sequence

```{.python}
my_seq[2:8]
```

> Seq('TACACT', IUPACUnambiguousDNA())

Note that it follows the normal conventions for Python strings. 
 - So the first element of the sequence is 0 (which is normal for computer science, but not so normal for biology).
 - The first item is included i.e. 2 in this case (3rd aa)
 - The last is excluded i.e 8 in this case (9th aa)

###### Exercice

**In one command could you extract the third codon positions of this DNA sequence ?**

[solution](solutions/_codonExtract.md)

#### concatenate sequence
You can in principle add any two Seq objects together just like you can with Python strings. But `Seq` object are made for biological data so you the concatenation method only accept to merge sequences with compatible alphabets. You are allowed to concatenate a protein sequence and a DNA sequence.

```{.python}
p_seq = Seq("EVRNAK", IUPAC.protein)
d_seq = Seq('TACACT', IUPAC.unambiguous_dna)
d_seq + my_seq
```

> Seq('TACACTAGTACACTGGT', IUPACUnambiguousDNA())


```{.python}
p_seq + my_seq
```

> Traceback (most recent call last):   
> ...   
> TypeError: Incompatible alphabets IUPACProtein() and IUPACUnambiguousDNA()

**If you __realy__ want to do that, how should you do ?** [solution](solutions/_concate.md)

####  Seq as Biological strings

The `Seq` object is more than a python string with a specific alphabet, it also offers methods specific to facilitate the biology oriented analysis.

#### Complement and reverse complement
DNA is double stranded but in most case sequence are represented as single stranded molecules. for many purpose , i.e alignment, we need to compare a query sequence to reference sequence. In this case, we need to know if the reference sequence contains the query sequence in one or the other strands. 

You can easily obtain the complement or reverse complement of a Seq object using its built-in methods:

```{.python}
my_seq.complement()

```

> Seq('TCATGTGACCA', IUPACUnambiguousDNA())

```{.python}
my_seq.reverse_complement()

```

> Seq('ACCAGTGTACT', IUPACUnambiguousDNA())


There is no specific methods in Biopython to only reverse your sequence


**Do you know why and how to proceed ?**   
[solution](solutions/_reverse.md)


Note that these methods only work for dna alphabet. Trying to (reverse)complement a protein sequence will raise you an error:

```{.python}
p_seq = Seq("EVRNAK", IUPAC.protein)
p_seq.reverse_complement()

 ```
 
> Traceback (most recent call last):__ 
> ...   
> ValueError: Proteins do not have complements!

#### Transcirption, reverse transcription and translation 
First we need to clarify the transcription strand issue !   

![Transcription strand](img/Transcription_strand.png)

Biologically the transcription do a reverse complement of the template strand while inserting Uracile instead of Thymine (TCAG → CUGA) to give the RNA.

However, in Biopython and bioinformatics in general, we typically work directly with the coding strand because this means we can get the mRNA sequence just by switching T → U

Let's do a simple transcription of our sequence:

```{.python}
r_seq=my_seq.transcribe()
r_seq

```

> Seq('AGUACACUGGU', IUPACUnambiguousRNA())

And a reverse transcription of the resulting sequence:


```{.python}
r_seq.back_transcribe()

```
> Seq('AGTACACTGGT', IUPACUnambiguousDNA())

As you can see, all this does is switch T -> U or U -> T and adjust the alphabet.


###### Exercice

**Could you generate the mRNA from this template strand sequence: __3'-TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC-5'__ ?**  

[solution](solutions/_mRNA.md)


Now let’s translate this mRNA into the corresponding protein sequence 

```{.python}
p_seq = r_seq.translate()
p_seq
```

> Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))

Note that Biopython also allow to directly translate from DNA sequence. In this case use the coding strand DNA sequence.

You should also notice in the above protein sequences that in addition to the end stop character, there is an internal stop as well. This is due to the use of the wrong translation table in this case.

The translation tables available in Biopython are based on those from the [NCBI](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By default, translation will use the standard genetic code (NCBI table id 1). 

In our case, we are dealing with a mitochondrial sequence. We need to tell the translation function to use the relevant genetic code instead:

```{.python}
p_seq = r_seq.translate(table="Vertebrate Mitochondrial")
p_seq
```

> Seq('MAIVMGRWKGAR*', HasStopCodon(IUPACProtein(), '*'))

Also you may want to translate the nucleotides up to the first in frame stop codon, and then stop. In this case you should use to_stop method:


```{.python}
r_seq.translate(to_stop=True)
```

> Seq('MAIVMGR', IUPACProtein())

```{.python}
r_seq.translate(table=2, to_stop=True)
```

> Seq('MAIVMGRWKGAR', IUPACProtein())

Notice that when you use the to_stop argument, the stop codon itself is not translated, Notice also that you can specify the table using the NCBI table number which is shorter.

## The SeqRecord object
Immediately “above” the `Seq` class is the Sequence Record or `SeqRecord` class, defined in the `Bio.SeqRecord` module. This class allows higher level features to be associated with the sequence.

[The detail API of the `SeqRecord` object](http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html)

A SeqRecord object holds a sequence and information about it.

**Main attributes:**

        * .id - Identifier such as a locus tag or an accesion number (string)
        * .seq - The sequence itself (Seq object or similar)

**Additional attributes:**

        * .name - Sequence name, e.g. gene name (string)
        * .description - A human readable description or expressive name for the sequence (string)
        * .letter_annotations - Per letter/symbol annotation (restricted dictionary). This holds Python sequences (lists, strings or tuples) whose length matches that of the sequence. A typical use would be to hold a list of integers representing sequencing quality scores, or a string representing the secondary structure.
        * .features - Any (sub)features defined (list of `SeqFeature` objects), i.e location, type or strand...
        * .annotations - Further information about the whole sequence (dictionary). Most entries are strings, or lists of strings. This allows the addition of more “unstructured” information to the sequence.
        * .dbxrefs - List of database cross references (list of strings)

Using a `SeqRecord` object is not very complicated, since all of the information is presented as attributes of the class. Usually you won’t create a `SeqRecord` “by hand”, but instead you a specific Class to read in a sequence file for you (presented in the next section). However, creating `SeqRecord` can be quite simple.


### Manually creating a seqRecord

To create a `SeqRecord` at a minimum you just need a `Seq` object:

```{.python}
from Bio.SeqRecord import SeqRecord
simple_seq_r = SeqRecord(my_seq)
simple_seq_r

```
> SeqRecord(seq=Seq('AGTACACTGGT', IUPACUnambiguousDNA()), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])

We can also manually pass the id, name and description to the object

```{.python}
simple_seq_r.id = "THX1138"
simple_seq_r.name = "THX 1138 4EB"
simple_seq_r.description = "Made up sequence I wish I could write a paper about"
simple_seq_r

```

> SeqRecord(seq=Seq('AGTACACTGGT', IUPACUnambiguousDNA()), id='THX1138', name='THX 1138 4EB', description='Made up sequence I wish I could write a paper about', dbxrefs=[])

Including an identifier is very important if you want to output your SeqRecord to a file.

The `SeqRecord` has an dictionary attribute annotations. This is used for any miscellaneous annotations that doesn’t fit under one of the other more specific attributes. Adding annotations is easy, and just involves dealing directly with the annotation dictionary:

```{.python}
simple_seq_r.annotations["evidence"] = "None. I just made it up."
simple_seq_r.annotations

```
> {'evidence': 'None. I just made it up.'}

Working with per-letter-annotations is similar, letter_annotations is a dictionary like attribute which will let you assign any Python sequence (i.e. a string, list or tuple) which has the same length as the sequence.

```{.python}
import random
simple_seq_r.letter_annotations["phred_quality"] =  random.sample(xrange(1, 50),len(simple_seq_r))
simple_seq_r.letter_annotations

```

> {'phred_quality': [22, 23, 3, 2, 29, 11, 34, 44, 5, 33, 16]}

The `format()` method of the `SeqRecord` class gives a string containing your record formatted using one of the output file formats supported. 

```{.python}
simple_seq_r.format('fasta')

```
> '>THX1138 Made up sequence I wish I could write a paper about\nAGTACACTGGT\n'

```{.python}
simple_seq_r.format('fastq')
```

> '@THX1138 Made up sequence I wish I could write a paper about\nAGTACACTGGT\n+\n78$#>,CM&B1\n'

## The SeqIO Class
The `SeqIO` Class provide a simple interface for working with assorted sequence file formats in a uniform way.

The “catch” is that you have to work with `SeqRecord` objects, which contain a `Seq` object with format specific annotation like an identifier and description.


[The detail API of the `SeqIO` module](http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html)

### Reading or parsing sequence files 

The main method of the `SeqIO` module is __.parse()__  is used to read in sequence data as `SeqRecord` objects. This function expects two arguments:

    1. A handle to read the data from, or a filename.
    2. A lower case string specifying sequence format. See http://biopython.org/wiki/SeqIO for a full listing of supported formats. 
    
There is an optional argument alphabet to specify the alphabet to be used. This is useful for file formats like FASTA where otherwise Bio.SeqIO will default to a generic alphabet.

The __.parse()__ methods is typically used with a for loop like this:

```{.python}
from Bio import SeqIO
for seq_record in SeqIO.parse("data/NC_000913.fna","fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    print(len(seq_record.features))

```

> gi|556503834|ref|NC_000913.3|   
> Seq('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAG...TTC', SingleLetterAlphabet())   
> 4641652   
> 0


If instead you wanted to load a GenBank format:

```{.python}
for seq_record in SeqIO.parse("data/NC_000913.gbk","genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    print(len(seq_record.features))

```

> gi|556503834|ref|NC_000913.3|   
> Seq('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAG...TTC', SingleLetterAlphabet())   
> 4641652   
> 22996 


Notice that in the 2 examples the `SeqRecord` object contain the same central information and that mostly the difference is in the amount of information provided in the annotations (__.features__) 

Similarly, if you wanted to read in a file in another file format, you would just need to change the format string as appropriate. For example “swiss” for SwissProt files or “embl” for EMBL text files. There is a full listing on the [Biopython wiki page](http://biopython.org/wiki/SeqIO)

Another very common way to use a Python iterator is within a list comprehension

```{.python}
identifiers = [seq_record.id for seq_record in SeqIO.parse("data/patato_pep.fasta","fasta")]
identifiers

```

> ['PGSC0003DMP400067339', 'PGSC0003DMP400027454', 'PGSC0003DMP400020381', 'PGSC0003DMP400022612', 'PGSC0003DMP400040011', 'PGSC0003DMP400032361', 'PGSC0003DMP400030628', 'PGSC0003DMP400028584', 'PGSC0003DMP400060824', 'PGSC0003DMP400037883']

We usually use a for loop to iterate over all the records one by one. The object returned by `SeqIO` is actually an iterator which returns `SeqRecord` objects. You get to access each record in turn, but once and only once. The plus point is that an iterator can save you memory when dealing with large files.

Sometime you need to be able to access the records in any order. In that case you could acces the records using a list or dictionary process. But be carefull because you could easily carsh or slowdown python when using this methods on large data sets.

#### The list approach
We can turn the record iterator into a list of SeqRecord objects using the built-in Python function __list()__ 


```{.python}
records = list(SeqIO.parse("data/patato_pep.fasta","fasta"))
records[3]
```

> SeqRecord(seq=Seq('MLPFESIEEASMSLGRNLTFGETLWFNYTADKSDFYLYCHNTIFLIIFYSLVPL...SE*', SingleLetterAlphabet()), id='PGSC0003DMP400022612', name='PGSC0003DMP400022612', description='PGSC0003DMP400022612 PGSC0003DMT400033236 Protein', dbxrefs=[])

#### The dictionary approach
`SeqIO` module has three related functions which allow dictionary like random access to a multi-sequence file. There is a trade off here between flexibility and memory usage. In summary:  

   - `.to_dict()` is the most flexible but also the most memory demanding option. This is basically a helper function to build a normal Python dictionary with each entry held as a SeqRecord object in memory, allowing you to modify the records.
   -  `.index()` is a useful middle ground, acting like a read only dictionary and parsing sequences into SeqRecord objects on demand.
   - `.index_db()` also acts like a read only dictionary but stores the identifiers and file offsets in a file on disk (as an SQLite3 database), meaning it has very low memory requirements, but will be a little bit slower. 


Let's test the different functions

You can use the function `.to_dict()` to make a SeqRecord dictionary (in memory). By default this will use each record’s identifier (i.e. the .id attribute) as the key

```{.python}
records_dict = SeqIO.to_dict(SeqIO.parse("data/patato_pep.fasta","fasta"))
records_dict.keys()

```

> ['PGSC0003DMP400030628', 'PGSC0003DMP400032361', 'PGSC0003DMP400027454', 'PGSC0003DMP400060824', 'PGSC0003DMP400040011', 'PGSC0003DMP400037883', 'PGSC0003DMP400022612', 'PGSC0003DMP400020381', 'PGSC0003DMP400067339', 'PGSC0003DMP400028584']

```{.python}
records_dict['PGSC0003DMP400020381']
```

> SeqRecord(seq=Seq('MLEKDSRDDRLDCVFPSKHDKDSVEEVSSLSSENTRTSNDCSRSNNVDSISSEV...KY*', SingleLetterAlphabet()), id='PGSC0003DMP400020381', name='PGSC0003DMP400020381', description='PGSC0003DMP400020381 PGSC0003DMT400029984 Protein', dbxrefs=[])

`.to_dict()` is very flexible because it holds everything in memory. The size of file you can work with is limited by your computer’s RAM. In general, this will only work on small to medium files.

---------------

For larger files you should consider `.index()`, which works a little differently. Although it still returns a dictionary like object, this does not keep everything in memory. Instead, it just records where each record is within the file and when you ask for a particular record, it then parses it on demand.


```{.python}
records_dict = SeqIO.index("data/patato_pep.fasta","fasta")
list(records_dict.keys())

```

> ['PGSC0003DMP400030628', 'PGSC0003DMP400032361', 'PGSC0003DMP400027454', 'PGSC0003DMP400060824', 'PGSC0003DMP400040011', 'PGSC0003DMP400037883', 'PGSC0003DMP400022612', 'PGSC0003DMP400020381', 'PGSC0003DMP400067339', 'PGSC0003DMP400028584']

Note that in this case the `.keys()` function return an iterator and we need to use the list function to get the key values

```{.python}
records_dict['PGSC0003DMP400020381']
```

> SeqRecord(seq=Seq('MLEKDSRDDRLDCVFPSKHDKDSVEEVSSLSSENTRTSNDCSRSNNVDSISSEV...KY*', SingleLetterAlphabet()), id='PGSC0003DMP400020381', name='PGSC0003DMP400020381', description='PGSC0003DMP400020381 PGSC0003DMT400029984 Protein', dbxrefs=[])

Note that `.index()` won’t take a handle, but only a filename. 

-----------------

`.index_db()` work on even extremely large files since it stores the record information as a file on disk (using an SQLite3 database) rather than in memory. Also, we can index multiple files together (providing all the record identifiers are unique).

`.index_db()` function takes three required arguments:

    - Index filename
    - List of sequence filenames to index (or a single filename)
    - File format (lower case string as used in the rest of the SeqIO module). 
    
```{.python}
patato_pep = SeqIO.index_db("patato_pep.idx", "data/patato_pep.fasta","fasta")
patato_pep.keys()
```

> ['PGSC0003DMP400020381', 'PGSC0003DMP400022612', 'PGSC0003DMP400027454', 'PGSC0003DMP400028584', 'PGSC0003DMP400030628', 'PGSC0003DMP400032361', 'PGSC0003DMP400037883', 'PGSC0003DMP400040011', 'PGSC0003DMP400060824', 'PGSC0003DMP400067339']

```{.python}
patato_pep['PGSC0003DMP400040011']
```

> SeqRecord(seq=Seq('MECDTEDSEDNSNIQADSNHRLVKFVIPGNNLLDQTKSSSTKVVLIFLESVEIL...NF*', SingleLetterAlphabet()), id='PGSC0003DMP400040011', name='PGSC0003DMP400040011', description='PGSC0003DMP400040011 PGSC0003DMT400059441 Protein', dbxrefs=[])


**So, which of these methods should you use and why ?** [solution](solutions/_seqIO1.md) 

### Writing sequence files

The `.write()` is used to output sequences (writing files). This function taking three arguments: 

 - Some SeqRecord objects. 
 - A handle or filename to write to. 
 - A sequence format.

First, let's write a sequence into the file __testOut.fa__ :

```{.python}
import os
SeqIO.write(simple_seq_r, "testOut.fa",  "fasta")
os.system("cat testOut.fa")

```
> \>THX1138 Made up sequence I wish I could write a paper about   
> AGTACACTGGT

then, let's write a set of patato sequences into the file __testOut.fa__ :

```{.python}
for seq_record in SeqIO.parse("data/patato_pep.fasta","fasta") : 
  SeqIO.write(seq_record, "testOut.fa",  "fasta")

os.system("cat testOut.fa")

```
> \>PGSC0003DMP400037883 PGSC0003DMT400056292 Protein   
> MMIGRDPEIWENPEEFIPERFLNSDIDYFKGQNFELIPFGAGRRGCPGIALGVATINLIL   
> SNLLYAFDWELPHGMIKEDIDTDGLPGLAMHKKNALCLVPKNYTHT*


**Do you notice something strange in testOut.fa ? and explain why ?** [solution](solutions/_seqIO2.md)


###### Exercice

**could you write the content data/patato_pep.fasta into testOut.fa ?** [solution](solutions/_seqIO3.md)


## The Blast Class
Everybody loves BLAST ! How can it get any easier to do comparisons between one of your sequences and every other sequence in the known world?

If you don't know Blast please explore http://blast.ncbi.nlm.nih.gov/Blast.cgi

Dealing with BLAST can be split up into two steps, both of which can be done from within Biopython:

 1. Running BLAST for your query sequence(s), and getting some output. 
 2. Parsing the BLAST output in Python for further analysis.

To do that Biopython provide the specific `Blast` module.

[The detail API of the `Blast` module](http://biopython.org/DIST/docs/api/Bio.Blast-module.html)
 
There are lots of ways you can run BLAST, especially you can run BLAST locally (on your own machine), or run BLAST remotely (on another machine, typically the NCBI servers).

In this workshop we will only focus on the remote way to run blast.


### Running BLAST on NCBI servers
The specific sub-module `Blast.NCBIWWW` allows to call the online version of BLAST through is main function `qblast()`.

This function has three mandatory arguments:

 1. The blast program to use for the search, as a lower case string. Currently qblast only works with blastn, blastp, blastx, tblast and tblastx.
 2. The databases to search against, as a lower case string.
 3. The query sequence as a string. This can either be the sequence itself, the sequence in fasta format, or an identifier like a GI number.
 
Note that the default settings on the NCBI BLAST website are not quite the same as the defaults on QBLAST. If you get different results, you’ll need to check the parameters (e.g., the expectation value threshold and the gap values).
 
```{.python}
from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastp", "nr", patato_pep['PGSC0003DMP400040011'].seq)
```

Supplying just the sequence means that BLAST will assign an identifier for your sequence automatically. You might prefer to use the SeqRecord object’s format method to make a FASTA string (which will include the existing identifier):

```{.python}
result_handle = NCBIWWW.qblast("blastn", "nr", patato_pep['PGSC0003DMP400040011'].format("fasta"))
```

### Writing and reading
Whatever arguments you give the qblast() function, you should get back your results in a handle object (by default in XML format). The next step would be to parse the XML output into Python objects representing the search results, but you might want to save a local copy of the output file first. 

To output the blast result we use the `read()` function:

```{.python}
save_file = open("my_blast.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()
os.system("head my_blast.xml")
```

> \<?xml version="1.0"?\>   
> \<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">   
> \<BlastOutput\>   
> &nbsp;&nbsp;\<BlastOutput_program\>blastp\</BlastOutput_program\>   
> &nbsp;&nbsp;\<BlastOutput_version\>BLASTP 2.3.1+\</BlastOutput_version\>   
> &nbsp;&nbsp;\<BlastOutput_reference\>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.\</BlastOutput_reference\>   
>  &nbsp;&nbsp;\<BlastOutput_db\>nr\</BlastOutput_db\>   
>  &nbsp;&nbsp;\<BlastOutput_query-ID\>Query_56965\</BlastOutput_query-ID\>   
>  &nbsp;&nbsp;\<BlastOutput_query-def\>PGSC0003DMP400040011 PGSC0003DMT400059441 Protein\</BlastOutput_query-def\>   
>  &nbsp;&nbsp;\<BlastOutput_query-len\>222\</BlastOutput_query-len\>   



## The Entrez Object

