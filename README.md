**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

# Concordia Python workshop - module 3 - Intorduction to Biopython   

 
by Mathieu Bourgey, _Ph.D_

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

Seq('AGTACACTGGT', Alphabet())

```{.python}
my_seq.alphabet
```

Alphabet()

In the above example, we haven't specified an alphabet so we end up with
a default generic alphabet. Biopython doesn't know if this is a
nucleotide sequence or a protein rich in alanines, glycines, cysteines
and threonines. If *you* know, you should supply this information

```{.python}
from Bio.Alphabet import IUPAC
my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
my_seq
```

``Seq('AGTACACTGGT', IUPACUnambiguousDNA())``

```{.python}
my_seq.alphabet
```

``IUPACUnambiguousDNA()``
