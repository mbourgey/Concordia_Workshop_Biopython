**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

# Concordia Python workshop - module 3 - Intorduction to Biopython   

by Mathieu Bourgey, _Ph.D_  

This Workshop is an adaptation of the [general Biopython tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)

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

#### complement and reverse complement
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


**Do you know why and how to preoceed ?** [solution](solutions/_reverse.md)


Note that these methods only works for dna alphabet property is maintained. Trying to (reverse)complement a protein sequence will raise you an error:

```{.python}
p_seq = Seq("EVRNAK", IUPAC.protein)
p_seq.reverse_complement()

 ```
 
> Traceback (most recent call last):__ 
> ...   
> ValueError: Proteins do not have complements!

