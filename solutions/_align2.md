Severals approaches could be used.

-----------------

The lazzy one: use another format

```{.python}
AlignIO.write(my_alignments, "mixed.fa", "fasta")
os.system("head mixed.fa")

```

> \>toto <unknown description>   
> ACTGCTAGCTAG   
> \>titi <unknown description>   
> ACT-CTAGCTAG   
> \>tata <unknown description>   
> \ACTGCTAGDTAG   
> \>PGSC0003DMP400022612   
> ------------------------------------------------------------   
> ------------------------------------------------------------   
> -------MLPFESIEEASMSLGRNLTFGETLWF---------------------------   



But the lazzy is usually not the best one (the fasta format is ambigous) and probably we had good reason to choose the 'phylip' format initialy! 

---------------------------
The reformat approach

Let's just extract the variable part of the id (the end) end paste it in the begining of the id

```{.python}
for seqR in my_alignments[1]:
   seqR.id = seqR.id[-10:] + "-" + seqR.id

print my_alignments[1][1]

```

> ID: P400060824-PGSC0003DMP400060824   
> Name: <unknown name>   
> Description: PGSC0003DMP400060824   
> Number of features: 0   
> Seq('---------------MQIFVKTLTGKTITLEVESSDTIDNVKAKIQDKEGIPPD...GGF', SingleLetterAlphabet())   


```{.python}
AlignIO.write(my_alignments, "mixed.phy", "phylip")
os.system("head mixed.phy")

```

> &nbsp; 3 &nbsp; 12   
> toto       ACTGCTAGCT AG   
> titi       ACT-CTAGCT AG   
> tata       ACTGCTAGDT AG   
> &nbsp; 10 &nbsp; 1294   
> P400022612 ---------- ---------- ---------- ---------- ----------   
> P400060824 ---------- -----MQIFV KTLTGKTITL EVESSDTIDN VKAKIQDKEG   
> P400067339 ---------- -----MGVWK DSNYGKGVII GVIDT--GIL PDHPSFSDVG   
> P400027454 MPHPTQVVAL LKAQQIRHVR LFNADRGMLL ALANTGIKVA VSVPNEQILG   
> P400028584 ---------- --------MS TSVEPNGAVL ---------- ----LDSTAG   


