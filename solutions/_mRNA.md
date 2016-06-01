If you do want to do a true biological transcription starting with the template strand, then this becomes a two-step process: 

##### Method 1

```{.python}
reverse_template_dna = Seq('TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC', IUPAC.unambiguous_dna)
template_dna = reverse_template_dna[::-1]
template_dna
```

> Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT', IUPACUnambiguousDNA())

```{.python}
template_dna.reverse_complement().transcribe()
```

> Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())

##### Method 2 - shortcut

```{.python}
reverse_template_dna = Seq('TACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATC', IUPAC.unambiguous_dna)
reverse_template_dna.complement().transcribe()
```

> Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())

