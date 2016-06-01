You’d have to first give both sequences generic alphabets:

```{.python}
from Bio.Alphabet import generic_alphabet
d_seq.alphabet = generic_alphabet
p_seq.alphabet = generic_alphabet
p_seq + d_seq
```

> Seq('EVRNAKTACACT', Alphabet())

