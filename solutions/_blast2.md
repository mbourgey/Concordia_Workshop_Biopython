The quantity of BLAST results. 

If our BLAST file is huge though, we may run into memory problems trying to save them all in a list.

Usually, we’ll be running one BLAST search at a time. Then, all we need to do is to pick up the first BLAST record in blast_records:

```{.python}
result_handle.seek(0)
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)
````

