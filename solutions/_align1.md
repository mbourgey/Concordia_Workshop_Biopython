The last line of the Trace back indicate:

> ValueError: Repeated name 'PGSC0003DM' (originally 'PGSC0003DMP400060824'), possibly due to truncation

Which means that for the 'phylip' format only a 10 characters string could be use as id. And if the `SeqRecord` id is longer only the first 10 charater will be retain. In our case the 10 first character of the peptide sequences names are identical which won't work ! 



