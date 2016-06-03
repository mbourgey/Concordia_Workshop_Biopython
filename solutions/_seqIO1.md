It depends on what you are trying to do (and how much data you are dealing with) !

In general picking `.index()` is a good starting point. But if you are dealing with millions of records, multiple files, or repeated analyses, then `.index_db()` is probably a better choice.

Reasons to choose `.to_dict()` over `.index()` or `.index_db()`:  
   
   - Very small dataset.
   - Flexibility - `SeqRecord` objects in memory is they can be changed, added to, or removed.
   - Speed - Despite its high memory needs having everything in RAm is way much faster

Reasons to choose `.index()` over `.index_db()`:  
  
   - Faster to build 
   - Faster to access as SeqRecord objects 
   - Can use any immutable Python object as the dictionary keys (e.g. a tuple of strings, or a frozen set) not just strings.
   - Don’t need to worry about the index database being out of date if the sequence file being indexed has changed. 

Reasons to choose `.index_db()` over `.index()`:  

 - Not memory limited - This very important for actual NGS data
 - Because the index is kept on disk, it can be reused. Although building the index database file takes longer, if you have a script which will be rerun on the same datafiles in future, this could save time in the long run. 
 - Indexing multiple files together 