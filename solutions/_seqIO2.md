Only the last sequence have been write in the file.

The reason is `.write()` does not work as the traditional python write() function. `SeqIO.write()` open, write and close the file at each call. So If you tell the `SeqIO.write()` function to write to a file that already exists, the old file will be overwritten without any warning.

