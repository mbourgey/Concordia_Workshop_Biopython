when you using a handler, it is, in fact, just a connectiowhich point to a specific position in an open file. So when the file have been read entirely or if you just want to point to the a specific position you could use the method seek to change the position of the pointer.

`.seek(0)` just replace the pointer to the begining of the file and is equivalent as doing `.close()` then `.open()`.

