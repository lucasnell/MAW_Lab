`convert_to_Glazer.R`
===========

R script to convert from older to newer Glazer stickleback assembly
---------



To this script you provide a single tab-delimited text file of positions on the old
assembly. The first column must be the chromosome, the second column the position.
It doesn't matter what you call the columns in the header, but there MUST be a header.

Things to keep in mind:

- The function will fail if the input file doesn't have a `.txt` extension. (I did 
  this intentionally so it didn't overwrite the input file!).
- The input file needs to have the chromosomes in the first column, positions in the 
  second. Column names don't matter, except for the fact that the output table will 
  have the same as the input.
- This script has to be located in the same folder as `convert_to_Glazer_table.RData`.
- You obviously have to have R installed on your system.


You would run this script as such for an input file named `in_positions.txt`:

```bash
./convert_to_Glazer.R in_positions.txt
```

It will output a file named `<input file name w/o extension>_NEW.txt` 
(`in_positions_NEW.txt` in this case) containing the positions on the Glazer assembly.

