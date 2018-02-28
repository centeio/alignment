# alignment
Sequence Alignment for Bioinformatica (Project 1)

## Objectives

1.  Implement global alignment algorithm
2.  Test against EBI REST Web Service to get the same result

## Installing libs

To install all libs required, run:

```
    $ pip install requirements.txt
```

## Testing alignment code

For testing the python code, run:

```
    $ python Local.py
```

If you want to save the output, run:

```
    $ python Local.py > outuput.txt
```

## Testing api consumption code

For testing this code, you must provide all the data needed inside the file, 
then you could run the code as above and also look for the "run\_on_ebi.log" file:
```
    $ python run_on_ebi.py
    or
    $ python run_on_ebi.py > ebi_result.txt
```

#### Resources:
[EMBOSS NEEDLE alignment tool by EBI](https://www.ebi.ac.uk/Tools/psa/emboss_needle/)

[Example of a protein Q9Y2Y8](http://www.uniprot.org/uniprot/Q9Y2Y8.fasta)

[Example of a protein B7Z8R9](http://www.uniprot.org/uniprot/B7Z8R9.fasta)