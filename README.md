# alignment
Sequence Alignment for Bioinformatica (Project 1)

## Objectives

1.  Implement global and local alignment algorithms (Needleman-Wunsch and Smith-Waterman) using PAM and BLOSUM matrices for scores calculations. (Linear GAP) **Points: 4/8**

2.  Consume EBI REST Web Service with same input and compare both alignments **Points: 3/8**

3. Implement affine GAP calculations **Points: 1/8** 

## Project Strucute

```
alignment
├── inputs
│   ├── B7Z8R9.fasta
│   ├── default1.fasta
│   ├── default2.fasta
│   ├── NM_002688.fasta
│   ├── NP_002679.fasta
│   ├── Q9Y2Y8.fasta
│   ├── XM_001166286.fasta
│   └── XP_009215056.fasta
├── log
│   └── run_on_ebi.log
├── outputs
│   ├── ebi_global_output.txt
│   ├── ebi_local_output.txt
│   ├── locally_global_affine_output.txt
│   ├── locally_global_linear_output.txt
│   ├── locally_local_affine_output.txt
│   └── locally_local_linear_output.txt
├── README.md
├── requirements.txt
└── src
    ├── AffineGlobal.py
    ├── AffineLocal.py
    ├── GlobalEBI.py
    ├── Global.py
    ├── LocalEBI.py
    ├── Local.py
    └── util
        ├── alignment.py
        ├── __init__.py
        ├── io.py
        └── matrices.py
```

## Installing libs

To install all libs required, run:

```
    $ pip install requirements.txt
```

## Testing locally

For testing our implementation:

```
    # cd src/
    $ python3 {file}.py
```

It will (by default), look for dummy1 and dummy2 fasta files on `inputs` directory. You can either replace the content of these two files with your sequences or change the code to read from another file. To do so, you must change these lines of code:

`{file}.py`
```python
153     
154        #Sequencias
155        seq1=io.read_fasta(io.read_file("../inputs/<YOUR_FILE>"))
156        seq2 =io.read_fasta(io.read_file("../inputs/<YOUR_FILE>"))
157    
```

This code will produce the alignment in the `{output}.txt` file inside outputs directory, as shown:

`{output}.txt`
```
[ALIGNMENT] 
        SEQ1: C---LF- 
        SEQ2: AAACLFF 
        #gaps: 4 	 #matches: 2 	 #mismatches: 1
```

## Testing remotelly using EBI API

For testing this code, you must provide all the data needed inside the file `run_on_ebi.py`. The input works exactly as above. If you want to change the parameters, just changes the following part:

`{ebi}.py`
```python
19
20      #EXAMPLE OF REQUEST
21      data = {"email":"up201711183@fc.up.pt",
22      "matrix":"EBLOSUM62",
23      "gapopen":"10.0",
24      "gapext":"0.5",
25      "endweight":"false",
26      "endopen":"10.0",
27      "endextend":"0.5",
28      "format":"pair",
29      "stype":"protein",
30      "asequence":io.read_file("../inputs/dummy1.fasta"),
31      "bsequence":io.read_file("../inputs/dummy2.fasta")
32       }
33    
```

This code will produce the alignment in the `{ebi_output}.txt` file inside outputs directory.

#### Resources:
[EMBOSS NEEDLE alignment tool by EBI](https://www.ebi.ac.uk/Tools/psa/emboss_needle/)

[Example of a protein Q9Y2Y8](http://www.uniprot.org/uniprot/Q9Y2Y8.fasta)

[Example of a protein B7Z8R9](http://www.uniprot.org/uniprot/B7Z8R9.fasta)