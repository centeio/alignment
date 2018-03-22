import os
fileDir = os.path.dirname(os.path.realpath('__file__'))

def write_file(path,data):
    filename = os.path.join(fileDir, path)
    f = open(filename, 'w+')
    f.write(str(data))
    f.close()

def read_file(path):
    filename = os.path.join(fileDir, path)
    f = open(filename, 'r+')
    result = f.read()
    f.close()
    return result

def read_fasta(seq):
    
    if seq[0] is not ">":
        return seq
    
    seq = seq.splitlines()[1:]
    seq = "".join(seq).replace(" ", "").replace("\r", "").replace("\n","")
    return seq
