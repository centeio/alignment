import requests
import logging as LOG
import re
from datetime import datetime
from util.alignment import Alignment
from util import io

"""
This script aims to consume EBI rest webservice for the needle
alignment algorithm

"""


#EBI URLS
EBI_RUN_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run/"
EBI_RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/"
REGEX = "(?<=\d )[AGVLYETDFQHICSMKPRNW\-]{1,50}"

#EXAMPLE OF REQUEST
data = {"email":"up201711183@fc.up.pt",
    "matrix":"EBLOSUM62",
    "gapopen":"10.0",
    "gapext":"0.5",
    "endweight":"false",
    "endopen":"10.0",
    "endextend":"0.5",
    "format":"pair",
    "stype":"protein",
    "asequence":io.read_file("../inputs/default1.fasta"),
    "bsequence":io.read_file("../inputs/default2.fasta")
      }
    

def get_alignment(string):
    """
    summary: This fuctions receives a string returned by EBI and parses it 
            to get just the two sequences aligned.
    return: An alignment object
    """
    
    string = re.sub(re.compile("#.*?\n" ) ,"" ,string) ##removing all comments
    p = re.compile(REGEX) #getting the REGEX variable
    m = p.findall(string) #returns all list of strings that match with that regex
    string1 = m[0::2] #gets all the even itens on the list (i.e the first sequence)
    string2 = m[1::2] #gets all the odds items on the list (i.e the second sequence)
    string1 = ''.join(string1) #join all items in a single string with '' sepring them
    string2 = ''.join(string2) #the same as above
    
    result = Alignment(string1,string2,"DEFAULT")
    result.calculate_mat_mis_gaps()
    
    return result


if __name__ == '__main__':
    #CONFIG
    LOG.basicConfig(filename='../log/run_on_ebi.log',level=LOG.INFO)
    LOG.info("----------\n")

    # RUN
    LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RUN_URL))
    r = requests.post(EBI_RUN_URL,data)
    LOG.info("[{0}] Request made and the code returned was {1}".format(datetime.now(),r.status_code))
    LOG.info("Nome do job gerado = {}".format(r.text))

    #RESULT
    LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RESULT_URL))
    #Addind the jobid
    result = requests.get(EBI_RESULT_URL+"{jobID}/aln".format(jobID=r.text))
    
    #Trying again until 200 is returned
    while (result.status_code is not 200):
        LOG.info("[{0}] Result not ready yet...".format(datetime.now()))
        LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RESULT_URL))
        result = requests.get(EBI_RESULT_URL+"{jobID}/aln".format(jobID=r.text))
    
    LOG.info("[{0}] Results got = {1}...".format(datetime.now(),result.text))

    alignment = get_alignment(result.text)

    io.write_file("../outputs/ebi_global_output.txt",str(alignment))
    
    #print(str(alignment))




