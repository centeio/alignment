import requests
import logging as LOG
import re
from datetime import datetime
from .util.alignment import Alignment
from .util import io

"""
This script aims to consume EBI rest webservice for the needle
alignment algorithm

"""


#EBI URLS
EBI_NEEDLE_RUN_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run/"
EBI_NEEDLE_RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/"
EBI_WATER_RUN_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_water/run/"
EBI_WATER_RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_water/result/"
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
    "asequence":"A",
    "bsequence":"A"
      }
    

class Ebi():
    def __init__(self,type="GLOBAL"):
        self.req_dict = {
            "LOCAL":EBI_WATER_RUN_URL,
            "GLOBAL":EBI_NEEDLE_RUN_URL
        }
        self.result_dict = {
            "LOCAL":EBI_WATER_RESULT_URL,
            "GLOBAL":EBI_NEEDLE_RESULT_URL
        }
        self.request_url = self.req_dict[type]
        self.result_url = self.result_dict[type]
        self.type_dict = {
            1:"protein",
            2:"dna"
        }
        self.matrix_dict = {
            1:"EBLOSUM62",
            2:"EPAM250",
            3:"EDNAFULL"
        }
        

    def run(self, seq1, seq2, seqtype, matrix, gapopen, gapext):
        #CONFIG
        # LOG.basicConfig(filename='../log/run_on_ebi.log',level=LOG.INFO)
        LOG.info("----------\n")

        data['stype'] = self.type_dict[seqtype]
        data['matrix'] = self.matrix_dict[matrix]
        data['gapopen'] = str(abs(gapopen))
        data['gapext'] = str(abs(gapext))
        data['asequence'] = seq1
        data['bsequence'] = seq2
        

        # RUN
        LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),self.request_url))
        r = requests.post(self.request_url,data)
        LOG.info("[{0}] Request made and the code returned was {1}".format(datetime.now(),r.status_code))
        print ("Nome do job gerado = {}".format(r.text))

        #RESULT
        # LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RESULT_URL))
        #Addind the jobid
        result = requests.get(self.result_url+"{jobID}/aln".format(jobID=r.text))
        
        #Trying again until 200 is returned
        while (result.status_code is not 200):
            # LOG.info("[{0}] Result not ready yet...".format(datetime.now()))
            # LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RESULT_URL))
            result = requests.get(self.result_url+"{jobID}/aln".format(jobID=r.text))
        
        LOG.info("[{0}] Results got = {1}...".format(datetime.now(),result.text))

        alignment = self.get_alignment(result.text)

        return alignment

    def get_alignment(self, string):
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
        
        result = Alignment(string1,string2)
        result.calculate_mat_mis_gaps()
        
        return result







