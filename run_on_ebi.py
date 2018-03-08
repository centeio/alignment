import xml.etree.ElementTree as ET
import requests
import logging as LOG
import re
from datetime import datetime

"""
This script aims to consume EBI rest webservice for the needle
alignment algorithm

"""


#EBI URLS
EBI_RUN_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run/"
EBI_RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/"
REGEX = "(?<=\d )[AGVLYETDFQHICSMKPRNW\-]{1,50}"
MATCH_REGEX = "(?:Identity:)      (\d+)"
GAP_REGEX = "(?:Gaps:)         (\d+)"
MISMATCH_REGEX = "(?:Similarity:)    (\d+)"

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
    "asequence":">sp|Q9Y2Y8|PRG3_HUMAN Proteoglycan 3 OS=Homo sapiens GN=PRG3 PE=1 SV=2\nMQCLLLLPFLLLGTVSALHLENDAPHLESLETQADLGQDLDSSKEQERDLALTEEVIQAE\nGEEVKASACQDNFEDEEAMESDPAALDKDFQCPREEDIVEVQGSPRCKICRYLLVRTPKT\nFAEAQNVCSRCYGGNLVSIHDFNFNYRIQCCTSTVNQAQVWIGGNLRGWFLWKRFCWTDG\nSHWNFAYWSPGQPGNGQGSCVALCTKGGYWRRAQCDKQLPFVCSF",
    "bsequence":">tr|B7Z8R9|B7Z8R9_HUMAN cDNA FLJ57474, highly similar to Homo sapiens plasticity related gene 3 (PRG-3), transcript variant 1, mRNA OS=Homo sapiens PE=2 SV=1\nMAGTVLLAYYFECTDTFQVHIQGFFCQDGDLMKPYPGTEEESFITPLVLYCVLAATPTAI\nIFIGEISMYFIKSTRESLIAQEKTILTGECCYLNPLLRRIIRFTGVFAFGLFATDIFVNA\nGQVVTGHLTPYFLTVCKPNYTSADCQAHHQFINNGNICTGDLEVIEKARRSFPSKHAALS\nIYSALYATMYITSTIKTKSSRLAKPVLCLGTLCTAFLTGLNRVSEYRNHCSDVIAGFILG\nTAVALFLGMCVVHNFKGTQGSPSKPKPEDPRGVPLMAFPRIESPLETLSAQNHSASMTEV\nT"}
    

class Alignment():
    def __init__(self, seq1, seq2, gaps, matches, mismatches):
        self.gaps = gaps
        self.matches = matches
        self.mismatches = mismatches
        self.seq1 = seq1
        self.seq2 = seq2
    def __str__(self):
        return "[ALIGNMENT] \n\
        SEQ1: {} \n\
        SEQ2: {} \n\
        #gaps: {} \t #matches: {} \t #mismatches: {}".format(self.seq1,self.seq2,self.gaps,self.matches,self.mismatches)

def get_alignment(string):
    """
    summary: This fuctions receives a string returned by EBI and parses it 
            to get just the two sequences aligned.
    return: An alignment object
    """
    
    matches = re.compile(MATCH_REGEX).findall(string)[0] #getting the number of matches
    gaps = re.compile(GAP_REGEX).findall(string)[0] #getting the number of gaps
    mismatches = re.compile(MISMATCH_REGEX).findall(string)[0] #getting the number of mismatches
    
    string = re.sub(re.compile("#.*?\n" ) ,"" ,string) ##removing all comments
    p = re.compile(REGEX) #getting the REGEX variable
    m = p.findall(string) #returns all list of strings that match with that regex
    string1 = m[0::2] #gets all the even itens on the list (i.e the first sequence)
    string2 = m[1::2] #gets all the odds items on the list (i.e the second sequence)
    string1 = ''.join(string1) #join all items in a single string with '' sepring them
    string2 = ''.join(string2) #the same as above
    
    result = Alignment(string1,string2,gaps,matches,mismatches)
    
    return result


if __name__ == '__main__':
    #CONFIG
    LOG.basicConfig(filename='run_on_ebi.log',level=LOG.INFO)
    LOG.info("----------\n")
    path = 'ebi_output.txt'

    # RUN
    LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RUN_URL))
    r = requests.post(EBI_RUN_URL,data)
    LOG.info("[{0}] Request made and the code returned was {1}".format(datetime.now(),r.status_code))
    print ("Nome do job gerado = {}".format(r.text))

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

    f = open(path, 'w+')
    f.write(alignment.seq1)
    f.write("\n")
    f.write(alignment.seq2)
    f.close()
    
    print(str(alignment))

    print("File with the result in: {}".format(path))



