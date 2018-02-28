import xml.etree.ElementTree as ET
import requests
import logging as LOG
from datetime import datetime

"""
This script aims to consume EBI rest webservice for the needle
alignment algorithm

"""


#EBI URLS
EBI_RUN_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run/"
EBI_RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/"

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
    

if __name__ == '__main__':
    #CONFIG
    LOG.basicConfig(filename='run_on_ebi.log',level=LOG.INFO)
    LOG.info("----------\n")

    # RUN
    LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RUN_URL))
    r = requests.post(EBI_RUN_URL,data)
    LOG.info("[{0}] Request made and the code returned was {1}".format(datetime.now(),r.status_code))
    print "Nome do job gerado = {}".format(r.text)

    #RESULT
    LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RESULT_URL))
    #Addind the jobid
    result = requests.get(EBI_RESULT_URL+"{jobID}/aln".format(jobID=r.text))
    
    #Trying again until 200 is returned
    while (result.status_code is not 200):
        LOG.info("[{0}] Result not ready yet...".format(datetime.now()))
        LOG.info("[{0}] Requesting the URL = {1}...".format(datetime.now(),EBI_RESULT_URL))
        result = requests.get(EBI_RESULT_URL+"{jobID}/aln".format(jobID=r.text))
    
    print result.text


