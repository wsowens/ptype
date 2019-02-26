'''Command line API for IDT's OligoAnalyzer tool
https://www.idtdna.com/calc/analyzer
'''
import requests
import json

TOOL_URLS = {
    "analyzer" : "https://www.idtdna.com/calc/analyzer/home/analyze",
    "hairpin" : "https://www.idtdna.com/calc/analyzer/home/hairpin",
    "selfDimer" : "https://www.idtdna.com/calc/analyzer/home/selfDimer",
    "heteroDimer" : "https://www.idtdna.com/calc/analyzer/home/heteroDimer",
    "tmMismatch" : "https://www.idtdna.com/calc/analyzer/home/tmMismatch"
}

def analyze(sequence, nucleotide_type, dnac1, Na, Mg, dNTPs, **kwargs):
    '''Function corresponding with the main "analyze" function
    Arguments:
    sequence - string containing the sequence of interest
    nucleotide_type - string with valid nucleotide type, options are 'DNA' and 'RNA'
    dnac1 - concenctration of the oligo (in uM)
    Na - concentration of Na+ (in uM)
    Mg - concentration of Mg2+ (in uM)
    dNTPs - concentration of dNTPs (in uM)

    Returns a dictionary with the following values:
    '''
    request_data = {
        'sequence': sequence,
        'nucleotideType': nucleotide_type,
        'oligoConc': dnac1,
        'naConc': Na,
        'mgConc': Mg,
        'dNTPsConc': dNTPs
    }
    result = requests.post(TOOL_URLS["analyzer"], data=request_data)
    #if result has error, raise with status code
    if not result.ok:
        raise Exception("Request Error: [%s]" % result.status_code )
    result = json.loads(result.text)
    # check for errors identified by the analyzer
    if "HasModelErrors" in result and result["HasModelErrors"]:
        Exception()
    # delete any irrelevant information
    for k in []:
        if k in result:
            del result[k]
    return result


HAIRPIN_DEFAULTS = {
    "SequenceType": "Linear",
    "MaxFoldings": 20,
    "StartPos": 0, 
    "StopPos": 0,
    "Temp": 25, 
    "Suboptimality": 50,  
}

# work in progress function
def hairpin(sequence, nucleotide_type, dnac1, Na, Mg, dNTPs, hairpin_settings=HAIRPIN_DEFAULTS):
    request_data = {
        'sequence': sequence,
        'nucleotideType': nucleotide_type,
        'oligoConc': dnac1,
        'naConc': Na,
        'mgConc': Mg,
        'dNTPsConc': dNTPs,
        'hairpinSettings': hairpin_settings
    }
    print("sending hairpin")
    result = requests.post(TOOL_URLS["hairpin"], data=request_data)
    #if result has error, raise with status code
    if not result.ok:
        raise Exception("Request Error: [%s]" % result.status_code )
    result = json.loads(result.text)
    # check for errors identified by the analyzer
    if "HasModelErrors" in result and result["HasModelErrors"]:
        raise Exception(result["ModelErrors"])
    # delete any irrelevant information
    for k in []:
        if k in result:
            del result[k]
    return result


#print(analyze("AACTAGTACTAGTAGTACA", "DNA", 0.25, 0.5, 2.5, 1))
#print(hairpin("AACTAGTACTAGTAGTACA", "DNA", 0.25, 0.5, 2.5, 1))
