'''Command line API for IDT's OligoAnalyzer tool
https://www.idtdna.com/calc/analyzer
'''
import requests
import json

TOOL_URLS = {
    "analyzer" : "https://www.idtdna.com/calc/analyzer/home/analyze",

}

def analyze(sequence, nucleotide_type, oligo_conc, na_conc, mg_conc, dNTPs_conc):
    '''Function corresponding with the main "analyze" function
    Arguments:
    sequence - string containing the sequence of interest
    nucleotide_type - string with valid nucleotide type, options are 'DNA' and 'RNA'
    oligo_conc - concenctration of the oligo (in uM)
    n_conc - concentration of Na+ (in uM)
    mg_conc - concentration of Mg2+ (in uM)
    dNTPs_conc - concentration of dNTPs (in uM)

    Returns a dictionary with the following values:
    '''
    request_data = {
        'sequence': sequence,
        'nucleotideType': nucleotide_type,
        'oligoConc': oligo_conc,
        'naConc': na_conc,
        'mgConc': mg_conc,
        'dNTPsConc': dNTPs_conc
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

print(analyze("", "DNA", 0.25, 0.5, 2.5, 1))