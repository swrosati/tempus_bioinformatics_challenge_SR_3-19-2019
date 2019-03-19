import requests
import io
import json
from pandas.io.json import json_normalize

def df_exac_variant_GET(chrom, loc, ref, alt):
    """Send API Get request to exac for variant info. 
    
    chrom = chromosome number only ("1", not "chr1"), type(str)
    loc = variant genomic location (eg. 123456), type(str)
    ref = reference base or bases ("A", or "ATTC")
    alt = altrenate allele base or bases ("A", or "ATTC")
    
    Returns json loaded as a py dictionary. 
    
    """
    url = "http://exac.hms.harvard.edu/rest/variant/"
    url = url + "-".join([str(chrom), str(loc), str(ref), str(alt)])
    urlData = requests.get(url)
    json_data = json.loads(urlData.text)

    return json_data

    
def parse_exac_json(json):
    """Parse addition, relevenat variant information from Exac variat API response.
    
    input type: json.loads(<api_req_response>.text)
    output: dict of following items:
    
    Always annotated:
    Exac variant allele frequency, RSID (if avail)
    
    If canonical transcript noted:
    Gene Symbol, variant Consequence, Transcript (ensemble), HGVS protein change, HGVS cDNA change, Sift, PoyPhen, loss of function (LOF), variant exon, strand
        
    """
    variant_dict = {} #empty dict to add values to
    
    if "variant" in json:
        if "allele_freq" in json['variant']:
            variant_dict['exac_AF'] = json['variant']["allele_freq"]
            
    if "vep_annotations" in json:
        if "RSID" in json["vep_annotations"]:
            variant_dict["RSID"] = json["vep_annotations"]["RSID"]
    
    #create a dictionary for mapping cannonical exac response to output dict value
    transcript_consequence_dict = {
    "SYMBOL":"gene_sym", 
    "Consequence":"variant_type",    
    "HGVSc":"hgvs_cdna_mut",
    "HGVSp":"hgvs_prot_mut",  
    "SIFT":"sift",    
    "PolyPhen":"polyphen",    
    "EXON":"var_exon_number",      
    "STRAND":"strand",            
    } 
    
    #want to loop through repsonse to get the cannonical transcript, then add canoncial variant_info
    if "consequence" in json and json["consequence"] != None : 
        for var_type in json["consequence"]:
#             print(var_type)
            for transcript in json["consequence"][var_type]:
                for annot_list in json["consequence"][var_type][transcript]:
                    if annot_list["CANONICAL"] == "YES": #looking for canonical transcript
                        #iterate through dict above to map json values to variant dictionary.
                        variant_dict['canonical_transcript'] = transcript # add transcript outside of loop. 
                        for tag in transcript_consequence_dict:
                            if tag in annot_list:
                                variant_dict[transcript_consequence_dict[tag]] = annot_list[tag]
    return variant_dict
