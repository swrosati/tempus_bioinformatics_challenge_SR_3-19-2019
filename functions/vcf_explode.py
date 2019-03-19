import pandas as pd
from pathlib import Path
import gzip
import os
import sys


#bring in vcf_explode
def readVCF(path, explode_info=True, explode_fmt=True):
    """Function for reading vcf files into pandas- removing the header.
    
    path = path to vcf file
    explode_info = bool, True explodes the info column so that VCF info ID is now a unique column.
    explode_fmt = bool, True explodes the fmt column so that VCF info ID is now a unique column.
        -IF there are multiple info files, str(<header> + "_") is appended to the beggining of the field.
    
    """
    import io
    if ".gzip" in str(path):
        with gzip.open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]

        #get df from lines
            vcfDF = pd.read_table(
                io.StringIO(str.join(os.linesep, lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str}
            ).rename(columns={'#CHROM': 'CHROM'})
        #     KEEPING THIS AROUND IF I EVER WANT TO MAKE A DF OF INDIVIDUAL VALUES:
        #     #Open again to get header lines:
            with open(path, 'r') as E:
                header = [l.strip("\n") for l in f if l.startswith('##')]
                headerDict = {}
                if explode_info ==True or explode_fmt == True:
                    for line in header:

                        if "INFO" in line or "FORMAT" in line:
                            lineSplit = line.split(";")
                            lineSplitInd = lineSplit[0].index("ID=")
                            ID = lineSplit[0][lineSplitInd+len("ID="):]
                            
    else:               
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]

        #get df from lines
            vcfDF = pd.read_table(
                io.StringIO(str.join(os.linesep, lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str}
            ).rename(columns={'#CHROM': 'CHROM'})
        #     KEEPING THIS AROUND IF I EVER WANT TO MAKE A DF OF INDIVIDUAL VALUES:
        #     #Open again to get header lines:
            with open(path, 'r') as E:
                header = [l.strip("\n") for l in f if l.startswith('##')] #list compress of items in header
                headerDict = {}
                if explode_info ==True or explode_fmt == True:
                    for line in header:

                        if "INFO" in line or "FORMAT" in line:
                            lineSplit = line.split(";")
                            lineSplitInd = lineSplit[0].index("ID=")
                            ID = lineSplit[0][lineSplitInd+len("ID="):]

       #             print(ID)
    return vcfDF

def VCFexplode(vcffile, explode_info=True, explode_fmt=True):
    """Parse VCF file into DataFrame. setting either explode_info or explode_fmt to true expands those columns into 
    sample ID specific columns with descrete terms.  
    """

    def VCFINFOcolExplodeToDict(rowString, rowKeys):
        """Return keys:values for each row in a vcfDF. """
        keys = [key for key in rowKeys if key in rowString] #get key of all values
        
        splitRowString = rowString.split(";")
        splitDict = dict((i.split("=")[0],i.split("=")[1]) for i in splitRowString)

        return splitDict 

    def VCFFMTColExplodeToDict(ColString, ColName, colKeys):
        """EXPLODE THE FORMAT COLUMN VALUES"""

        Values = ColString.split(":")
        keys = [ColName +"_" + key for key in colKeys] #get key of all values
          #get starting loc of all values

        return dict(zip(keys,Values))
    
    #read the file:
    vcfDF = readVCF(vcffile, explode_info, explode_fmt)
    
    if explode_info == True or explode_fmt == True:
        vcf = readVCF(vcffile) #read a second copy to work on


        with open(vcffile, "r") as vcftext: #open vcf to parse the col headers:
            vcftext = vcftext.readlines()
            vcfINFOColumns = [i[i.index("##INFO=<ID=")+len("##INFO=<ID="):i.index(",")] for i in vcftext if "##INFO=<ID=" in i]
            vcfFORMATColumns = [i[i.index("##FORMAT=<ID=")+len("##FORMAT=<ID="):i.index(",")] for i in vcftext if "##FORMAT=<ID=" in i]

    
    if explode_info == True:
        vcfDF = vcfDF.drop("INFO", axis=1)#chop down info
        vcf["INFO"] = vcf["INFO"].apply(lambda x: VCFINFOcolExplodeToDict(x,vcfINFOColumns)) #explode columns in info
        info = vcf["INFO"].apply(pd.Series)
        vcfDF = pd.concat([vcfDF, info], axis=1) #add to vcf df


        
#     vcfDF.drop("INFO") 
    if explode_fmt == True:
        vcfDF = vcfDF.drop(vcfDF.columns[9], axis=1)#chop format
        for column in vcf.iloc[:,9:]: #explode format columns higher than iformat columns
            vcf[column] = vcf[column].apply(lambda x: VCFFMTColExplodeToDict(x, column, vcfFORMATColumns))
            FormatLine = vcf[column].apply(pd.Series) #make pd scereise
            vcfDF = pd.concat([vcfDF, FormatLine], axis=1) # add to line

    def positiontoInt(positionval): #make the pos an int
        if positionval != "":
            return int(positionval)
        else:
            return positionval
        
    vcfDF["POS"] = vcfDF["POS"].apply(lambda x: positiontoInt(x))

    vcfDF = vcfDF.fillna("")
            
    return vcfDF