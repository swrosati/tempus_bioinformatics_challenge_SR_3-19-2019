# tempus_bioinformatics_challenge_SR_3-19-2019

The project contained in this folder is for the Tempus Bioinformatics Challenge.  Challenge was completed in JupytperLab v34.9 for maximum tracebility.

**DEPENDENCIES:**
See ./dependencies.txt

**FILE STRUCTURE:**
./challenge_docs = Files provided by tempus for coding challenge. Includes PDF and .vcf file.

./functions: Custom functions used for this project.
    -explode_df.py: credit to github user: MaxU (https://stackoverflow.com/questions/12680754/split-explode-pandas-dataframe-string-entry-to-separate-rows)
    -parse_variant_function: API get functions for this project (ExAC).  Both functions authored by Stefano Rosati.
    -VCF_explode: Combination of multiple pandas scripts, Stefano Rosati currated.  
    
./outputs:
    -/annotated_Tempus_BIIN_Challenge_SR_3-19-2019.csv: contains all ExAC annotated information for all variants in the VCF file.  
    All variants were kept in order to prevent loss of important variant information, as it is possible to have a SNP that is more
    deleterious than an INDEL at the same location.
    -/annotated_Tempus_BIIN_Challenge_most_deleterious_SR_3-19-2019.csv: Contains only one variant for 
    each site-- highest magnitude variant following: "ins" > "del" > "complex" > "mnp" > "snp"
    
**INFORMATION:**
- The python notebook contained in this file is authored by Stefano Rosati. 
