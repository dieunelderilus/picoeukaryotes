#Author: dieunel Derilus
# This python script will takes as input a two columns file containing gene ID and GI  number , just to link each gene id to its corresponding  Number
 import pandas as pd

##open the table delimited.txt file linking gene ID and GI number(two columns with gene_id and GI as header)
##this file could be obtained after the blast search of the proteome fasta file  against the non redundant protein database
id_gi=pd.read_table('ref_gi.txt')

##open the table delimited.txt file database  linking GI Number to  UniProt number
## this file could be obtained  from the  :  grep '  GI '  idmapping.dat | cut -f1,3 | awk '{print $2 "\t" $1}' | sort > GI-uniprot.tab.sorted 
## the idmapping.dat file could be downloaded here :ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
gudb=pd.read_table("GI-uniprot.tab.sorted")
gudb.columns=['GI','uniprot']

##open the table delimited.txt file database  linking UniProt number to  K knumber
## this file could be obtained  from   with  :  grep '  KO  ' idmapping.dat | cut -f1,3 | sort > uniprot-KO.tab.sorted 
ukdb=pd.read_table("uniprot-KO.tab.sorted")
ukdb.columns=['uniprot','knumber']

# create dictionary between Gi and uniprot
d1=dict(zip(gudb.GI,gudb.uniprot))
# map   GI  number to uniprot
id_gi['uniprot']=id_gi.GI.map(d1)
id_gi1=id_gi[pd.notnull(id_gi['uniprot'])]
# create dictionary between UniProt to Knumber
d2=dict(zip(ukdb.uniprot,ukdb.knumber))
# map uniprot to knumer 
id_gi3=id_gi2['knumber']=id_gi.uniprot.map(d1)
id_gi=id_gi[pd.notnull(id_gi['knumber'])]
