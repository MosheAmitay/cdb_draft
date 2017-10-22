# in: file with pdb codes
# out: file with only the complexes
# with macromolecule names, chains, aa length, resulution of complex, year  
#
 
        
import time       
from pypdb import *
from find_seperate_prot import findSeperate,seqSeperateProtein

def get_chains(dict):
    final=""
    str1=str(dict) 
    for i in str1:
        if i.isupper(): final+=i
    return final 

"""
def get_pdbs(f):
    buffer,final=[],[]
    for line in f:
        buffer+=line.split(",")
        for pdb in buffer:
            final.append(pdb.strip())
    			
    return final   
"""
pdbs=[]
#pdbs="5MS2 2JEL 1JWH 5N6U".split()
#print (pdbs)

file=open("complex2file.txt","r")
pdbs=file.read().split("\n")
print (pdbs)
results=open("complexes_2_proteins.txt","w")
complex2file=open("complex2file.txt","w")
counter=0
comlex2=[]
for i in range(30):
    print (counter,pdbs[i])
    all_info=get_all_info(pdbs[i])
    counter+=1
  #  if (len(all_info["polymer"][0]['chain'])==1) and len(all_info["polymer"][1]['chain'])==1:

    """if :
            comlex2.append(pdbs[i])
            """
    pdb=pdbs[i]
    title=describe_pdb(pdb)["title"].lower()
    prot1_name=all_info["polymer"][0]["polymerDescription"]["@description"]
    prot2_name=all_info["polymer"][1]["polymerDescription"]["@description"]
    prot1_chain=get_chains(all_info["polymer"][0]["chain"])
    prot2_chain=get_chains(all_info["polymer"][1]["chain"])
    if len(prot1_chain)>1 or len(prot2_chain)>1:
        continue
    prot1_length=int(all_info["polymer"][0]["@length"].strip())
    prot2_length=int(all_info["polymer"][1]["@length"].strip())
    if ("complex" in title) or ("bound" in title) and (prot1_length>=30) and (prot2_length>=30) and (prot1_name!=prot2_name):
        protein1=findSeperate(pdb,prot1_chain)
        protein2=findSeperate(pdb, prot2_chain)
        if (protein1 =='-1') or (protein2=='-1'):
            continue
        protein1Des=describe_pdb(protein1)
        protein2Des = describe_pdb(protein2)
        prot1_resolution=protein1Des['resolution']
        prot2_resolution = protein2Des['resolution']
        prot1_pub_year=protein1Des['release_date'][0:4]
        prot2_pub_year = protein2Des['release_date'][0:4]
        prot1_seq=seqSeperateProtein(pdb,prot1_chain,protein1)
        prot2_seq = seqSeperateProtein(pdb, prot2_chain,protein2)

        summary=(pdb+",  "+prot1_name+",  "+prot1_chain+", \
        "+str(prot1_length)+",  "+prot2_name+", "+prot2_chain+",  "+str(prot2_length)+"\n")
        results.write(summary)
        print(summary)



