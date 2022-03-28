import pubchempy as pcp
c.synonyms
c.isomeric_smiles

import sys
import pubchempy as pcp
data=sys.argv[1]
cids=data.split("_")
for i in cids:
 compound=pcp.Compound.from_cid(i)
 for name in compound.synonyms:
  print(name,end=" | ")
  
 print()
 
 
 
 
 
 
 
 
 #i.isomeric_smiles
 #print (end,sep)
 
