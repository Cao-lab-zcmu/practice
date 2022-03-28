import sys
import pubchempy as pcp
data=sys.argv[1]
cids=data.split("_")
for i in cids:
 print(i)
 compound=pcp.Compound.from_cid(i)
 for name in compound.synonyms:
  print(name,end=" | ")
  
 print()
