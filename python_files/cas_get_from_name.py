import sys
import csv
import pubchempy as pcp
data="cas.csv"
n=1
with open(data) as file:
 with open("cas","w") as f:
  reader = csv.reader(file,delimiter="\t")
  for i in reader:
   n+=1
   print(n,i)
   print(i,file=f)
   results=pcp.get_synonyms(i, "name")
   for name in results:
    print("BEGIN_compound",file=f)
    print(name,file=f)
    print("END_compound",file=f)
