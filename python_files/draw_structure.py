import sys
from rdkit import Chem
from rdkit.Chem import Draw
file=sys.argv[1]
with open(file, "r") as f:  # 打开文件
    data=f.read()
    
id_smiles=data.split("|||")
for combi in id_smiles:
 df=combi.split("_")
 id=df[0]
 smile=df[1]
 print("Draw structure. Feature name is",id,sep=" ")
 name=str(id),".svg"
 savename="".join(name)
 mol=Chem.MolFromSmiles(smile)
 Draw.MolToFile(mol,savename)
