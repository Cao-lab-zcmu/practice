# vconda activate STOUT
echo "INFO: make sure you have run: conda activate STOUT"
workdir=$(dirs | sed 's/~/\/home\/'$USER'/g')
read -p "Please input file name to convert >>> " i
awk -F $'\t' '
{
  if(FNR==1)
    {
      for(i=1; i<=NF; i++)
        {
          if($i~/^smiles$/)
            {
              col_smiles=i
            }
        }
    }
  if(FNR>=2)
    {
      printf $col_smiles"\n" > "SMILES_temp"
    }
}' $i
sed -i s#\'##g SMILES_temp
file=SMILES_temp
cd /home/$USER/Downloads/codes/github_tool/SMILES-to-IUPAC-Translator
python STOUT_V_2.1.py --STI ${workdir}/$file ${workdir}/iupuc_${i} #SMILES to IUPAC
## python STOUT_V_2.1.py --STI test out #SMILES to IUPAC
cd $workdir
rm SMILES_temp
