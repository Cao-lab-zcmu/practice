# location="~/operation/re_fecal_neg/kegg"”
cd ~/operation/re_fecal_neg/kegg
savepath="cnumber_cid.tsv"
awk -F "[:][ ]" '
{
  if(FNR==1)
    {
      printf "kegg number\t" "pubchem id" > "'$savepath'"
    }
  if($0~/BEGIN_COMPOUND/)
    {
      getline
      if($0~/C[0-9]*/)
        {
          cnum=$0
          next_sub=1
          print "Catch C number: ", $0
          printf "\n"cnum > "'$savepath'"
        }
    else
      {
        next_sub=0
        print "Invalid sublist."
      }
    }
  if($1~/PubChem/ && next_sub==1)
    {
      printf "\t"$2 > "'$savepath'"
    }
}' dblink.list 
## filter the blank and output the list
data="cnumber_cid.tsv"
nlimit=100
list=$(awk -F $'\t' '
{
  if($1~/C[0-9]*/)
    {
      if($2=="")
        {
          print "Escape CID of",$1
        }
      if($2!="")
        {
          gsub(/ /,"",$2)
          if(n>'$nlimit')
            {
              n=0
              print cid_set
            }
          n+=1
          if(n==1)
            {
              cid_set=$2
            }
        else
          {
            cid_set=cid_set","$2
          }
      }
  }
}
END{
print cid_set
}' $data)
## nrow
nrow=$(awk 'END{print NR}' <(echo "$list"))
##
for i in $(seq $nrow)
do
  cids=$(awk '{if(NR=='$i'){print $0}}' <(echo "$list"))
  check=$(echo "$cids" | awk '{if($0~/^[0-9]/){printf 0}else{printf 1}}')
  while [ $check == 0 ]
  do
    echo "Try catch pubchem API (${i}/${nrow})..."
    if [ -f ${i}_smiles.csv ]
    then
      check=$(awk '
      {
        if(FNR==1)
          {
            if($0~/CID/)
              {
                printf "1"
              }
          else
            {
              printf "0"
            }
        }
    }
  END{
  if(FNR==0)
    {
      printf "0"
    }
}' ${i}_smiles.csv)
    fi;
if [ $check == 0 ]
then
  curl --connect-timeout 20 --retry 100 --retry-delay 30 https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cids}/property/CanonicalSMILES,IsomericSMILES/CSV > ${i}_smiles.csv
fi;
done;
done;
## merge and reformat
awk -F "[,]" '
{
  if(NR==1)
    {
      printf $1
      for(i=2; i<=NF; i++)
        {
          printf "\t"$i
        }
      printf "\n"
    }
  if(FNR>=2)
    {
      printf $1
      for(i=2; i<=NF; i++)
        {
          printf "\t"$i
        }
      printf "\n"
    }
}' *_smiles.csv | sed 's/\"//g' | awk '{if($2!=""){print $0}}' > merge_smiles.tsv
## as sirius db
Rscript ~/Downloads/codes/sirius_db.R

