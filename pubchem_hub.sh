data="fingerid_first_score.tsv"
## data='lignans_and_iridoids.tsv'
awk -F $'\t' '
 	{
 	if(NR==FNR)
 	 	{
	 	if(FNR==1)
	 	 	{
	 	 	for(i=1; i<=NF; i++)
	 	 	 	{
	 	 	 	if($i~/^id$/)
	 	 	 	 	{
	 	 	 	 	col_id=i
	 	 	 	 	}
	 	 	 	if($i~/smiles/)
	 	 	 	 	{
	 	 	 	 	col_smilse=i
	 	 	 	 	}
	 	 	 	if($i~/links/)
	 	 	 	 	{
	 	 	 	 	col_links=i
	 	 	 	 	}
	 	 	 	}
	 	 	}
	 	if(FNR>=2)
	 	 	{
	 	 	# print FNR," >>> ",$col_id
	 	 	if($col_links~/^PubChem/)
	 	 	 	{
	 	 	 	split($col_links, a, "[:][(]||[,][ ]||[)][;]||[)]")
	 	 	 	cid[$col_id]=a[2]
	 	 	 	# print a[2]
	 	 	 	}
	 	 	}
	 	}
	}
 	END{ 	
 	 	printf "id\t"  "cid\n" > "cid_metadata.tsv"
 	 	for(i in cid)
 	 	 	{
 	 	 	n+=1
 	 	 	printf i"\t"  cid[i]"\n" >> "cid_metadata.tsv"
 	 	 	if(n==1)
 	 	 	 	{
 	 	 	 	cid_set=cid[i]
 	 	 	 	}
 	 	 	if(n>=2)
 	 	 	 	{
 	 	 	 	cid_set=cid_set","cid[i]
 	 	 	 	}
 	 	 	}
 	 	# print cid_set
 	 	}' $data

list=$(
  awk -F $'\t' '
  {
    if(NR==2)
      {
        printf $2
      }
    if(NR>2)
      {
        n+=1
        if(n<=100)
          {
            printf ","$2
          }
      else
        {
          n=1
          printf "\n"$2
        }
      }
  }' cid_metadata.tsv
)

nrow=$(awk 'END{print NR}' <(echo "$list"))

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
  curl --connect-timeout 20 --retry 100 --retry-delay 30 https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cids}/property/IUPACName,CanonicalSMILES,IsomericSMILES/CSV > ${i}_smiles.csv
fi;
done;
done;

## run pubchem_hub_merge.R
