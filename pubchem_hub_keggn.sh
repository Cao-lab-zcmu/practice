mkdir synonyms
data=fingerid_candidate_top.tsv
# data=fingerid_first_score.tsv
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
          if(id==$col_id)
            {
              rankn+=1
            }
        else
          {
            rankn=1
          }
          id=$col_id
          # print FNR," >>> ",$col_id
          if($col_links~/^PubChem/)
            {
              l=split($col_links, a, "[;]")
              if(a[1]~/PubChem/)
                {
                  n=split( a[1], b, "[(]||[,][ ]||[)]" )
                  for(i=2; i<=n-1; i++)
                    {
                      cid_set[b[i]]=$col_id
                      rank_compound[$col_id, b[i]]=rankn
                    }
                }
            }
          if($col_links~/KEGG/)
            {
              l=split($col_links, a, "[;]")
              for(j=1; j<=l; j++)
                {
                  if(a[j]~/KEGG/)
                    {
                      n=split( a[j], b, "[(]||[,][ ]||[)]" )
                      for(i=2; i<=n-1; i++)
                        {
                          kegg_set[b[i]]=$col_id
                          rank_compound[$col_id, b[i],"kegg"]=rankn
                        }
                    }
                }
            }
        }
    }
}
END{ 	
printf "cid\t" "id\t" "kegg\t" "score_rank\n" > "synonyms/kegg_complement.tsv"
for(i in kegg_set)
  {
    printf "null\t"  kegg_set[i]"\t"  i"\t"  rank_compound[kegg_set[i], i, "kegg"]"\n" > "synonyms/kegg_complement.tsv"
  }
printf "id\t"  "cid\t" "score_rank\n" > "synonyms/cid_metadata.tsv"
for(i in cid_set)
  {
    if(i!="null")
      {
        n+=1
        printf cid_set[i]"\t"  i"\t" rank_compound[cid_set[i], i]"\n" >> "synonyms/cid_metadata.tsv"
        if(n==1)
          {
            cid_sets=i
          }
        if(n>=2)
          {
            cid_sets=cid_sets","i
          }
      }
  }
# print cid_set
}' $data
## oganize cid dataset
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
  }' synonyms/cid_metadata.tsv
)
nrow=$(awk 'END{print NR}' <(echo "$list"))
## search pubchem
for i in $(seq $nrow)
do
  cids=$(awk '{if(NR=='$i'){print $0}}' <(echo "$list"))
  check=$(echo "$cids" | awk '{if($0~/^[0-9]/){printf 0}else{printf 1}}')
  while [ $check == 0 ]
  do
    echo "Try catch pubchem API (${i}/${nrow})..."
    ## check file
    if [ -f synonyms/${i}_synonyms.XML ]
    then
      check=$(awk '{if(FNR==1){if($0~/xml version/){printf "1"}else{printf "0"}}}
      END{if(FNR==0){printf "0"}}' synonyms/${i}_synonyms.XML)
    fi;
    ## catch data via curl
    if [ $check == 0 ]
    then
      curl --connect-timeout 20 --retry 100 --retry-delay 30 https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cids}/synonyms/XML > synonyms/${i}_synonyms.XML
    fi;
  done
done
## data reformat
awk -F "[>]||[<][/]" '
{
  if(NR==1)
    {
      target_file="synonyms/reformat_kegg.tsv"
      printf "cid\t" "kegg\n" > target_file
    }
  if($1~/<CID$/)
    {
      cid=$2
    }
  if($2~/^C[0-9]*[0-9]$/)
    {
      printf cid"\t" $2"\n" > target_file
    }
}' synonyms/*_synonyms.XML

## excute next step
## history
## serum_neg.R
