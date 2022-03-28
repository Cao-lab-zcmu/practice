### marvin
## location ~/operation/back/0703_all
data="results/fingerid_first_score.tsv"
savepath="results/structure_2d/smiles_draw"
awk -F $'\t' '
{
  if(FNR==1)
    {
      for(i=1; i<=NF; i++)
        {
          if($i~/^id/)
            {
              col_id=i
            }
          if($i~/smiles/)
            {
              col_smiles=i
            }
        }
    }
  if(FNR>=2)
    {
      system("molconvert mol \""  $col_smiles  "\" -o '$savepath'/"  $col_id  ".mol")
      close("molconvert mol \""  $col_smiles  "\" -o '$savepath'/"  $col_id  ".mol")
      printf "Info: the ID is "$col_id". Row number:"FNR".\n"
    }
}' $data
####################################
### use obabel to convert mol format
ls results/structure_2d/smiles_draw/[0-9]*.mol | awk -F $'\t' '
{
  split($0, a, "[.][m][o][l]")
  system("obabel "  $0  " -imol -osvg -O "  a[1]  ".mol.svg")
  close("obabel "  $0  " -imol -osvg -O "  a[1]  ".mol.svg")
  print a[1]".svg"
}'
###################################
sed -i -e 's/white/transparent/g; s/stroke-width="2.0"/stroke-width="4.0"/g;' results/structure_2d/smiles_draw/[0-9]*.mol.svg
###################################
ls results/structure_2d/smiles_draw/[0-9]*.mol.svg | awk '
{
  if($0~/cairo/)
    {
      next
    }
else
  {
    system("cairosvg "  $0  " -o "  $0 ".cairo.svg")
    close("cairosvg "  $0  " -o "  $0 ".cairo.svg")
    printf "Info: the filename is "$0"\n"
  }
}'
##################################
##################################              candidate structure
#mkdir results/structure_2d/candidate
data="results/structure_2d/candidate/*_class.tsv"
id_set=$(awk -F $'\t' '
{
  if(FNR==1)
    {
      for(i=1; i<=NF; i++)
        {
          if($i~/"group_sub"/ || $i~/^group_sub$/)
            {
              col_group_sub=i
              #print i,NF
            }
          if($i~/"group"/ || $i~/^group$/)
            {
              col_group=i
            }
        }
    }
  if(FNR>=2)
    {
      id[$col_group_sub]=$col_group_sub
      id[$col_group]=$col_group
    }
}
END{
for(i in id)
  {
    n+=1
    if(n==1)
      {
        printf i
      }
  else
    {
      printf "_"i
    }
}
}' $data)
##################################
path="/media/wizard/back/0703_all"
awk '
BEGIN{
path="'$path'"
file_set="'$id_set'"
n=split(file_set, a, "[_]")
for(i=1; i<=n; i++)
  {
    "ls -d "  path  "/*_" a[i] | getline
    close("ls -d "  path  "/*_" a[i] | getline)
    name=name""$0"/structure_candidates.tsv "
  }
printf name > "tmp"
}'
#################################
savepath="results/structure_2d/candidate"
can_num=30
awk -F $'\t' '
{
  if(FNR==1)
    {
      n=split(FILENAME, a, "[/]")
      for(i=1; i<=n; i++)
        {
          if(a[i]~/[0-9](.*)_(.*)_(.*)[0-9]/)
            {
              file=a[i]
            }
        }
      n=split(file, b, "[_]")
      id=b[n]
      for(i=1; i<=NF; i++)
        {
          if($i~/smiles/)
            {
              col_smiles=i
            }
        }
    }
  if(FNR>=2 && FNR <='$can_num')
    {
      num[id]+=1
      system("molconvert mol \"" $col_smiles "\" -o '$savepath'/" id "_can_" num[id] ".mol")
      close("molconvert mol \"" $col_smiles "\" -o '$savepath'/" id "_can_" num[id] ".mol")
      printf id " >>> " num[id] "\n"
    }
}' $(cat tmp)
#############################
ls results/structure_2d/candidate/3918_*.mol | awk -F $'\t' '
{
  split($0, a, "[.][m][o][l]")
  system("obabel "  $0  " -imol -osvg -O "  a[1]  ".mol.svg")
  close("obabel "  $0  " -imol -osvg -O "  a[1]  ".mol.svg")
  print a[1]".svg"
}'
############################
sed -i -e 's/white/transparent/g' results/structure_2d/candidate/3918_*.mol.svg
###################################
ls results/structure_2d/candidate/3918_*.mol.svg | awk '
{
  if($0~/cairo/)
    {
      next
    }
else
  {
    system("cairosvg "  $0  " -o "  $0 ".cairo.svg")
    close("cairosvg "  $0  " -o "  $0 ".cairo.svg")
    printf "Info: the filename is "$0"\n"
  }
}'
###################################

