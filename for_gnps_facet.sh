plimit=0
until [[ "$plimit" > "0.5" ]] && [[ "$plimit" < "0.99" ]]
do
  read -p "If the minimum alignment similarity is greater than 0.3, this step will automatic calculate the discrepancy fingerprint between alignment features(in beta). Please enter the minimum posterior probability of the molecular fingerprint(plimit_1) to be filtered. 0.9-0.99 is recommended. >>> " plimit;
done;
plimit2=0
until [[ "$plimit2" > "0.1" ]] && [[ "$plimit2" < "$plimit" ]]
do
  read -p "The minimum posterior probability of the molecular fingerprint(plimit_2) to be controled. 0.1-0.5 may work well. >>> " plimit2;
done;
tlimit=$(ls temp/ftaligntemp/filter_net_* | awk -F _ '{print $NF}')
check_rep=$(awk 'END{print NR}' <(echo "$tlimit"))
if [[ "$check_rep" > "1" ]]
then
  tlimit=0
  until [ -f temp/ftaligntemp/filter_net_$tlimit ]
  do
    read -p "Plural and different fragment_tree_network files were found to exist locally. Please select a tlimit parameter you have utilized to run fragment_tree_network_delta. >>> " tlimit;
  done;
fi;
if [ -f temp/ftaligntemp/refilter_net_$tlimit ]
then rm temp/ftaligntemp/refilter_net_$tlimit
fi;
if [ -f temp/ftaligntemp/fpsample ]
then rm temp/ftaligntemp/fpsample
fi;
echo "Running delta_fingerprint computation..."
echo "Aquiring data from sirius index..."
data1="*_*/compound.info"
data2="temp/ftaligntemp/filter_net_$tlimit"
echo "Run fragment_tree_network_delta."
savepath="temp/ms_data_$tlimit"
echo "step 1: ms data"
if ! [ -f $savepath ]
then
  awk -F "[\t]||[:]||[_]" -v OFS=$'\t' '
  BEGIN{  
  printf "..."
  printf "id\t"  "m/z\t"  "rt\n" > "results/mz_and_rt.tsv"
}
{
  if(FILENAME~/compound.info/)
    {
      if(FNR==1)
        {
          printf "Info: catch >>> "FILENAME"\n"
        }
      if($1=="name")
        {
          i+=1;
          id[i]=$NF;
          n=split($NF,a,"[_]")
          printf a[n]"\t" >> "results/mz_and_rt.tsv"
        }
      if($1=="ionMass")
        {
          mz[id[i]]=$NF
          printf sprintf("%.4f",$NF)"\t" >> "results/mz_and_rt.tsv"
        }
      if($1=="ionType")
        {
          type[id[i]]=$NF
        }
      if($1=="rt")
        {
          rt[id[i]]=$2;
          printf sprintf("%.2f",$2/60)"\n" >> "results/mz_and_rt.tsv"
        }
    }
  if(FILENAME~/ftaligntemp/)
    {
      #source  target  ftalign  delta_mz  delta_rt  source_iontype  target_iontype;
      if(FNR==1)
        {
          printf "" > "'$savepath'"
        }
      print $1,$2,$3,sprintf("%.3f",mz[$2]-mz[$1]),sprintf("%.2f",rt[$2]-rt[$1]),type[$1],type[$2]  >> "'$savepath'"
    }
}' $data1 $data2
fi;
########################
echo "step 2: data path"
source_file="temp/Mo_filename"
data="temp/ms_data_$tlimit"
savepath="temp/datapath_$tlimit"
if ! [ -f $savepath ]
then
  awk -F $'\t' -v OFS=$'\t' '
  {
    if(NR==FNR)
      {
        n=split($2, a, "[_]")
        file[a[n]]=$2
        formu_type[a[n]]=$3
      }
    if(NR!=FNR)
      {
        #<path>sourceFormula  <path>targetFormula
        path1=file[$1]"/fingerprints/"formu_type[$1]".fpt"
        path2=file[$2]"/fingerprints/"formu_type[$2]".fpt"
        data[path1]=path1
        data[path2]=path2
      }
  }
END{
for(i in data)
  {
    n+=1
    print "Check file:",n
    if(getline<i==-1)
      {
        printf "Escape filename: " i "\n"
      }
  else
    {
      printf i" " > "'$savepath'"
    }
  close(i)
}
print "Sum:",n
}' $source_file $data
fi;
datapath=$(cat temp/datapath_$tlimit)
echo "step 3: fingerprint"
data_allfp="temp/data_allfp_$tlimit"
if ! [ -f $data_allfp ]
then
  awk -F $'\t' '
  BEGIN{
  n=0
}
{
  if(FNR==1)
    {
      if(n>1)
        {
          close(file)
        }
      file=FILENAME;
      print "Get fingerprints: ",FILENAME
      n+=1;
      printf FILENAME"\n"$0"\n" > "'$data_allfp'"
    }
else
  {
    print $0 > "'$data_allfp'"
  }
}' $datapath
fi;
###########################
awk -F $'\t' -v OFS=$'\t' '
{
  if(NR==1)
    {
      for(i=1; i<=NF; i++)
        {
          if($i~/absolute/)
            {
              col_index=i
            }
          if($i~/description/)
            {
              col_description=i
            }
        }
    }
  if(NR>=2)
    {
      print $col_index,$col_description
    }
}' csi_fingerid.tsv > temp/ftaligntemp/pos_fingerprint_index;
awk -F $'\t' -v OFS=$'\t' '
{
  if(NR==1)
    {
      for(i=1; i<=NF; i++)
        {
          if($i~/absolute/)
            {
              col_index=i
            }
          if($i~/description/)
            {
              col_description=i
            }
        }
    }
  if(NR>=2)
    {
      print $col_index,$col_description
    }
}' csi_fingerid_neg.tsv > temp/ftaligntemp/neg_fingerprint_index;
posindex="temp/ftaligntemp/pos_fingerprint_index"
negindex="temp/ftaligntemp/neg_fingerprint_index"
data="temp/ms_data_$tlimit"
echo "step 4: merge data"
awk -F $'\t' '
BEGIN{
file=0
x=0
count=0
f=0
posnum=0
negnum=0
}
{
  if(FNR==1)
    {
      file+=1
    }
  if(NR==FNR)
    {
      if(($1~/fingerprint/))
        {
          if(x+0>f+0)
            {
              f=x  # calculate the max index.
            }
          count+=1;  # calculate the all fingerprints file number.
          split($1,a,"[/]"); e=split(a[1],b,"[_]"); id=b[e]; #catch the id.
          x=0;
        }
    else
      {
        x+=1;
        fp[id,x]=$1;
      }
  }
if(file==2)
  {
    pos[FNR]=$1
    posnum+=1
  }
if(file==3)
  {
    neg[FNR]=$1
    negnum+=1
  }
if(file==4)
  {
    if("'$tlimit'"+0 >= 0.3)
      {
        if(fp[$1,1]!="" && fp[$2,1]!="")
          {
            n1=split($6, g, "[]]");
            n2=split($7, h, "[]]");
            if(g[n1]=="+" && h[n2]=="+")
              {
                for(x=1; x<=posnum; x++)
                  {
                    if(fp[$1,x]+0>='$plimit'+0 && fp[$2,x]+0<='$plimit2'+0)
                      {
                        data_s[FNR,x]=pos[x]
                      }
                  else if(fp[$2,x]+0>='$plimit'+0 && fp[$1,x]+0<='$plimit2'+0)
                    {
                      data_t[FNR,x]=pos[x]
                    }
                }
            }
          if(g[n1]=="-" && h[n2]=="-")
            {
              for(x=1; x<=negnum; x++)
                {
                  if(fp[$1,x]+0>='$plimit'+0 && fp[$2,x]+0<='$plimit2'+0)
                    {
                      data_s[FNR,x]=neg[x]
                    }
                else if(fp[$2,x]+0>='$plimit'+0 && fp[$1,x]+0<='$plimit2'+0)
                  {
                    data_t[FNR,x]=neg[x]
                  }
              }
          }
        if(g[n1]=="-" && h[n2]=="+" || h[n2]=="-" && g[n1]=="+")
          {
            for(i=1; i<=posnum; i++)
              {
                for(j=1; j<=negnum; j++)
                  {
                    if(pos[i]==neg[j])
                      {
                        mirror[i]=j  #if pos=i, the identical index of neg is mirror[i]
                      }
                  }
              }
          }
        if(g[n1]=="-" && h[n2]=="+") 
          {
            for(x=1; x<=f; x++)
              {
                if(fp[$1,mirror[x]]+0>='$plimit'+0 && fp[$2,x]+0<='$plimit2'+0)
                  {
                    data_s[FNR,x]=neg[mirror[x]]
                  }
              else if(fp[$2,x]+0>='$plimit'+0 && fp[$1,mirror[x]]+0<='$plimit2'+0)
                {
                  data_t[FNR,x]=neg[mirror[x]]
                }
            }
        }
      if(h[n2]=="-" && g[n1]=="+")
        {
          for(x=1; x<=f; x++)
            {
              if(fp[$1,x]+0>= '$plimit'+0 && fp[$2,mirror[x]]+0<='$plimit2'+0)
                {
                  data_s[FNR,x]=neg[mirror[x]]
                }
            else if(fp[$2,mirror[x]]+0>='$plimit'+0 && fp[$1,x]+0<='$plimit2'+0)
              {
                data_t[FNR,x]=neg[mirror[x]]
              }
          }   
      }
  }
}
if(FNR==1)
  {
    printf "source\t"  "target\t"  "ftalign_similarity\t"  "delta_m/z\t"  "source_fp_uniq\t"  "target_fp_uniq\n"
  }
else
  {
    printf $1"\t"  $2"\t"  $3"\t"  $4"\t";
    if('$tlimit'+0 >= 0.3)
      {
        if(fp[$1,1]=="" || fp[$2,1]=="")
          {
            printf "NA@NA\n"
          }
      else
        {
          printf "source:"  #the source fingerprint uniq.
          for(x=1; x<=f; x++)
            {
              if(data_s[FNR,x]!="")
                {
                  printf data_s[FNR,x]","
                }
            };
          printf "@target:"  #the target fingerprint uniq.
          for(x=1; x<=f; x++)
            {
              if(data_t[FNR,x]!="")
                {
                  printf data_t[FNR,x]","
                }
            }
          printf "\n"
        }
    }
else
  {
    printf "NA@NA\n"
  }
}
}
}' $data_allfp $posindex $negindex $data | sed -e 's/,@/\t/g; s/,$//g; s/@/\t/g'  > results/source_target_tree_$tlimit.tsv;
echo "All instances have written into <results/source_target_tree_$tlimit.tsv>."
#####################################
echo "step 5: separate child-nebula from parent-nebula."
data="results/stat_classification.tsv"
savepath="temp/filter_0_class.tsv"
awk -F $'\t' '
{
  if(FNR==1)
    {
      for(i=1; i<=NF; i++)
        {
          if($i~/^definition$/)
            {
              col_class=i
            }
        }
    }
  if(FNR>=2)
    {
      class[$col_class]=$col_class
    }
}
END{
for(i in class)
  {
    printf class[i]"\n" > "'$savepath'"
  }
}' $data
######################################
data1="temp/filter_0_class.tsv" 
data2="results/canopus_pp.tsv"
data3="results/fingerid_first_score.tsv"
until [[ "$similarity_limit" > 0 ]] && [[ "$similarity_limit" < 1 ]]
do
  read -p "Please set the Tanimoto similarity threshold contribute to child-nebula. >>> " similarity_limit;
done;
until [[ "$definition_limit" > 0.01 ]] && [[ "$definition_limit" < 1 ]]
do
  read -p "Please enter the classes posterior probabilities limition. 0.5~0.99 may work well. >>> " definition_limit;
done;
class_pp_limit=$definition_limit
savepath="temp/idenfication_filter_$class_pp_limit.tsv"
awk -F $'\t' '
{
  if(NR==FNR)
    {
      filter_class[$1]=$1
    }
  if(FILENAME~/canopus/)
    {
      if(FNR==1)
        {
          col_id=1
          for(i=2; i<=NF; i++)
            {
              for(j in filter_class)
                {
                  if(j==$i)
                    {
                      col_class[j]=i
                      print i
                    }
                }
            }
        }
      if(FNR>=2)
        {
          for(i in col_class)
            {
              if(sprintf("%.3f",$col_class[i])+0>="'$class_pp_limit'"+0)
                {
                  class_set[i,$col_id]=i
                  print "ID: ",$1,i,$col_class[i]
                }
            }
        }
    }
  if(FILENAME~/fingerid_first_score/)
    {
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i~/^id/)
                {
                  col_id=i
                }
              if($i~/similarity/)
                {
                  col_similarity=i
                }
            }
          printf "class_nebula_facet\t"  $0"\n" > "'$savepath'"
        }
      if(FNR>=2)
        {
          if($col_similarity+0 >= "'$similarity_limit'"+0)
            {
              for(i in class_set)
                {
                  if(i~"\034"$col_id"$")
                    {
                      printf class_set[i]"\t"  $0"\n" > "'$savepath'"
                    }
                }
            }
        }
    }
}' $data1 $data2 $data3
###################################
data="temp/idenfication_filter_$class_pp_limit.tsv"
savepath="temp/idenfication_filter2_${class_pp_limit}.tsv"
until [[ "$num_limit_1" -gt 0 ]]
do
  read -p "Please enter the features number threshold contribute to child-nebula (min number). >>> " num_limit_1;
done;
until [[ "$num_limit_2" -gt "$num_limit_1" ]]
do
  read -p "Please enter the features number threshold contribute to child-nebula (max number). >>> " num_limit_2;
done;
awk -F $'\t' '
{
  if(FNR==1)
    {
      printf $0"\n" > "'$savepath'"
      for(i=1; i<=NF; i++)
        {
          if($i~/class_nebula_facet/)
            {
              col_class=i
            }
        }
    }
  if(FNR>=2)
    {
      num[$col_class]+=1
      data[FNR]=$0
      class[FNR]=$col_class
    }
}
END{
for(i in num)
  {
    if(num[i]+0 >= '$num_limit_1'+0 && num[i]+0 <= '$num_limit_2')
      {
        printf "The nodes number of the child-nebula is " num[i] ".\n"
        for(j in class)
          {
            if(class[j]==i)
              {
                printf data[j]"\n" >> "'$savepath'"
              }
          }
      }
  }
}' $data
##################################
mkdir results/network_facet_$class_pp_limit
data1="temp/idenfication_filter2_${class_pp_limit}.tsv"
data2="results/source_target_tree_$tlimit.tsv" # "results/source_target_tree_0.4.tsv"
save_class="results/filter_child_class.tsv"
savepath="results/network_facet_$class_pp_limit/"
awk -F $'\t' '
{
  if(NR==FNR)
    {
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i~/class_nebula_facet/)
                {
                  col_class=i
                }
              if($i~/^id$/)
                {
                  col_id=i
                }
            }
        }
      if(FNR>=2)
        {
          class[$col_class]=$col_class
          class_id[$col_class,$col_id]=$col_id
          stat_id[$col_class,$col_id]=$col_id
          belong[$col_class,$col_id]=$col_class
        }
    }
  if(NR>FNR)
    {
      if(FNR==1)
        {
          for(i in class)
            {
              printf i"\n" > "'$save_class'"
              printf $0"\t"  "facet\n" > "'$savepath'" i ".tsv"
            }
        }
      if(FNR>=2)
        {
          for(i in class)
            {
              if( class_id[i,$1]==$1 && class_id[i,$2]==$2 )
                {
                  printf $0"\t"  i"\n" >> "'$savepath'" i ".tsv"
                  delete stat_id[i,$1]
                  delete stat_id[i,$2]
                }
            }
        }
    }
}
END{
for(i in stat_id)
  {
    ## source target similarity delta_mz fp fp class
    printf stat_id[i]"\t" stat_id[i]"\t"  "1\t"  "0\t"  "null\t"  "null\t"  belong[i]"\n"  >> "'$savepath'" belong[i] ".tsv"
  }
}' $data1 $data2
#####################
data1="canopus.tsv"
data2="results/filter_child_class.tsv"
data3="results/canopus_pp.tsv"
savepath="results/canopus_pp_filter.tsv"
awk -F $'\t' '
{
  if(FILENAME~/canopus.tsv/)
    {
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i~/absolute/)
                {
                  col_index=i
                }
              if($i~/^id/)
                {
                  col_chemid=i
                }
              if($i~/name/)
                {
                  col_name=i
                }
              if($i~/description/)
                {
                  col_des=i
                }
            }
        }
      if(FNR>2)
        {
          ab_index[$col_name]=$col_index
          chemid[$col_name]=$col_chemid
          des[$col_name]=$col_des
        }
    }
  if(FILENAME~/filter_child_class/)
    {
      class[$1]=$1
      if(FNR==1)
        {
          printf "index\t"  "chem_id\t"  "name\t"  "description\n" > "results/child_class.tsv"
        }
      printf ab_index[$1]"\t"  chemid[$1]"\t"  $1"\t"  des[$1]"\n" >> "results/child_class.tsv"
    }
  if(FILENAME~/canopus_pp/)
    {
      if(FNR==1)
        {
          printf $1 > "'$savepath'"
          for(i=2; i<=NF; i++)
            {
              if(class[$i]!="")
                {
                  n+=1
                  printf "\tC"ab_index[$i] >> "'$savepath'"
                  col_set[n]=i
                }
            }
          printf "\n" >> "'$savepath'"
        }
      if(FNR>=2)
        {
          printf $1 >> "'$savepath'"
          for(i=1; i<=n; i++)
            {
              printf "\t"$col_set[i] >> "'$savepath'"
            }
          printf "\n" >> "'$savepath'"
        }
    }
}' $data1 $data2 $data3
mv results/network_facet_$class_pp_limit results/gnps_network_facet_$class_pp_limit
