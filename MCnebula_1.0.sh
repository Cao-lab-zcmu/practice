## bin/bash

echo "We are all in the gutter, 
but some of us are looking at the stars.";
PS3='Please select the workflow to be executed. >>> '
select command in \
  "MCnebula_workflow" \
  "structure_extract" \
  "classification_extract_sum" \
  "classification_extract_filter" \
  "fragment_tree_network" \
  "fragment_tree_network_delta" \
  "structure_candidate_top10" \
  "(beta)" \
  "exit" 
do
  if [[ $command == "MCnebula_workflow" ]]
  then
    confirm=0
    until [[ $confirm == "yes" ]] || [[ $confirm == "no" ]]
    do
      read -p "Running all proccesses sequentially? [yes/no] >>> " confirm
    done;
    if [[ $confirm == "no" ]]
    then exit
    fi;
    default="structure_extract classification_extract_sum classification_extract_filter fragment_tree_network fragment_tree_network_delta"
    list=$(echo $default);
  else list=$( echo $command)
  fi;
  for option in $(echo $list)
  do
    case $option in
      ###################################
      ###################################
      ###################################
      ###################################
      structure_extract)
      echo "Run structure_extract."
      projectpath=0
      until [ -d $projectpath ] && [ -f $projectpath/canopus_summary.tsv ] && [ -f $projectpath/.format ]
      do
        read -p "Please input the path of the sirius project >>> " projectpath;
      done;
      cd $projectpath;
      mkdir results;
      mkdir temp;
      mkdir temp/fintemp;
      if [ -f temp/Mo_filename ]
      then rm temp/Mo_filename
      fi;
      data="*_*/fingerid/*.tsv"
      awk -F $'\t' '
      BEGIN{
      all_id=0
      max_id=-1
      p_id="null"
      printf "Info: loading the data..."
    }
  {
    if(FNR==1)
      {
        if(NR>FNR)
          {
            close(pfile)
          }
        f=split(FILENAME,a,"[/]");
        n=split(a[1],b,"[_]");
        id=b[n];
        if(id!=p_id)
          {
            all_id+=1
            if(id>max_id)
              {
                max_id=id
              }
          }
        file[id]=a[1]
        split(a[f], g, "[.]")
        formu_type[id]=g[1]
        if(all_id==1)
          {
            for(i=1; i<=NF; i++)
              {
                if(($i~/inchikey2D/))
                  {
                    col_2D=i
                  }
                if($i=="inchi")
                  {
                    col_inchi=i
                  }
                if(($i~/Formula/))
                  {
                    col_formu=i
                  }
                if($i=="score")
                  {
                    col_score=i
                  }
                if($i=="name")
                  {
                    col_name=i
                  }
                if($i=="smiles")
                  {
                    col_smiles=i
                  }
                if($i=="xlogp")
                  {
                    col_x=i
                  }
                if(($i~/imilarity/))
                  {
                    col_simi=i
                  }
                if($i~/links/)
                  {
                    col_links=i
                  }
              }
          }
        pfile=FILENAME
        printf g[1]"\t"formu_type[id]"\n"
        printf "Info: the data of " id " (ID) has input. The filename is " FILENAME ".\n"
      }
    if(FNR>=2)
      {
        row_score=$col_score
        row_simi=$col_simi
        if(max["score",id]=="" || max["score",id]+0<row_score+0)
          {
            max["score",id]=row_score
            name["score",id]=$col_name
            formula["score",id]=$col_formu
            formu_type["score",id]=formu_type[id]
            simi["score",id]=$col_simi
            smiles["score",id]=$col_smiles
            inchi["score",id]=$col_inchi
            in2D["score",id]=$col_2D
            score["score",id]=$col_score
            xlogp["score",id]=$col_x
            links["score",id]=$col_links
          }
        if(max["simi",id]=="" || max["simi",id]+0<row_simi+0)
          {
            max["simi",id]=row_simi
            name["simi",id]=$col_name
            formula["simi",id]=$col_formu
            formu_type["simi",id]=formu_type[id]
            simi["simi",id]=$col_simi
            smiles["simi",id]=$col_smiles
            inchi["simi",id]=$col_inchi
            in2D["simi",id]=$col_2D
            score["simi",id]=$col_score
            xlogp["simi",id]=$col_x
            links["simi",id]=$col_links
          }
      }
  }
END{
printf "id\t"  "name\t"  "formula\t"  "similarity\t"  "smiles\t"  "inchi\t"  "inchikey2D\t"  "score\t"  "xlogp\t"  "links\n" > "results/fingerid_sum.tsv"
for(l in file)
  {
    if(max["score",l]!="")
      {
        printf l"\t"  name["score",l]"\t"  formula["score",l]"\t"  simi["score",l]"\t"  smiles["score",l]"\t"  inchi["score",l]"\t"  in2D["score",l]"\t"  score["score",l]"\t"  xlogp["score",l]"\t"  links["score",l]"\n" >> "results/fingerid_sum.tsv"
        printf formula["score",l]"\t"  file[l]"\t"  formu_type["score",l]"\n" > "temp/Mo_filename"
      }
    if(max["simi",l]!="" && max["simi",l]!=max["score",l])
      {
        printf l"\t"  name["simi",l]"\t"  formula["simi",l]"\t"  simi["simi",l]"\t"  smiles["simi",l]"\t"  inchi["simi",l]"\t"  in2D["simi",l]"\t"  score["simi",l]"\t"  xlogp["simi",l]"\t"  links["simi",l]"\n" >> "results/fingerid_sum.tsv"
      }
  }
}' $data
data_sum="results/fingerid_sum.tsv"
awk -F $'\t' '
{
  if(NR==1)
    {
      for(i=1; i<=NF; i++)
        {
          if($i~/id/)
            {
              col_id=i
            }
          if($i~/score/)
            {
              col_score=i
            }
        }
    }
  if(NR>=2)
    {
      id=$col_id
      if(max_id+0<id || max_id=="")
        {
          max_id=id
        }
      if(score[id]+0<$col_score+0 || score[id]=="")
        {
          score[id]=$col_score
          data[id]=$0
        }
    }
}
END{
printf "id\t"  "name\t"  "formula\t"  "similarity\t"  "smiles\t"  "inchi\t"  "inchikey2D\t"  "score\t"  "xlogp\t"  "links\n" > "results/fingerid_first_score.tsv"
for(i=1; i<=max_id; i++)
  {
    if(data[i]!="")
      {
        printf data[i] "\n" >> "results/fingerid_first_score.tsv"
      }
  }
}' $data_sum
echo "structure_extract results have been successfully assembled into <results/fingerid_sum.tsv> and <results/fingerid_first_score.tsv>";
;;
###################################
###################################
###################################
###################################
classification_extract_sum)
echo "Run classification_extract_sum.";
if [ -f temp/Mo_filename ] && [ -f canopus.tsv ] && [ -f canopus_neg.tsv ] && [ -f .format ]
then
  echo "Project path acknowledged."
else
  until [ -d $projectpath ] && [ -f $projectpath/canopus.tsv ] && [ -f $projectpath/.format ]
  do
    read -p "Please input the path of the sirius project >>> " projectpath;
  done;
  cd $projectpath;
fi;
data1="temp/Mo_filename"
data2="canopus.tsv"
data3="canopus_neg.tsv"
datas=$(awk -F $'\t' '
{
  x=$2"/canopus/"$3".fpt"
  if(getline < x == 1)
    {
      printf x" "
    }
  close(x)
}' $data1)
##################
awk -F $'\t' '
{
  if(NR==FNR)
    {
      i=split($2,s,"[_]")
      the_id[FNR]=s[i]
      n=FNR
    }
  if(FILENAME ~ /canopus.tsv/)
    {
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i ~ /^name/)
                {
                  col_class=i
                }
              if($i ~ /absolute/)
                {
                  col_abindex=i
                }
            }
        }
      if(FNR>=2)
        {
          abindex[1,FNR]=$col_abindex
          indexset[$col_abindex]=$col_class
        }
    }
  if(FILENAME ~ /canopus_neg.tsv/)
    {
      p_filename=FILENAME
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i ~ /name/)
                {
                  col_class=i
                }
              if($i ~ /absolute/)
                {
                  col_abindex=i
                }
            }
        }
      if(FNR>=2)
        {
          abindex[2,FNR]=$col_abindex
          indexset[$col_abindex]=$col_class
        }
    }
  if(FILENAME ~ /.fpt/)
    {
      if(FNR==1)
        {
          close(p_filename)
          printf "Info: data_file name of " p_filename " has been input.\n"
          p_filename=FILENAME
          split(FILENAME,a,"[/]");
          m=split(a[1],b,"[_]");
          id=b[m];
          if(FILENAME ~ /\+.fpt/)
            {
              ion=1
            }
          if(FILENAME ~ /\-.fpt/)
            {
              ion=2
            }
        }
      pp[id,abindex[ion,FNR+1]]=sprintf("%.3f",$1)
    }
}
END{
printf "id" > "results/canopus_pp.tsv"
for(i in indexset)
  {
    ord+=1
    orderlist[ord]=i
    printf "\t"indexset[i] >> "results/canopus_pp.tsv"
  }
printf "\n" >> "results/canopus_pp.tsv"
for(i=1; i<=n; i++)
  {
    printf the_id[i] >> "results/canopus_pp.tsv"
    for(j=1; j<=ord; j++)
      {
        if(pp[the_id[i],orderlist[j]]=="")
          {
            pp[the_id[i],orderlist[j]]=0
          }
        printf "\t"pp[the_id[i],orderlist[j]] >> "results/canopus_pp.tsv"
      }
    printf "\n" >> "results/canopus_pp.tsv"
  }
}' $data1 $data2 $data3 $datas
echo "classiication_extract_sum have been successfully written into <results/canopus_pp.tsv>"
;;
###################################
###################################
###################################
classification_extract_filter)
echo "Run classification_extract_filter.";
if [ -f results/canopus_pp.tsv ]
then
  echo "Project path acknowledged."
else
  until [ -d $projectpath ] && [ -f $projectpath/results/canopus_pp.tsv ] && [ -f $projectpath/results/canopus_summary.tsv ]
  do
    read -p "Please input the path of the sirius project >>> " projectpath;
  done;
  cd $projectpath;
fi; 
check=0
until [[ "$check" == "yes" ]] || [[ "$check" == "no" ]]
do
  read -p "Collate the posterior probabilities of the classification of each feature? [yes/no] >>> " check;
done;
if [[ $check == "yes" ]]
then
  definition_limit=0
  until [[ "$definition_limit" > "0.01" ]] && [[ "$definition_limit" < "1" ]]
  do
    read -p "Please enter the classes posterior probabilities limition. 0.5~0.99 may work well. >>> " definition_limit;
  done;
  data1="canopus_summary.tsv"
  data2="canopus.tsv"
  data3="results/canopus_pp.tsv"
  awk -F $'\t' '
  {
    if(NR==FNR)
      {
        if(FNR==1)
          {
            p_file=FILENAME
            for(i=1; i<=NF; i++)
              {
                if($i~/name/)
                  {
                    col_id=i
                  }
                if($i~/specific/)
                  {
                    col_specific=i
                  }
                if($i~/level/)
                  {
                    col_level=i
                  }
                if($i~/subclass/)
                  {
                    col_subclass=i
                  }
                if($i~/^class/)
                  {
                    col_class=i
                  }
                if($i~/superclass/)
                  {
                    col_superclass=i
                  }
              }
          }
        if(FNR>=2)
          {
            n=split($col_id,a,"[_]")
            id=a[n]
            specific[id]=$col_specific
            level[id]=$col_level
            subclass[id]=$col_subclass
            class[id]=$col_class
            superclass[id]=$col_superclass
          }
      }
    if(FILENAME~/canopus.tsv/)
      {
        if(FNR==1)
          {
            close(p_file)
            p_file=FILENAME
            for(i=1; i<=NF; i++)
              {
                if($i~/name/)
                  {
                    col_name=i
                  }
                if($i~/description/)
                  {
                    col_description=i
                  }
              }
          }
        if(FNR>=2)
          {
            description[$col_name]=$col_description
          }
      }
    if(FILENAME~"'$data3'")
      {
        if(FNR==1)
          {
            close(p_file)
            for(i=1; i<=NF; i++)
              {
                if($i~/^id$/)
                  {
                    col_id=i
                  }
                if(i>=2)
                  {
                    col_class_name[i]=$i
                  }
              }
            printf "id\t"  "definition_source\t"  "definition\t"  "definition_pp\t"  "definition_description\t"  "specific\t"  "specific_pp\t"  "level_5\t"  "level_5_pp\t"  "subclass\t"  "subclass_pp\t"  "class\t"  "class_pp\t"  "superclass\t"  "superclass_pp\n"  > "results/stat_classification.tsv"
          }
        if(FNR>=2)
          {
            for(i=2; i<=NF; i++)
              {
                c_pp[col_class_name[i]]=sprintf("%.4f",$i)
              }
            if(level[$col_id] != "" && c_pp[level[$col_id]] >= "'$definition_limit'")
              {
                definition_source="level_5"
                definition=level[$col_id]
              }
          else if(subclass[$col_id] != "" && c_pp[subclass[$col_id]] >= "'$definition_limit'")
            {
              definition_source="subclass"
              definition=subclass[$col_id]
            }
        else if(class[$col_id] != "" && c_pp[class[$col_id]] >= "'$definition_limit'")
          {
            definition_source="class"
            definition=class[$col_id]
          }
      else if(superclass[$col_id] != "" && c_pp[superclass[$col_id]] >= "'$definition_limit'")
        {
          definition_source="superclass"
          definition=superclass[$col_id]
        }
    else
      {
        definition_source="null"
        definition="null"
        c_pp[definition]="null"
        description[definition]="null"
      }
    printf $col_id"\t"  definition_source"\t"  definition"\t"  c_pp[definition]"\t"  description[definition]"\t"  specific[$col_id]"\t"  c_pp[specific[$col_id]]"\t"  level[$col_id]"\t"  c_pp[level[$col_id]]"\t"  subclass[$col_id]"\t"  c_pp[subclass[$col_id]]"\t"  class[$col_id]"\t"  c_pp[class[$col_id]]"\t"  superclass[$col_id]"\t"  c_pp[superclass[$col_id]]"\n" >> "results/stat_classification.tsv"
  }
}
}' $data1 $data2 $data3
echo "classification_extract_filter have been successfully written into <results/stat_classification.tsv>"
fi;
;;
###################################
###################################
###################################
###################################
fragment_tree_network)
echo "Run fragment_tree_network.";
if [ -d temp ] && [ -f ftalign.tsv ] && [ -f .format ]
then
  echo "Project path acknowledged."
else
  until [ -d $projectpath ] && [ -f $projectpath/ftalign.tsv ] && [ -f $projectpath/.format ]
  do
    read -p "Please input the path of the sirius project. Make sure you have moved the fragment tree alignment results to the directory where the project is located. >>> " projectpath;
  done;
  cd $projectpath;
fi;
echo "Please enter the minimum alignment similarity(tlimit) that you want to filter."
tlimit=0
until [[ "$tlimit" > "0.01" ]] && [[ "$tlimit" < "0.99" ]]
do
  read -p "0.4-0.7 is recommended >>> " tlimit;
done;
if ! [ -d temp/ftaligntemp ]
then 
  mkdir temp/ftaligntemp
fi;
data="ftalign.tsv"
savepath="temp/ftaligntemp/tmp"
awk -F $'\t' -v OFS=$'\t' '
{
  printf "Info: NR = " NR ". FNR = " FNR ".\n"
  if(NR==FNR)
    {
      for(i=1; i<=NF; i++)
        {
          if(NR==1 && i!=1 || NR!=1 && i==1)
            {
              n=split($i,x,"[_]")
              raw[NR,i]=x[n]
            }
        else
          {
            if(NR==i)
              {
                raw_norm[NR,i]=$i
              }
          }
      }
  }
if(NR>FNR && FNR>=2)
  {
    for(i=2; i<=NF; i++)
      {
        norm1=$i/raw_norm[i,i]
        norm2=$i/raw_norm[FNR,FNR]
        norms=sprintf("%.2f", ((norm1+norm2)/2))
        if((norms+0 > '$tlimit'+0))
          {
            if(raw[FNR,1]+0>=raw[1,i]+0)
              {
                print raw[FNR,1], raw[1,i], norms > "'$savepath'"
              }
            if(raw[FNR,1]+0<raw[1,i]+0)
              {
                print raw[1,i], raw[FNR,1], norms > "'$savepath'"
              }
          }
      }
  }
}' $data $data
sort -u $savepath > temp/ftaligntemp/tmp2
savepath="temp/ftaligntemp/filter_net_$tlimit"
awk -F $'\t' '
{
  if(NR==FNR)
    {
      if($1!=$2)
        {
          replink[$1]=$1
          replink[$2]=$2
        }
    }
  if(NR!=FNR)
    {
      if($1!=$2)
        {
          printf $0"\n" > "'$savepath'" 
        } 
      if($1==$2)
        {
          if(replink[$1]=="")
            {
              printf $0"\n" > "'$savepath'"
            }
        }
    }
}' temp/ftaligntemp/tmp2 temp/ftaligntemp/tmp2
echo "fragment_tree_network have been successfully written into <temp/ftaligntemp/filter_net_$tlimit>"
;;
###################################
###################################
###################################
###################################
fragment_tree_network_delta)
if [ -d temp ] && [ -f ftalign.tsv ] && [ -f .format ]
then
  echo "Project path acknowledged."
else
  until [ -d $projectpath ] && [ -f $projectpath/ftalign.tsv ] && [ -f $projectpath/.format ]
  do
    read -p "Please input the path of the sirius project >>> " projectpath;
  done;
  cd $projectpath;
fi;
echo "The following module attempts to compute the differential fingerprints of connected clusters based on the molecular fingerprint results of sirius 4."
###################
###################
###################
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
;;
###################################
###################################
structure_candidate_top10)
echo "structure_candidate_top10 (awk version ≥ 3.1)."
projectpath=0
until [ -d $projectpath ] && [ -f $projectpath/.format ]
do
  read -p "Please input the path of the sirius project >>> " projectpath;
done;
cd $projectpath;
if [ -f results ]
then echo "Expand into target dir <results>"
else mkdir results
  mkdir temp
  mkdir temp/fintemp;
fi;
## gather all structure candidate
data="*_*/fingerid/*.tsv"
awk -F $'\t' '
BEGIN{
all_id=0
max_id=-1
p_id="null"
printf "Info: loading the data..."
}
{
  if(FNR==1)
    {
      if(NR>FNR)
        {
          close(pfile)
        }
      f=split(FILENAME,a,"[/]");
      n=split(a[1],b,"[_]");
      id=b[n];
      if(id!=p_id)
        {
          all_id+=1
          if(id>max_id)
            {
              max_id=id
            }
        }
      file[id]=a[1]
      split(a[f], g, "[.]")
      formu_type[id]=g[1]
      if(all_id==1)
        {
          printf "id\t"  "name\t"  "formula\t"  "similarity\t"  "smiles\t"  "inchi\t"  "inchikey2D\t"  "score\t"  "xlogp\t"  "links\n" > "results/fingerid_candidate_all.tsv"
          targetfile="results/fingerid_candidate_all.tsv"
          for(i=1; i<=NF; i++)
            {
              if(($i~/inchikey2D/))
                {
                  col_2D=i
                }
              if($i=="inchi")
                {
                  col_inchi=i
                }
              if(($i~/Formula/))
                {
                  col_formu=i
                }
              if($i=="score")
                {
                  col_score=i
                }
              if($i=="name")
                {
                  col_name=i
                }
              if($i=="smiles")
                {
                  col_smiles=i
                }
              if($i=="xlogp")
                {
                  col_x=i
                }
              if(($i~/imilarity/))
                {
                  col_simi=i
                }
              if($i~/links/)
                {
                  col_links=i
                }
            }
        }
      pfile=FILENAME
      printf g[1]"\t"formu_type[id]"\n"
      printf "Info: the data of " id " (ID) has input. The filename is " FILENAME ".\n"
    }
  if(FNR>=2)
    {
    ##  "id\t"  "name\t"  "formula\t"  "similarity\t"  "smiles\t"  "inchi\t"  "inchikey2D\t"  "score\t"  "xlogp\t"  "links\n" > targetfile
    printf id"\t" $col_name"\t" $col_formu"\t" $col_simi"\t" $col_smiles"\t" $col_inchi"\t" $col_2D"\t" $col_score"\t" $col_x"\t" $col_links"\n" > targetfile
    }
} ' $data
## sort the candidates to get top 10
data="results/fingerid_candidate_all.tsv"
version=0
version=$(awk 'BEGIN{a[1]=1; asort(a, b); print b[1]}')
if [ $version != 1 ]
then
  echo "Awk version ≤ 3.1. Function <asort> not available."
fi
awk -F $'\t' '
{
  if(FNR==1)
    {
      printf "id\t"  "name\t"  "formula\t"  "similarity\t"  "smiles\t"  "inchi\t"  "inchikey2D\t"  "score\t"  "xlogp\t"  "links\n" > "results/fingerid_candidate_top10.tsv"
      targetfile="results/fingerid_candidate_top10.tsv"
      for(i=1; i<=NF; i++)
        {
          if($i~/^id$/)
            {
              col_id=i
            }
          if($i~/^score/)
            {
              col_score=i
            }
        }
    }
  if(FNR>=2)
    {
      if(id!=$col_id)
        {
          ## output data top10
          if(id!="")
            {
              j=asort(score, a)
              for(i=j; i>=(j-9); i--)
                {
                  if(i==0)
                    {
                      break
                    }
                  if(data[a[i]]=="")
                    {
                      print "yes"
                    }
                  printf data[a[i]]"\n" > targetfile
                }
              delete data
              delete score
            }
          ## gather first row of the id
          n=1
          score[n]=$col_score
          if(data[$col_score]!="")
            {
              $col_score+=0.00000001
            }
          data[$col_score]=$0
        }
    else
        {
          n+=1
          score[n]=$col_score
          if(data[$col_score]!="")
            {
              $col_score+=0.00000001
            }
          data[$col_score]=$0
        }
      id=$col_id
    }
}
END{
{
  j=asort(score, a)
  for(i=j; i>=j-9; i--)
    {
      if(i==0)
        {
          break
        }
      printf data[a[i]]"\n" > targetfile
    }
}
}' $data
exit
;;
###################################
###################################
"(beta)")
exit;
;;
###################################
###################################
###################################
###################################
exit)
echo "The mystery of creation is like the darkness of night--it is great.
Delusions of knowledge are like the fog of the morning."
exit;
;;
###################################
###################################
###################################
###################################
*)
echo "error"
exit;
;;
esac;
done;#(for)
done;#(select)
