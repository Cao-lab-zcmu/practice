## new algorithm (transform principle)
######## here the low pp will be filtered
data="results/stat_classification.tsv"
stat=$(awk -F $'\t' '
{
  if(NR==1)
    {
      for(i=1; i<=NF; i++)
        {
          if($i~/^id$/)
            {
              col_id=i
            }
          if($i~/^definition$/)
            {
              col_definition=i
            }
        }
    }
  if(NR>=2)
    {
      if($col_definition!="null")
        {
          id_class[$col_id]=$col_definition
          class_num[$col_definition]+=1
        }
    }
}
END{
for(i in class_num)
  {
    printf i"\t"  class_num[i]"\t"
    for(j in id_class)
      {
        if(id_class[j]==i)
          {
            printf j"@"
          }
      }
    printf "\n"
  }
}' $data | sed 's/@$//g')
data2="results/re_neg_RT.tsv"
echo "" > results/neg_RT.tsv
data3="results/neg_RT.tsv"
compar1="Raw"
compar2="Pro"
log=10
log_to=2
awk -F ["\t"@] -v OFS=$'\t' '
BEGIN{
maxNF=0
rows=0
compar1n=0
compar2n=0
}
{
  if(NR==FNR)
    {
      rows+=1
      if(maxNF<NF)
        {
          maxNF=NF
        }
      class[FNR]=$1
      num[FNR]=$2
      for(j=3; j<=NF; j++)
        {
          id[FNR,j]=$j
        }
    }
else
  {
    if(FNR==1 && FILENAME~/re_/)
      {
        for(i=1; i<=NF; i++)
          {
            if(($i~/'$compar1'/) && $i~/area/)
              {
                colum1[i]=i
                compar1n+=1
              }
            if(($i~/'$compar2'/) && $i~/area/)
              {
                colum2[i]=i
                compar2n+=1
              }
            if(($i~/retention/))
              {
                rtcolum=i
              }
            if(($i~/m\/z/))
              {
                mzcolum=i
              }
          }
      }
    if(FNR>=2 && $0!="")
      {
        if(sum[1,$1]=="" && sum[2,$1]=="") ######## revise data 
          {
            sum[1,$1]=0
            sum[2,$1]=0
            rt[$1]=$rtcolum
            mz[$1]=$mzcolum
            for(i=1; i<=NF; i++)
              {
                if(i==colum1[i] && colum1[i]!="")
                  {
                    sum[1,$1]+=$i
                  }
                if(i==colum2[i] && colum2[i]!="")
                  {
                    sum[2,$1]+=$i
                  }
              }
          }
      }
  }
}
END{
printf "class\t"  "id\t"  "log'$log'_raw\t"  "log'$log'_pro\t" \
  "log'$log'_delta_area\t"  "pro_to_raw\t"  "log'$log_to'_pro_to_raw\t"  "variety\t"  "number\t" > "algorithm.tsv"
  printf "rt\t"  "m/z\n" >> "algorithm.tsv"
  for(i=1; i<=rows; i++)
    {
      for(j=3; j<=maxNF; j++)
        {
          if(id[i,j]!="")
            {
              raw=(sum[1,id[i,j]]/compar1n)
              pro=(sum[2,id[i,j]]/compar2n)
              log_raw=log(raw)/log('$log')
              log_pro=log(pro)/log('$log')
              delta_area=pro-raw
              if(raw!=0)
                {
                  to_raw=pro/raw
                  norm_to_raw=log(to_raw)/log('$log_to')
                }
            else
              {
                to_raw="infinity"
                norm_to_raw="infinity"
              }
            if(delta_area+0>0)
              {
                norm_delta=log(delta_area)/log('$log')
                variety="increase"
              }
          else if(delta_area<0)
            {
              norm_delta=-1*(log(-1*delta_area)/log('$log'))
              variety="decrease"
            }
        else
          {
            norm_delta=0
          }
        printf class[i]"\t"  id[i,j]"\t"  log_raw"\t"  log_pro"\t" \
          norm_delta"\t"  to_raw"\t"  norm_to_raw"\t"  variety"\t"  num[i]"\t" \
          rt[id[i,j]]"\t"  mz[id[i,j]]"\n" >> "algorithm.tsv"
        }
    }
}
}' <(echo "$stat") $data2 $data3
#######################
data="algorithm.tsv"
awk -F $'\t' -v OFS=$'\t' '
{
  if(NR==1)
    {
      print "id","rt","m/z","classification","variety","pro/raw","log10_raw","log10_pro","norm_delta" > "compound.tsv"
      for(i=1; i<=NF; i++)
        {
          if($i~/^id/)
            {
              col_id=i
            }
          if($i~/^rt/)
            {
              col_rt=i
            }
          if($i~/m\/z/)
            {
              col_mz=i
            }
          if($i~/class/)
            {
              col_class=i
            }
          if($i~/variety/)
            {
              col_variety=i
            }
          if($i~/^pro_to_raw/)
            {
              col_ratio=i
            }
          if($i~/log10_raw/)
            {
              col_log10_raw=i
            }
          if($i~/log10_pro/)
            {
              col_log10_pro=i
            }
          if($i~/delta_area/)
            {
              col_delta=i
            }
        }
    }
  if(NR>=2)
    {
      print $col_id, $col_rt, sprintf("%.5f",$col_mz), $col_class, $col_variety, $col_ratio, $col_log10_raw, $col_log10_pro, \
        $col_delta \
        >> "compound.tsv"
      }
  }' $data
########################
data="results/fingerid_first_score.tsv"
list="compound.tsv"
awk -F $'\t' '	
{
  if(NR==FNR)
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
                  col_simi=i
                }
              if($i~/name/)
                {
                  col_name=i
                }
              if($i~/formula/)
                {
                  col_formula=i
                }
              if($i~/^inchi$/)
                {
                  col_inchi=i
                }
              if($i~/smiles/)
                {
                  col_smiles=i
                }
              if($i~/score/)
                {
                  col_score=i
                }
              if($i~/xlogp/)
                {
                  col_xlogp=i
                }
              if($i~/inchikey2D/)
                {
                  col_inchikey2D=i
                }
              if($i~/links/)
                {
                  col_links=i
                }
            }
        }
      if(FNR>=2)
        {
          simi[$col_id]=$col_simi
          name[$col_id]=$col_name
          formula[$col_id]=$col_formula
          smiles[$col_id]=$col_smiles
          inchi[$col_id]=$col_inchi
          inchikey2D[$col_id]=$col_inchikey2D
          score[$col_id]=$col_score
          xlogp[$col_id]=$col_xlogp
          links[$col_id]=$col_links
        }
    }
  if(NR!=FNR)
    {
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i~/^id/)
                {
                  col_list_id=i
                }
              printf $i"\t" > "com_compound.tsv"
            }
          printf "similarity\t"  "name\t"  "formula\t"  "xlogp\t"  "smiles\t" \
            "inchi\t"  "inchikey2D\t"  "links\n" >> "com_compound.tsv"
          }
        if(FNR>=2)
          {
            for(i=1; i<=NF; i++)
              {
                printf $i"\t" >> "com_compound.tsv"
              }
            printf simi[$col_list_id]"\t"  name[$col_list_id]"\t"  formula[$col_list_id]"\t"  xlogp[$col_list_id]"\t" \
              smiles[$col_list_id]"\t" inchi[$col_list_id]"\t"  inchikey2D[$col_list_id]"\t" \
              links[$col_list_id]"\n" >> "com_compound.tsv"
            }
        }
    }' $data $list
  #############################
  #search 
  data="com_compound.tsv"
  savename="com_lignans_and_iridoids.tsv"
  awk -F $'\t' '
  {
    if(FNR==1)
      {
        for(i=1; i<=NF; i++)
          {
            if($i~/class/)
              {
                col_classification=i
              }
          }
        printf $0"\n" > "'$savename'"
      }
    if(FNR>=2)
      {
        if(tolower($col_classification)~/iridoid/)
          {
            printf $0"\n" >> "'$savename'"
          }
        if(tolower($col_classification)~/lignan/)
          {
            printf $0"\n" >> "'$savename'"
          }
      }
  }' $data
### filter_class
###########################
###########################
###################

 #########################
 ###sun.tsv
 ###matrix
 #	 id	classification	log10_raw	log10_pro
 #	147	Lignans, neolignans and related compounds	6.11895	6.51289
 #	147	Furofuran lignans	6.11895	6.51289
 #	147	O-methylated flavonoids	6.11895	6.51289
 #	147	Coumaric acids and derivatives	6.11895	6.51289
 #	147	Amino acids and derivatives	6.11895	6.51289
 #	147	Terpene glycosides	6.11895	6.51289
 data="for_sun.tsv"
 awk -F $'\t' '
 {
   if(NR==1)
     {
       for(i=1; i<=NF; i++)
         {
           if($i~/^id/)
             {
               col_id=i
             }
           if($i~/class/)
             {
               col_class=i
             }
           if($i~/log10_raw/)
             {
               col_log_raw=i
             }
           if($i~/log10_pro/)
             {
               col_log_pro=i
             }
         }
     }
   if(NR>=2)
     {
       id[FNR]=$col_id
       class[FNR]=$col_class
       uniq_class[$col_class]=$col_class
       log_raw[FNR]=$col_log_raw
       log_pro[FNR]=$col_log_pro
     }
 }
END{
printf "name\t"  "group\t"
for(i=2; i<=NR; i++)
  {
    if(i<NR)
      {
        printf "raw_"id[i]"\t"  "pro_"id[i]"\t"
      }
    if(i==NR)
      {
        printf "raw_"id[i]"\t"  "pro_"id[i]"\n"
      }
  }
for(i in uniq_class)
  {
    printf uniq_class[i]"\t"  "NA\t"
    {
      for(j=2; j<=NR; j++)
        {
          if(j<NR)
            {
              if(class[j]==uniq_class[i])
                {
                  if(log_raw[j]!="-inf")
                    {
                      printf log_raw[j]*(-1)"\t"  log_pro[j]"\t"
                    }
                else
                  {
                    printf 0"\t"  log_pro[j]"\t"
                  }
              }
          else
            {
              printf "0\t"  "0\t"
            }
        }
    else
      {
        if(class[j]==uniq_class[i])
          {
            if(log_raw[j]!="-inf")
              {
                printf log_raw[j]*(-1)"\t"  log_pro[j]"\n"
              }
          else
            {
              printf 0"\t"  log_pro[j]"\n"
            }
        }
    else
      {
        printf "0\t"  "0\n"
      }
  }
}
}
}
}
' $data > sun.tsv
##########
sort -t $'\t' -k 7 -n com_lignans_and_iridoids.tsv > test.tsv
awk -F $'\t' '
{
  if(NR==1)
    {
      printf "rank\t"  "log10_pro\t"  $0"\n"
      rank=255
      for(i=1; i<=NF; i++)
        {
          if($i~/pro\/raw/)
            {
              col_ratio=i
            }
          if($i~/log10_raw/)
            {
              col_log_raw=i
            }
        }
    }
  if(NR>=2)
    {
      printf rank"\t"  $col_log_raw+log($col_log_raw)/log(10)"\t"  $0"\n"
      rank-=1
    }
}' test.tsv > rank.tsv
######################################
##xcms
data="com_lignans_and_iridoids.tsv"
savepath="../thermo_mzML_0518/EIC_metadata.tsv"
awk -F $'\t' '
{
  if(NR==1)
    {
      for(i=i; i<=NF; i++)
        {
          if($i~/^id$/)
            {
              col_id=i
            }
          if($i~/m\/z/)
            {
              col_mz=i
            }
        }
      printf "id\t"  "m/z\n" > "'$savepath'"
    }
  if(NR>=2)
    {
      printf $col_id"\t"  $col_mz"\n" >> "'$savepath'"
    }
}' $data
######################################
######################################
#peak during time Correction
data1="results/neg_RT.tsv"
data2="results/0924_neg_RT.tsv"
savepath="results/re_neg_RT.tsv"
mz_tolerance=0.005
rt_tolerance=0.1
awk -F $'\t' '
{
  if(NR==FNR)
    {
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i~/ID/)
                {
                  col_id=i
                }
              if($i~/retention/)
                {
                  col_rt=i
                }
              if($i~/m\/z/)
                {
                  col_mz=i
                }
            }
        }
      if(FNR>=2)
        {
          data1_mz[$col_id]=$col_mz
          data1_rt[$col_id]=$col_rt
          set[$col_id]=$col_id
          dataset[$col_id]=$0
        }
    }
  if(NR>FNR)
    {
      if(FNR==1)
        {
          printf $0"\n" > "'$savepath'"
          for(i=1; i<=NF; i++)
            {
              if($i~/ID/)
                {
                  col_id=i
                }
              if($i~/retention/)
                {
                  col_rt=i
                }
              if($i~/m\/z/)
                {
                  col_mz=i
                }
            }
        }
      if(FNR>=2)
        {
          data2_mz[$col_id]=$col_mz
          data2_rt[$col_id]=$col_rt
          for(i in data1_mz)
            {
              if(data1_mz[i]<=$col_mz+"'$mz_tolerance'" && data1_mz[i]>=$col_mz-"'$mz_tolerance'")
                {
                  if(data1_rt[i]<=$col_rt+"'$rt_tolerance'" && data1_rt[i]>=$col_rt-"'$rt_tolerance'")
                    {
                      data1_num[i]+=1
                      delete set[i];
                      print "data1",i,">>>","data2",$col_id,">>>",data1_num[i]
                      printf i"\t" data1_mz[i]"\t" data1_rt[i] >> "'$savepath'"
                      for(j=4; j<=NF; j++)
                        {
                          printf "\t"$j >> "'$savepath'"
                        }
                      printf "\n" >> "'$savepath'"
                    }
                }
            }
        }
    }
}
END{
for(i in set)
  {
    printf dataset[i]"\n" >> "'$savepath'"
    printf set[i]"\n"
  }
}' $data1 $data2
######################################
datapath="/media/wizard/back/thermo_mzML_0518/EIC"
mkdir $datapath/EIC_merge
echo "" > $datapath/file.tsv
data1="$datapath/../metadata.tsv"
data2="$datapath/EIC*.mzML/*.tsv"
awk -F $'\t' '
{
  if(NR==FNR)
    {
      if(NR>=2)
        {
          total_id[FNR]=$1
        }
    }
  if(FILENAME~/intensity/)
    {
      if(FNR==1)
        {
          n=split(FILENAME,f,"[/]")
          split(f[n-1], g,"[_]")
          samplename=g[2]
          split(f[n], a, "[_]")
          id=a[1]
          if(samplename!=p_samplename)
            {
              num_sample+=1
              sample[num_sample]=samplename
            }
          p_samplename=samplename
        }
      if(FNR>=1)
        {
          data_scan[samplename,id,FNR]=$1
          if($2!="NA")
            {
              data_intensity[samplename,id,FNR]=$2
            }
        else
          {
            data_intensity[samplename,id,FNR]="0"
          }
      }
  }
if(FILENAME~/rt.tsv/)
  {
    if(FNR==1)
      {
        n=split(FILENAME,f,"[/]")
        split(f[n-1], g,"[_]")
        samplename=g[2]
        printf samplename"\n"
      }
    if(FNR>=1)
      {
        data_scan[samplename,FNR]=$1
        data_rt[samplename,FNR]=$2
        max_rows[samplename]=FNR
      }
  }
}
END{
for(i in total_id)
  {
    printf "rt\t"  "intensity\t"  "sample\n" > "'$datapath'/EIC_merge/" total_id[i] ".tsv"
    for(j in sample)
      {
        for(k=1; k<=max_rows[sample[j]]; k++)
          {
            printf data_rt[sample[j],k]"\t"  data_intensity[sample[j],total_id[i],k]"\t"  sample[j]"\n" \
              >> "'$datapath'/EIC_merge/" total_id[i] ".tsv"
            }
        }
    }
}' $data1 $data2
######################################
mkdir results/EIC_rt_during
data1="/media/wizard/back/thermo_mzML_0518/EIC/metadata.tsv"
data2="results/re_neg_RT.tsv"
data3="/media/wizard/back/thermo_mzML_0518/EIC/EIC_merge/*.tsv"
savepath="results/EIC_rt_during/"
excess_time="0.1"
awk -F $'\t' '
{
  if(NR==FNR)
    {
      group_name[$1]=$2
    }
  if(FILENAME~/results/)
    {
      if(FNR==1)
        {
          p_file=FILENAME
          for(i=1; i<=NF; i++)
            {
              if($i~/ID/)
                {
                  col_id=i
                }
              if($i~/m\/z/)
                {
                  col_mz=i
                }
              if($i~/retention/)
                {
                  col_rt=i
                }
              if($i~/start$/)
                {
                  split($i,a,"[ ]")
                  samplename=a[1]
                  # print samplename
                  col_start[samplename]=i
                }
              if($i~/end$/)
                {
                  split($i,a,"[ ]")
                  samplename=a[1]
                  col_end[samplename]=i
                }
            }
        }
      if(FNR>=2)
        {
          mz[$col_id]=$col_mz
          # print $col_mz
          center_rt[$col_id]=$col_rt
          for(i in col_start)
            {
              if($col_start[i]!="0")
                {
                  rt_start[$col_id,i]=$col_start[i]
                  # print $col_start[i]
                }
            else if(reference_sample[$col_id]=="")
              {
                for(j in col_start)
                  {
                    if($col_start[j]!="0")
                      {
                        rt_start[$col_id,i]=$col_start[j]
                        reference_sample[$col_id]=j
                        break;
                      }
                  }
              }
          else
            {
              rt_start[$col_id,i]=$col_start[reference_sample[$col_id]]
            }
        }
      for(i in col_end)
        {
          if($col_end[i]!="0")
            {
              rt_end[$col_id,i]=$col_end[i]
            }
        else
          {
            rt_end[$col_id,i]=$col_end[reference_sample[$col_id]]
          }
      }
  }
}
if(FILENAME~/EIC_merge/)
  {
    if(FNR==1)
      {
        close(p_file)
        p_file=FILENAME
        close("'$savepath'" p_id ".tsv")
        num_id+=1
        n=split(FILENAME,a,"[/]||[.]")
        id=a[n-1]
        # print id
        p_id=id
        if(num_id==1)
          {
            for(i=1; i<=NF; i++)
              {
                if($i~/^rt/)
                  {
                    col_rt=i
                  }
                if($i~/intensity/)
                  {
                    col_intensity=i
                  }
                if($i~/sample/)
                  {
                    col_sample=i
                  }
              }
          }
        printf $0"\t" "group\t" "label\t" "color\t" "mz\t" "center_rt\n" > "'$savepath'" id ".tsv"
      }
    if(FNR>=2)
      {
        rt_min=($col_rt)/60
        if(threshold[id,$col_sample]=="")
          {
            threshold[id,$col_sample]=rt_start[id,$col_sample]+(rt_end[id,$col_sample]-rt_start[id,$col_sample])*(1/2)
          }
        if(rt_min>=rt_start[id,$col_sample]-"'$excess_time'" && rt_min<=rt_end[id,$col_sample]+"'$excess_time'")
          {
            if(rt_min+0>=threshold[id,$col_sample] && end_sig[id,$col_sample]!="1")
              {
                label_sig[id,$col_sample]=1
                end_sig[id,$col_sample]=1
              }
          else
            {
              label_sig[id,$col_sample]=0
            }
          if(rt_min+0>=rt_start[id,$col_sample]+0 && rt_min+0<=rt_end[id,$col_sample]+0)
            {
              color[id,$col_sample]=group_name[$col_sample]
            }
        else
          {
            color[id,$col_sample]="Non feature"
          }
        printf sprintf("%.2f",rt_min)"\t"  $col_intensity"\t"  $col_sample"\t"  group_name[$col_sample]"\t" label_sig[id,$col_sample]"\t"  color[id,$col_sample] >> "'$savepath'" id ".tsv"
        start_FNR[id]+=1
        if(start_FNR[id]=="1")
          {
            printf "\t"sprintf("%.4f",mz[id])  "\t"sprintf("%.2f",center_rt[id]) >> "'$savepath'" id ".tsv"
          }
        printf "\n" >> "'$savepath'" id ".tsv"
      }
  }
}
}' $data1 $data2 $data3
##############################
#### violin plot
mkdir multi_pp_class
###################################
data1="filter_0_class.tsv" #lignans and iridoids
data2="results/canopus_pp.tsv"
data3="com_compound.tsv"
ex_export="com_1011.tsv"
similarity_limit="0.4"
###################################
for i in 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99
do
  class_pp_limit=$i
  #class_pp_limit="0.9"
  savename="multi_pp_class/for_sun_$class_pp_limit.tsv"
  awk -F $'\t' '
  {
    if(NR==FNR)
      {
        filter_class[$1]=$2
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
                    print i,$col_class[i]
                  }
              }
          }
      }
    if(FILENAME~/'$data3'/)
      {
        if(FNR==1)
          {
            for(i=1; i<=NF; i++)
              {
                if($i~/^id/)
                  {
                    col_id=i
                  }
                if($i~/log10_raw/)
                  {
                    col_log_raw=i
                  }
                if($i~/log10_pro/)
                  {
                    col_log_pro=i
                  }
                if($i~/similarity/)
                  {
                    col_similarity=i
                  }
              }
            printf "id\t"  "classification\t"  "log10_raw\t"  "log10_pro\n" > "'$savename'"
            printf $0"\n" > "'$ex_export'"
          }
        if(FNR>=2)
          {
            if($col_similarity+0 >= "'$similarity_limit'"+0)
              {
                for(i in class_set)
                  {
                    if(i~"\034"$col_id"$")
                      {
                        printf $col_id"\t"  class_set[i]"\t"  $col_log_raw"\t"  $col_log_pro"\n" >> "'$savename'"
                        printf $0"\n" >> "'$ex_export'"
                      }
                  }
              }
          }
      }
  }' $data1 $data2 $data3
done
###################################
#### stat num
for i in 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99
do
  class_pp_limit=$i
  data="multi_pp_class/for_sun_$class_pp_limit.tsv"
  savepath="multi_pp_class/for_violin_${class_pp_limit}_pattern.tsv"
  awk -F $'\t' '
  {
    if(FNR==1)
      {
        printf $0"\n" > "'$savepath'"
      }
    num[$2]+=1
    data[FNR]=$0
    class[FNR]=$2
  }
END{
for(i in num)
  {
    print i,num[i]
    if(num[i]+0>=50)
      {
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
done;
#################
#####################################  the former is network_facet
for i in 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99
do
  class_pp_limit=$i
  mkdir results/network_facet_ladder2_$class_pp_limit
  data1="multi_pp_class/for_violin_${class_pp_limit}_pattern.tsv"
  data2="results/source_target_tree_0.4.tsv" # "results/source_target_tree_0.4.tsv"
  save_class="results/filter_child_class.tsv"
  savepath="results/network_facet_ladder2_$class_pp_limit/"
  awk -F $'\t' '
  {
    if(NR==FNR)
      {
        if(FNR==1)
          {
            for(i=1; i<=NF; i++)
              {
                if($i~/classification/)
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
            printf "" > "'$save_class'"
            for(i in class)
              {
                printf i"\n" >> "'$save_class'"
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
    printf stat_id[i]"\t" stat_id[i]"\t"  "1\t"  "0\t"  "null\t"  "null\t"  belong[i]"\n" >> "'$savepath'" belong[i] ".tsv"
  }
}' $data1 $data2
done;
#################################
### for ring plot
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
#################################
## Eucommia peak erea normalized
# accoding: com_lignans_and_iridoids.tsv 
data=com_lignans_and_iridoids.tsv
awk -F $'\t' '
{
}' $data
#### ############  instance for 3d plot
### step1 RT ~ intensity
Rscript ~/Downloads/codes/instance_3d_xcms.R
#################  json tree 
mkdir json_tree
data1="com_lignans_and_iridoids.tsv"
id=2268
formula=C17H24O10
data2="/media/wizard/back/0703_all/*_$id/trees/$formula*.json" ####### in media
cp $data2 json_tree/tree_$id.json
# id	rt	m/z	classification	variety	pro/raw
# 3918	13.2883753333333	701.22926	Terpene glycosides	increase	52.5192 
# 2529	11.5328771666667	699.24926	Lignan glycosides	increase	4.73482
# 3674	7.13257936666667	551.16124	Iridoid O-glycosides	increase	3.91387
# 3380	12.6588768333333	613.21304	Terpene glycosides	increase	3.9104
data="json_tree/tree_$id.json"
savepath="json_tree/"
awk -F "[ ][:][ ]||[,]" '
{
  if($0~/"root"/)
    {
      root=$2
    }
  if($0~/"id"/)
    {
      id[$2]=$2;
      the_id=$2
      getline;
      formula[the_id]=$2
    }
  if($0~/"source"/)
    {
      source[$2]=$2
      link[$2]+=1
      the_source=$2
      getline;
      target[the_source, link[the_source]]=$2
      getline;
      formula_edge[the_source, link[the_source]]=$2
    }
}
END{
printf "id\t" "label\n" > "'$savepath'"  "nodes_'$id'.tsv"
printf "from\t" "to\t" "label\n" > "'$savepath'"  "edges_'$id'.tsv"
#####
printf "root\t" root"\n" >> "'$savepath'"  "nodes_'$id'.tsv"
for(i in id)
  {
    printf i"\t" formula[i]"\n" >> "'$savepath'"  "nodes_'$id'.tsv"
  }
for(i in source)
  {
    for(j=1; j<=link[i]; j++)
      {
        printf i"\t"  target[i,j]"\t" formula_edge[i,j]"\n" >> "'$savepath'"  "edges_'$id'.tsv"
      }
  }
}' $data
######################
Rscript ~/Downloads/codes/json_tree.R
######################
## Image reshape
######################
metadata="canopus_neg.tsv"
#data="/media/wizard/back/0703_all/490_initial_8_neg_495/canopus/C17H24O10_[M-H]-.fpt"
savepath="canopus_parent_index.tsv"
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
          if($i~/^parent/)
            {
              col_parent=i
            }
        }
    }
  if(FNR>=2)
    {
      parent[$col_id]=$col_parent
      rows[FNR]=$col_id
      #print rows[FNR]
    }
}
END{
print "END"
printf "id\t"  "parentid\t"  "num\n" > "'$savepath'"
for(i=2; i<=FNR; i++)
  {
    root=rows[i]
    num[root]+=1
    index_id[rows[i]]=root
    while(root!="")
      {
        if(parent[root]!="")
          {
            index_id[rows[i]]=parent[root]"-"index_id[rows[i]]
            num[parent[root]]+=1
          }
        root=parent[root]
      }
  }
for(i=2; i<=FNR; i++)
  {
    printf rows[i]"\t"  index_id[rows[i]]"\t"  num[rows[i]]"\n" > "'$savepath'"
  }
}' $metadata
######################
#### 1028 violin plot
###################################
data1="results/filter_child_class.tsv"
data2="results/canopus_pp.tsv"
data3="com_compound.tsv"
ex_export="com_1011.tsv"
similarity_limit="0.4"
###################################
class_pp_limit="0.5"
savename="for_sun_$class_pp_limit.tsv"
awk -F $'\t' '
{
  if(NR==FNR)
    {
      filter_class[$1]=$2
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
              if(sprintf("%.3f",$col_class[i])+0 > "'$class_pp_limit'"+0)
                {
                  class_set[i,$col_id]=i
                  print i,$col_class[i]
                }
            }
        }
    }
  if(FILENAME~/'$data3'/)
    {
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i~/^id/)
                {
                  col_id=i
                }
              if($i~/log10_raw/)
                {
                  col_log_raw=i
                }
              if($i~/log10_pro/)
                {
                  col_log_pro=i
                }
              if($i~/similarity/)
                {
                  col_similarity=i
                }
            }
          printf "id\t"  "classification\t"  "log10_raw\t"  "log10_pro\n" > "'$savename'"
          printf $0"\n" > "'$ex_export'"
        }
      if(FNR>=2)
        {
          if($col_similarity+0 > "'$similarity_limit'"+0)
            {
              for(i in class_set)
                {
                  if(i~"\034"$col_id"$")
                    {
                      printf $col_id"\t"  class_set[i]"\t"  $col_log_raw"\t"  $col_log_pro"\n" >> "'$savename'"
                      printf $0"\n" >> "'$ex_export'"
                    }
                }
            }
        }
    }
}' $data1 $data2 $data3
###################################
#### stat num
class_pp_limit=0.5
data="for_sun_$class_pp_limit.tsv"
savepath="for_violin_${class_pp_limit}.tsv"
awk -F $'\t' '
{
  if(FNR==1)
    {
      printf $0"\n" > "'$savepath'"
    }
  num[$2]+=1
  data[FNR]=$0
  class[FNR]=$2
}
END{
for(i in num)
  {
    print i,num[i]
    if(num[i]+0>=10)
      {
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
#################
#################
data1="filter_class.csv" #lignans and iridoids
data2="results/canopus_pp.tsv"
data3="com_compound.tsv"
ex_export="results/com_lignans_and_iridoids.tsv"
similarity_limit="0.4"
class_pp_limit="0.5"
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
                  print i,$col_class[i]
                }
            }
        }
    }
  if(FILENAME~/'$data3'/)
    {
      if(FNR==1)
        {
          for(i=1; i<=NF; i++)
            {
              if($i~/^id/)
                {
                  col_id=i
                }
              if($i~/log10_raw/)
                {
                  col_log_raw=i
                }
              if($i~/log10_pro/)
                {
                  col_log_pro=i
                }
              if($i~/similarity/)
                {
                  col_similarity=i
                }
            }
          printf $0"\n" > "'$ex_export'"
        }
      if(FNR>=2)
        {
          if($col_similarity+0 >= "'$similarity_limit'")
            {
              for(i in class_set)
                {
                  if(i~"\034"$col_id"$")
                    {
                      printf $0"\n" >> "'$ex_export'"
                    }
                }
            }
        }
    }
}' $data1 $data2 $data3

