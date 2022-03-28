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
