list=$(awk '{print NR}' filter_class.csv)
for i in $list 
do
  data1=$(awk '{if(NR=='$i'){printf $0}}' filter_class.csv) #lignans and iridoids
  data2="results/canopus_pp.tsv"
  data3="com_compound.tsv"
  ex_export="results/com_${i}.tsv"
  similarity_limit="0.4"
  class_pp_limit="0.5"
  awk -F $'\t' '
  {
    if(NR==FNR)
      {
        filter_class[$1]=$1
        index_class=$1
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
            printf "Index\t"$0"\n" > "'$ex_export'"
          }
        if(FNR>=2)
          {
            if($col_similarity+0 >= "'$similarity_limit'")
              {
                for(i in class_set)
                  {
                    if(i~"\034"$col_id"$")
                      {
                        printf index_class"\t"$0"\n" >> "'$ex_export'"
                      }
                  }
              }
          }
      }
  }' <(echo "$data1") $data2 $data3
done;
