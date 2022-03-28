
 similarity_limit=0.4
 num_limit_1=30
 num_limit_2=500
 tlimit=0.4
 for i in 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99
 do
 definition_limit=$i
 data1="temp/filter_0_class.tsv" 
 data2="results/canopus_pp.tsv"
 data3="results/fingerid_first_score.tsv"
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
       printf stat_id[i]"\t" stat_id[i]"\t"  "1\t"  "0\t"  "null\t"  "null\t"  belong[i]"\n" \
       \
       >> "'$savepath'" belong[i] ".tsv"
       }
     }' $data1 $data2
 done;
