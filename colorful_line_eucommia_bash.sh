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

