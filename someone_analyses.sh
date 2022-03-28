mzmine_path="0630neg.csv"

savepath="fecal_pos_mzmine.tsv"

awk -F , '
 	{
 	if(NR==1)
 	 	{
 	 	for(i=4; i<=NF; i++)
 	 	 	{
 	 	 	col_sample[$i]=i
 	 	 	}
 	 	l=asorti(col_sample,b)
 	 	\
 	 	printf $1 "\t" $2 "\t" $3 "\t" > "'$savepath'"
 	 	\
 	 	for(i=1; i<l; i++)
 	 	 	{
 	 	 	printf b[i] "\t" >> "'$savepath'"
 	 	 	}
 	 	printf b[l] "\n" >> "'$savepath'"
 	 	}
 	if(NR>=2)
 	 	{
 	 	printf $1"\t" $2"\t" $3 "\t" >> "'$savepath'"
 	 	\
 	 	for(i=1; i<l; i++)
 	 	 	{
 	 	 	printf $col_sample[b[i]] "\t" >> "'$savepath'"
 	 	 	}
 	 	printf $col_sample[b[l]] "\n" >> "'$savepath'"
 	 	}
 	}' $mzmine_path
 	
 ####
 
 data="fecal_pos_mzmine.tsv"

savepath="fecal_pos_area.tsv"

awk -F $'\t' '
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
 	 	 	if($i~/area/)
 	 	 	 	{
 	 	 	 	split($i,a,"[ ]")
 	 	 	 	
 	 	 	 	num+=1
 	 	 	 	
 	 	 	 	max_num=num
 	 	 	 	
 	 	 	 	col_area[num]=i
 	 	 	 	
 	 	 	 	col_name[num]=a[1]
 	 	 	 	}
 	 	 	}
 	 	}
 	if(FNR>=2)
 	 	{
 	 	n+=1
 	 	
 	 	id_seq[n]=$col_id
 	 	
 	 	rt[$col_id]=$col_rt
 	 	
 	 	mz[$col_mz]=$col_mz
 	 	
 	 	for(i=1; i<=max_num; i++)
 	 	 	{
 	 	 	area[$col_id,i]=$col_area[i]
 	 	 	}
 	 	}
 	}
 	END{
 	 	printf "sample" > "'$savepath'"
 	 	
 	 	for(i=1; i<=n; i++)
 	 	 	{
 	 	 	printf "\t"id_seq[i] >> "'$savepath'"
 	 	 	}
 	 	printf "\n" >> "'$savepath'"
 	 	
 	 	for(i=1; i<=max_num; i++)
 	 	 	{
 	 	 	printf col_name[i] >> "'$savepath'"
 	 	 	
 	 	 	for(j=1; j<=n; j++)
 	 	 	 	{
 	 	 	 	printf "\t"area[id_seq[j],i] >> "'$savepath'"
 	 	 	 	}
 	 	 	printf "\n" >> "'$savepath'"
 	 	 	}
 	 	}' $data	
 
 #######
 	
data1="fecal_pos_area.tsv"

data2="metadata.tsv"

savepath="re_metadata.tsv"

awk -F $'\t' '
 	{
 	if(NR==FNR && FNR>=2)
 	 	{
 	 	sort_sample[FNR]=$1
 	 	}
 	if(NR>FNR && FNR==1)
 	 	{
 	 	printf $0"\n" > "'$savepath'"
 	 	}
 	if(NR>FNR && FNR>=2)
 	 	{
 	 	split($1,a,"[ ]")
 	 	
 	 	name=a[1]
 	 	
 	 	data[name]=name"\t"$2"\t"$3
 	 	}
 	}
 	END{
 	 	for(i=2; i<=FNR; i++)
 	 	 	{
 	 	 	printf data[sort_sample[i]]"\n" >> "'$savepath'"
 	 	 	}
 	 	}' $data1 $data2

mv $savepath "metadata.tsv"	

######

mkdir re_pca_plsda

 	 	
 	 	compare1="pro"
 	 	
 	 	compare2="raw"
 	 	
	 	data1="fecal_pos_area.tsv"

	 	data2="metadata.tsv"
	 	
	 	savepath1="re_pca_plsda/${compare1}_${compare2}_pca.tsv"

	 	savepath2="re_pca_plsda/metadata_${compare1}_${compare2}.tsv"
	 	
	 	awk -F $'\t' '
	 	 	{
	 	 	if(FNR==1)
	 	 	 	{
	 	 	 	if(FILENAME=="'$data1'")
	 	 	 	 	{
	 	 	 	 	print $0 > "'$savepath1'"
	 	 	 	 	}
	 	 	 	else
	 	 	 	 	{
	 	 	 	 	print $0 > "'$savepath2'"
	 	 	 	 	}
	 	 	 	}
	 	 	if(FNR>=2)
	 	 	 	{
	 	 	 	if($1~/^'$compare1'/)
	 	 	 	 	{
	 	 	 	 	if(FILENAME=="'$data1'")
	 	 	 	 	 	{
	 	 	 	 	 	print $0 >> "'$savepath1'"
	 	 	 	 	 	}
	 	 	 	 	else
	  	 	 	 	 	{
	 	 	 	 	 	print $0 >> "'$savepath2'"
	 	 	 	 	 	}
	 	 	 	 	}
	 	 	 	if($1~/^'$compare2'/)
	 	 	 	 	{
	 	 	 	 	if(FILENAME=="'$data1'")
	 	 	 	 	 	{
	 	 	 	 	 	print $0 >> "'$savepath1'"
	 	 	 	 	 	}
	 	 	 	 	else
	  	 	 	 	 	{
	 	 	 	 	 	print $0 >> "'$savepath2'"
	 	 	 	 	 	}
	 	 	 	 	}
	 	 	 	}
	 	 	}' $data1 $data2
	 	 	
################
################	 	 	
	 	 		 	
data1="fecal_pos_area.tsv"

data2="metadata.tsv"

compare1="blank"

compare2="std"

compare3="raw"

compare4="pro"

#data1="fecal_pos_area.tsv"

#data2="metadata.tsv"

savepath1="re_pca_plsda/multi_pca.tsv"

savepath2="re_pca_plsda/metadata_multi.tsv"

awk -F $'\t' '
 	{
 	if(FNR==1)
 	 	{
 	 	if(FILENAME=="'$data1'")
 	 	 	{
 	 	 	print $0 > "'$savepath1'"
 	 	 	}
 	 	else
 	 	 	{
 	 	 	print $0 > "'$savepath2'"
 	 	 	}
 	 	}
 	if(FNR>=2)
 	 	{
 	 	if($1~/^'$compare1'/)
 	 	 	{
 	 	 	if(FILENAME=="'$data1'")
 	 	 	 	{
 	 	 	 	print $0 >> "'$savepath1'"
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	print $0 >> "'$savepath2'"
 	 	 	 	}
 	 	 	}
 	 	if($1~/^'$compare2'/)
 	 	 	{
 	 	 	if(FILENAME=="'$data1'")
 	 	 	 	{
 	 	 	 	print $0 >> "'$savepath1'"
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	print $0 >> "'$savepath2'"
 	 	 	 	}
 	 	 	}
 	 	if($1~/^'$compare3'/)
 	 	 	{
 	 	 	if(FILENAME=="'$data1'")
 	 	 	 	{
 	 	 	 	print $0 >> "'$savepath1'"
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	print $0 >> "'$savepath2'"
 	 	 	 	}
 	 	 	}
 	 	if($1~/^'$compare4'/)
 	 	 	{
 	 	 	if(FILENAME=="'$data1'")
 	 	 	 	{
 	 	 	 	print $0 >> "'$savepath1'"
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	print $0 >> "'$savepath2'"
 	 	 	 	}
 	 	 	}
 	 	}
 	}' $data1 $data2
 	

wd="re_pca_plsda"

Rscript pca_ggbiplot.R $wd


wd="re_pca_plsda"

Rscript pca_prcomp.R $wd

mkdir re_pca_plsda/opls_da

data="re_pca_plsda/*.tsv"

savepath="re_pca_plsda/opls_da/"

awk -F $'\t' '
 	{
 	if(FNR==1)
 	 	{
 	 	n=split(FILENAME,a,"[/]")
 	 	
 	 	file=a[n]
 	 	}
 	if(file~/multi/)
 	 	{
 	 	nextfile;
 	 	}
 	if(file~/pca.tsv/ || file~/metadata/)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	printf $0"\n" > "'$savepath'"file
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	printf $0"\n" >> "'$savepath'"file
 	 	 	}
 	 	}
 	}' $data
 	
 #########
 	
 wd="re_pca_plsda/opls_da"

 Rscript ~/Downloads/opls_da.R $wd	
 
 
 
 
 
 

