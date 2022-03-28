### Similar-fragmentation-network based multidimensional data analysis
mzmine_path="../0924_neg_RT.tsv.csv"
savepath="results/0924_neg_RT.tsv"
version=$(echo | awk '{l=asort(a,b);print l}')
 	if [ $version == 0 ]
 	then
 	mzmine_data=$(awk -F , '
 	 	{
 	 	if(NR==1)
 	 	 	{
 	 	 	for(i=4; i<=NF; i++)
 	 	 	 	{
 	 	 	 	col_sample[$i]=i
 	 	 	 	}
 	 	 	l=asorti(col_sample,b)
 	 	 	\
 	 	 	printf $1 "\t" $2 "\t" $3
 	 	 	\
 	 	 	for(i=1; i<l; i++)
 	 	 	 	{
 	 	 	 	printf b[i] "\t"
 	 	 	 	}
 	 	 	printf b[l] "\n"
 	 	 	}
 	 	if(NR>=2)
 	 	 	{
 	 	 	printf $1"\t" $2"\t" $3
 	 	 	\
 	 	 	for(i=1; i<l; i++)
 	 	 	 	{
 	 	 	 	printf $col_sample[b[i]] "\t"
 	 	 	 	}
 	 	 	printf $col_sample[b[l]] "\n"
 	 	 	}
 	 	}' $mzmine_path)
 	echo "$mzmine_data" > $savepath
 	fi;
#############
climit="0.95"
numlimit="0.2"
data1="results/canopus_pp_filter.tsv"
data2="r_network/mzmine_table.tsv"
compar1="Raw"
compar2="Pro"
log=10
log_to=2
stat=$(awk -F $'\t' -v OFS=$'\t' '
 	{
 	if(NR==1)
 	 	{
 	 	for(i=2; i<=NF; i++)
 	 	 	{
 	 	 	class[i]=$i
 	 	 	num[i]=0
 	 	 	}
 	 	}
 	if(NR>=2)
 	 	{
 	 	for(i=2; i<=NF; i++)
 	 	 	{
 	 	 	if($i>='$climit')
 	 	 	 	{
 	 	 	 	num[i]+=1;
 	 	 	 	id[i,NR]=$1;
 	 	 	 	}
 	 	 	}
 	 	}
 	}
 	END{
 	for(i=2; i<=NF; i++)
 	 	{
 	 	p=num[i]/(NR-1)
 	 	\
 	 	if(p<='$numlimit') #
 	 	 	{
 	 	 	printf class[i]"\t"  num[i]"\t"
 	 	 	\
 	 	 	for(j=2; j<=NR; j++)
 	 	 	 	{
 	 	 	 	if(id[i,j]!="")
 	 	 	 	 	{
 	 	 	 	 	printf id[i,j]"@";
 	 	 	 	 	}
 	 	 	 	}
 	 	 	printf "\n"
 	 	 	}
 	 	}
 	}' $data1 | sed 's/@$//g')
awk -F ["\t"@] -v OFS=$'\t' '
 	BEGIN{
 	 	maxNF=0
 	 	\
 	 	rows=0
 	 	\
 	 	compar1n=0
 	 	\
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
 	 	\
 	 	num[FNR]=$2
 	 	\
 	 	for(j=3; j<=NF; j++)
 	 	 	{
 	 	 	id[FNR,j]=$j
 	 	 	}
 	 	}
 	else
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if(($i~/'$compar1'/))
 	 	 	 	 	{
 	 	 	 	 	colum1[i]=i
 	 	 	 	 	compar1n+=1
 	 	 	 	 	}
 	 	 	 	if(($i~/'$compar2'/))
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
 	 	if(FNR>=2)
 	 	 	{
 	 	 	sum[1,$1]=0
 	 	 	\
 	 	 	sum[2,$1]=0
 	 	 	\
 	 	 	rt[$1]=$rtcolum
 	 	 	\
 	 	 	mz[$1]=$mzcolum
 	 	 	\
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
 	END{
 	 	printf "class\t"  "id\t"  "log'$log'_raw\t"  "log'$log'_pro\t" \
 	 	\
 	 	"log'$log'_delta_area\t"  "pro_to_raw\t"  "log'$log_to'_pro_to_raw\t"  "variety\t"  "number\t" > "boxplot.tsv"
 	 	\
 	 	printf "rt\t"  "m/z\n" >> "boxplot.tsv"
 	 	\
 	 	for(i=1; i<=rows; i++)
 	 	 	{
 	 	 	for(j=3; j<=maxNF; j++)
 	 	 	 	{
 	 	 	 	if(id[i,j]!="")
 	 	 	 	 	{
 	 	 	 	 	raw=(sum[1,id[i,j]]/compar1n)
 	 	 	 	 	\
 	 	 	 	 	pro=(sum[2,id[i,j]]/compar2n)
 	 	 	 	 	\
 	 	 	 	 	log_raw=log(raw)/log('$log')
 	 	 	 	 	\
 	 	 	 	 	log_pro=log(pro)/log('$log')
 	 	 	 	 	\
 	 	 	 	 	delta_area=pro-raw
 	 	 	 	 	\
 	 	 	 	 	if(raw!=0)
 	 	 	 	 	 	{
 	 	 	 	 	 	to_raw=pro/raw
 	 	 	 	 	 	\
 	 	 	 	 	 	norm_to_raw=log(to_raw)/log('$log_to')
 	 	 	 	 	 	}
 	 	 	 	 	else
 	 	 	 	 	 	{
 	 	 	 	 	 	to_raw="infinity"
 	 	 	 	 	 	\
 	 	 	 	 	 	norm_to_raw="infinity"
 	 	 	 	 	 	}
 	 	 	 	 	if(delta_area>0)
 	 	 	 	 	 	{
 	 	 	 	 	 	norm_delta=log(delta_area)/log('$log')
 	 	 	 	 	 	\
 	 	 	 	 	 	variety="increase"
 	 	 	 	 	 	}
 	 	 	 	 	else if(delta_area<0)
 	 	 	 	 	 	{
 	 	 	 	 	 	norm_delta=-1*(log(-1*delta_area)/log('$log'))
 	 	 	 	 	 	\
 	 	 	 	 	 	variety="decrease"
 	 	 	 	 	 	}
 	 	 	 	 	else
 	 	 	 	 	 	{
 	 	 	 	 	 	norm_delta=0
 	 	 	 	 	 	}
 	 	 	 	 	printf class[i]"\t"  id[i,j]"\t"  log_raw"\t"  log_pro"\t" \
 	 	 	 	 	\
 	 	 	 	 	norm_delta"\t"  to_raw"\t"  norm_to_raw"\t"  variety"\t"  num[i]"\t" \
 	 	 	 	 	\
 	 	 	 	 	rt[id[i,j]]"\t"  mz[id[i,j]]"\n" >> "boxplot.tsv"
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}' <(echo "$stat") $data2
##############
keywords_1="lignans"
keywords_2="iridoids"
awk -F $'\t' -v OFS=$'\t' '
 	{
 	if(NR==1)
 	 	{
 	 	print "id","rt","m/z","classification","variety","pro/raw","log10_raw","norm_delta" > "lignans_and_iridoids.tsv"
 	 	\
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
 	 	 	if($i~/delta_area/)
 	 	 	 	{
 	 	 	 	col_delta=i
 	 	 	 	}
 	 	 	}
 	 	}
 	if((tolower($col_class)~/'$keywords_1'/) || (tolower($col_class)~/'$keywords_2'/))
 	 	{
 	 	print $col_id, $col_rt, $col_mz, $col_class, $col_variety, $col_ratio, $col_log10_raw, $col_delta >> "lignans_and_iridoids.tsv"
 	 	}
  	}' boxplot.tsv
data="results/fingerid_first_score.tsv"
list="lignans_and_iridoids.tsv"
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
 	 	 	\
 	 	 	name[$col_id]=$col_name
 	 	 	\
 	 	 	formula[$col_id]=$col_formula
 	 	 	\
 	 	 	smiles[$col_id]=$col_smiles
 	 	 	\
 	 	 	inchi[$col_id]=$col_inchi
 	 	 	\
 	 	 	inchikey2D[$col_id]=$col_inchikey2D
 	 	 	\
 	 	 	score[$col_id]=$col_score
 	 	 	\
 	 	 	xlogp[$col_id]=$col_xlogp
 	 	 	\
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
 	 	 	 	printf $i"\t" > "com_lignans_and_iridoids.tsv"
 	 	 	 	}
 	 	 	printf "similarity\t"  "name\t"  "formula\t"  "xlogp\t"  "smiles\t" \
 	 	 	\
 	 	 	"inchi\t"  "inchikey2D\t"  "links\n" >> "com_lignans_and_iridoids.tsv"
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	printf $i"\t" >> "com_lignans_and_iridoids.tsv"
 	 	 	 	}
 	 	 	printf simi[$col_list_id]"\t"  name[$col_list_id]"\t"  formula[$col_list_id]"\t"  xlogp[$col_list_id]"\t" \
 	 	 	\
 	 	 	smiles[$col_list_id]"\t" inchi[$col_list_id]"\t"  inchikey2D[$col_list_id]"\t" \
 	 	 	\
 	 	 	links[$col_list_id]"\n" >> "com_lignans_and_iridoids.tsv"
 	 	 	}
 	 	}
 	}' $data $list
 ##################
 #sort 
 data_o="com_lignans_and_iridoids.tsv"
 data_r=$(cat <(head -n 1 $data_o) <(sort -t $'\t' -k 6 -n -r <(sed '1d' $data_o)))
 echo "$data_r" > com_lignans_and_iridoids.tsv
 ##################
 #data (MS1) idenfication
data1="../data_chi.tsv"
data2="../data_eng.tsv"
data3="com_lignans_and_iridoids.tsv"
tolerance=0.02 # mass_window
awk -F $'\t' '
 	{
 	if(max=="" || max<FNR)
 	 	{
 	 	max=FNR
 	 	}
 	if(NR==FNR)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if($i~/Num/)
 	 	 	 	 	{
 	 	 	 	 	num=i
 	 	 	 	 	}
 	 	 	 	if($i~/Compound/)
 	 	 	 	 	{
 	 	 	 	 	compound=i
 	 	 	 	 	}
 	 	 	 	if($i~/Precursor/)
 	 	 	 	 	{
 	 	 	 	 	precursor=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	chi_num[FNR]=$num
 	 	 	\
 	 	 	chi_com[FNR]=$compound
 	 	 	\
 	 	 	chi_pre[FNR]=$precursor
 	 	 	}
 	 	}
 	if(FILENAME~/eng/)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	close("'$data1'")
 	 	 	\
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if($i~/Num/)
 	 	 	 	 	{
 	 	 	 	 	num=i
 	 	 	 	 	}
 	 	 	 	if($i~/Compound/)
 	 	 	 	 	{
 	 	 	 	 	compound=i
 	 	 	 	 	}
 	 	 	 	if($i~/Precursor/)
 	 	 	 	 	{
 	 	 	 	 	precursor=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	eng_num[FNR]=$num
 	 	 	\
 	 	 	eng_com[FNR]=$compound
 	 	 	\
 	 	 	eng_pre[FNR]=$precursor
 	 	 	}
 	 	}
 	if(FILENAME~/com/)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	close("'$data2'")
 	 	 	\
 	 	 	printf $0"\t"  "ms1_candidate\n" > "identi_1_'$data3'"
 	 	 	\
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if($i~/m\/z/)
 	 	 	 	 	{
 	 	 	 	 	mass=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	printf $0 "\t" >> "identi_1_'$data3'"
 	 	 	\
 	 	 	for(i=1; i<=max; i++)
 	 	 	 	{
 	 	 	 	if($mass>=chi_pre[i]-'$tolerance' && $mass<=chi_pre[i]+'$tolerance')
 	 	 	 	 	{
 	 	 	 	 	printf chi_num[i]  "_"  chi_com[i]  "_"  chi_pre[i]  " | " >> "identi_1_'$data3'"
 	 	 	 	 	}
 	 	 	 	if($mass>=eng_pre[i]-'$tolerance' && $mass<=eng_pre[i]+'$tolerance')
 	 	 	 	 	{
 	 	 	 	 	printf eng_num[i]  "_"  eng_com[i]  "_"  eng_pre[i]  " | " >> "identi_1_'$data3'"
 	 	 	 	 	}
 	 	 	 	}
 	 	 	printf "\n" >> "identi_1_'$data3'"
 	 	 	}
 	 	}
 	}' $data1 $data2 $data3
 ##############
 ##database idenfication 
 data1="../identi_1.tsv"
 data2="r_network/ms2_figures/*.tsv"
 data_origin="origin"
 database="../MSMS-Public-Neg-VS15.msp"
 data3="../database.msp"
 tolerance=0.01
 ms2_tolerance=0.2
 in_tolerance=20
 weight1=80
 weight2=20
 point_limit=80
 sed 's/\r//g' $database > ../database.msp
 awk -F $'\t' '
 	{
 	if(NR==FNR)
 	 	{
 	 	if(NR==1)
 	 	 	{
 	 	 	p_file=FILENAME
 	 	 	\
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if($i~/^id/)
 	 	 	 	 	{
 	 	 	 	 	col_id=i
 	 	 	 	 	}
 	 	 	 	if($i~/m\/z/)
 	 	 	 	 	{
 	 	 	 	 	col_mz=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(NR>=2)
 	 	 	{
 	 	 	pre_mz[$col_id]=$col_mz
 	 	 	\
 	 	 	printf "info: datafile1 is ready. The data is " $col_id " to " $col_mz "\n"  #info
 	 	 	}
 	 	}
 	if(FILENAME~/ms2_figures/)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	close(p_file)
 	 	 	\
 	 	 	p_file=FILENAME
 	 	 	\
 	 	 	n=split(FILENAME,a,"[/]")
 	 	 	\
 	 	 	split(a[n],b,"[.]")
 	 	 	\
 	 	 	id=b[1]
 	 	 	\
 	 	 	printf "info: datafile2 is ready. The id number is " id "\n"  #info
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	if("'$data_origin'"=="sirius")
 	 	 	 	{
 	 	 	 	if($2+0<0)
 	 	 	 	 	{
 	 	 	 	 	count[id]+=1
 	 	 	 	 	\
 	 	 	 	 	data2_mz[id,count[id]]=$1
 	 	 	 	 	\
 	 	 	 	 	data2_in[id,count[id]]=$2*(-1)
 	 	 	 	 	\
 	 	 	 	 	printf "info: The ms2 peak is " count[id] " " $1 " " data2_in[id,count[id]] "\n"  #info
 	 	 	 	 	}
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	if($2+0>0)
 	 	 	 	 	{
 	 	 	 	 	count[id]+=1
 	 	 	 	 	\
 	 	 	 	 	data2_mz[id,count[id]]=$1
 	 	 	 	 	\
 	 	 	 	 	data2_in[id,count[id]]=$2
 	 	 	 	 	\
 	 	 	 	 	printf "info: The ms2 peak is " count[id] " " $1 " " data2_in[id,count[id]] "\n"  #info
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}
 	if(FILENAME~/.msp/)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	msp_num+=1
 	 	 	\
 	 	 	FS="[:][ ]||[\t]"
 	 	 	\
 	 	 	close(p_file)
 	 	 	\
 	 	 	p_file=FILENAME
 	 	 	\
 	 	 	n=split($0,c,"[:][ ]")
 	 	 	\
 	 	 	name=c[n]
 	 	 	\
 	 	 	printf "info: datafile3 is ready. The first name is " name "\n"  #info
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	if($1~/NAME/)
 	 	 	 	{
 	 	 	 	if(belong_id[msp_num]+0>=1)
 	 	 	 	 	{
 	 	 	 	 	for(i=1; i<=database_count[msp_num]; i++)
 	 	 	 	 	 	{
 	 	 	 	 	 	mz=data3_mz[msp_num,i]
 	 	 	 	 	 	\
 	 	 	 	 	 	intensity=(data3_in[msp_num,i]/max_in[msp_num])*100
 	 	 	 	 	 	\
 	 	 	 	 	 	if(intensity+0>=5)
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	point='$weight1'
 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	full_marks[msp_num]+=point
 	 	 	 	 	 	 	}
 	 	 	 	 	 	else
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	point='$weight2'
 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	full_marks[msp_num]+=point
 	 	 	 	 	 	 	}
 	 	 	 	 	 	for(j=1; j<=belong_id[msp_num]; j++)
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	for(k=1; k<=count[subdirectory[msp_num,j]]; k++)
 	 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	 	l1=data2_mz[subdirectory[msp_num,j],k]-'$ms2_tolerance'
 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	m1=data2_mz[subdirectory[msp_num,j],k]+'$ms2_tolerance'
 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	if(mz+0>=l1+0 && mz+0<=m1+0)
 	 	 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	 	 	l2=data2_in[subdirectory[msp_num,j],k]-'$in_tolerance'
 	 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	 	m2=data2_in[subdirectory[msp_num,j],k]+'$in_tolerance'
 	 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	 	if(intensity+0>=l2+0 && intensity+0<=m2+0)
 	 	 	 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	 	 	 	assign_point[subdirectory[msp_num,j],msp_num]+=point;
 	 	 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	 	 	break;
 	 	 	 	 	 	 	 	 	 	}
 	 	 	 	 	 	 	 	 	}
 	 	 	 	 	 	 	 	}
 	 	 	 	 	 	 	}
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	msp_num+=1
 	 	 	 	\
 	 	 	 	name=$2
 	 	 	 	}
 	 	 	if($1~/PRECURSORMZ/)
 	 	 	 	{
 	 	 	 	for(i in pre_mz)
 	 	 	 	 	{
 	 	 	 	 	if($2+0>=pre_mz[i]-'$tolerance' && $2+0<=pre_mz[i]+'$tolerance')
 	 	 	 	 	 	{
 	 	 	 	 	 	printf "info: the MS1 are " $2 " vs " pre_mz[i] "\n"
 	 	 	 	 	 	data3_mz[msp_num]=$2
 	 	 	 	 	 	belong_id[msp_num]+=1
 	 	 	 	 	 	subdirectory[msp_num,belong_id[msp_num]]=i
 	 	 	 	 	 	assign[i,msp_num]=name
 	 	 	 	 	 	sep[i,msp_num]=msp_num
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	if(belong_id[msp_num]+0>=1)
 	 	 	 	{
 	 	 	 	if($1~/PRECURSORTYPE/)
 	 	 	 	 	{
 	 	 	 	 	data3_type[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_formula[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_ontology[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_inchikey[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_smiles[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_rt[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_ccs[msp_num]=$2;
 	 	 	 	 	}
 	 	 	 	if($1~/Num Peaks/)
 	 	 	 	 	{
 	 	 	 	 	data3_peaks[msp_num]=$2
 	 	 	 	 	}
 	 	 	 	if($1~/^[0-9]/)
 	 	 	 	 	{
 	 	 	 	 	database_count[msp_num]+=1
 	 	 	 	 	\
 	 	 	 	 	data3_mz[msp_num,database_count[msp_num]]=$1
 	 	 	 	 	\
 	 	 	 	 	data3_in[msp_num,database_count[msp_num]]=$2
 	 	 	 	 	\
 	 	 	 	 	if(max_in[msp_num]<$2)
 	 	 	 	 	 	{
 	 	 	 	 	 	max_in[msp_num]=$2
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}
 	if(msp_num!="" && (FILENAME~"'$data1'"))
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	FS="\t"
 	 	 	\
 	 	 	close(p_file)
 	 	 	\
 	 	 	printf $0"\t" "custom_idenfication\n" > "../identi_2_'$data_origin'.tsv"
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	printf $0 >> "../identi_2_'$data_origin'.tsv"
 	 	 	\
 	 	 	for(i in assign_point)
 	 	 	 	{
 	 	 	 	if(i~"^"$1"\034" && assign_point[i]+0>='$point_limit')
 	 	 	 	 	{
 	 	 	 	 	score[i]=assign_point[i]
 	 	 	 	 	}
 	 	 	 	}
 	 	 	for(n=asort(score,sort_score); n>=1; n--)
 	 	 	 	{
 	 	 	 	for(i in score)
 	 	 	 	 	{
 	 	 	 	 	if(score[i]==sort_score[n])
 	 	 	 	 	 	{
 	 	 	 	 	 	printf "\t" score[i] "("  full_marks[sep[i]]  ")"  " | " \
 	 	 	 	 	 	\
 	 	 	 	 	 	data3_type[sep[i]] " | "  \
 	 	 	 	 	 	\
 	 	 	 	 	 	data3_formula[sep[i]] " | " \
 	 	 	 	 	 	\
 	 	 	 	 	 	data3_ontology[sep[i]] " | " \
 	 	 	 	 	 	\
 	 	 	 	 	 	assign[i] " | " \
 	 	 	 	 	 	\
 	 	 	 	 	 	>> "../identi_2_'$data_origin'.tsv"
 	 	 	 	 	 	\
 	 	 	 	 	 	delete score[i]; break;
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	delete score
 	 	 	\
 	 	 	printf "\n" >> "../identi_2_'$data_origin'.tsv"
 	 	 	}
 	 	}
 	}' $data1 $data2 $data3 $data1
 ##################
 ##find name from pubchem
 data1="../identi_2_origin.tsv"
 awk -F $'\t' '
  	{
  	if(NR==1)
  	 	{
  	 	FS="[\t]||[;][ ]"
  	 	for(i=1; i<=NF; i++)
  	 	 	{
  	 	 	if($i~/^id/)
  	 	 	 	{
  	 	 	 	col_id=i
  	 	 	 	}
  	 	 	if($i~/links/)
  	 	 	 	{
  	 	 	 	col_links=i
  	 	 	 	}
  	 	 	}
  	 	printf "" > "id_links.tsv"
  	 	}
  	if(NR>=2)
  	 	{
  	 	printf $col_id"\t" $col_links"\n" >> "id_links.tsv"
  	 	}
  	}' $data
 data2="id_links.tsv"
 awk -F "[\t]||[;][ ]" '
 	BEGIN{
 	 	printf "" > "id_links_sep1.tsv"
 	 	}
  	{
  	for(i=1; i<=NF; i++)
  	 	{
  	 	if(i!=NF)
  	 	 	{
  	 	 	printf $i"\t" >> "id_links_sep1.tsv"
  	 	 	}
  	 	else
  	 	 	{
  	 	 	printf $i"\n" >> "id_links_sep1.tsv"
  	 	 	}
  	 	}
  	}' $data2
 data3="id_links_sep1.tsv"
 awk -F $'\t' '
 	BEGIN{
 	 	printf "" > "id_cid.tsv"
 	 	}
 	{
 	if($2~/PubChem/)
 	 	{
 	 	split($2,a,"[(]||[,]||[)]")
 	 	\
 	 	if(a[2]~/^[0-9]/)
 	 	 	{
 	 	 	printf $1 "\t" a[2] "\n" >> "id_cid.tsv"
 	 	 	}
 	 	}
 	}' $data3
 data4="id_cid.tsv"
 arg=$(awk -F $'\t' '{printf $1"_"$2"@"}' $data4 | sed 's/@$//g')
script=$(echo '
import sys
import pubchempy as pcp
data=sys.argv[1]
id_cids=data.split("@")
for i in id_cids:
 data_file=open("synonyms.tsv","a+")
 essemble=i.split("_")
 id=essemble[0]
 cid=essemble[1]
 compound=pcp.Compound.from_cid(cid)
 smiles=compound.isomeric_smiles
 print(str(id)+"\t",end="",file=data_file)
 for name in compound.synonyms:
  print(name,end=" | ",file=data_file)
 print("\t"+smiles,file=data_file)
 print("Info: the data of "+str(id)+"_"+str(cid)+" has been collected.")
 data_file.close()')
python <(echo "$script") $arg
 ################################################
 ################################################
 ## transform principle
 data="../identi_2_origin.tsv"
 tolerance="0.1"
 ratio_limit="1.2"
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
 	 	 	if($i~/norm_delta/)
 	 	 	 	{
 	 	 	 	col_delta=i
 	 	 	 	}
 	 	 	if($i~/class/)
 	 	 	 	{
 	 	 	 	col_class=i
 	 	 	 	}
 	 	 	if($i~/pro\/raw/)
 	 	 	 	{
 	 	 	 	col_ratio=i
 	 	 	 	}
 	 	 	}
 	 	}
 	if(NR>=2)
 	 	{
 	 	id[FNR]=$col_id
 	 	\
 	 	norm_delta[$col_id]=$col_delta
 	 	\
 	 	class[$col_id]=$col_class
 	 	\
 	 	ratio[$col_id]=$col_ratio
 	 	}
 	}
 	END{
 	 	printf "source\t" "target\t" "source_delta\t" "target_delta\t" "trans_probability\t" \
 	 	\
 	 	"ratio\t" "classification\n" > "principle_'$tolerance'_'$ratio_limit'.tsv"
 	 	\
 	 	for(i in id)
 	 	 	{
 	 	 	if(ratio[id[i]]+0>="'$ratio_limit'"+0)
	 	 	 	{
	 	 	 	for(j in norm_delta)
	 	 	 	 	{
	 	 	 	 	if(norm_delta[id[i]]+0>0 && norm_delta[j]+0<0 || norom_delta[id[i]]+0<0 && norm_delta[j]+0>0)
	 	 	 	 	 	{
	 	 	 	 	 	l=(norm_delta[j])*(-1)-"'$tolerance'"
	 	 	 	 	 	\
	 	 	 	 	 	m=(norm_delta[j])*(-1)+"'$tolerance'"
	 	 	 	 	 	\
	 	 	 	 	 	if(norm_delta[id[i]]+0>=l && norm_delta[id[i]]<=m && class[id[i]]==class[j])
	 	 	 	 	 	 	{
	 	 	 	 	 	 	trans=norm_delta[id[i]]+norm_delta[j]
	 	 	 	 	 	 	\
	 	 	 	 	 	 	if(trans+0<0)
	 	 	 	 	 	 	 	{
	 	 	 	 	 	 	 	trans*=(-1)
	 	 	 	 	 	 	 	}
	 	 	 	 	 	 	if(id[i]+0>j)
	 	 	 	 	 	 	 	{
	 	 	 	 	 	 	 	printf id[i]"\t" j"\t" norm_delta[id[i]]"\t" norm_delta[j]"\t" \
	 	 	 	 	 	 	 	\
	 	 	 	 	 	 	 	trans"\t" id[i]":"ratio[id[i]]" | "j":"ratio[j]"\t" \
	 	 	 	 	 	 	 	\
	 	 	 	 	 	 	 	class[id[i]]"\n" >> "principle_'$tolerance'_'$ratio_limit'.tsv"
	 	 	 	 	 	 	 	}
	 	 	 	 	 	 	if(id[i]+0<j)
	 	 	 	 	 	 	 	{
	 	 	 	 	 	 	 	printf j"\t" id[i]"\t" norm_delta[j]"\t" norm_delta[id[i]]"\t" \
	 	 	 	 	 	 	 	\
	 	 	 	 	 	 	 	trans"\t" j":"ratio[j]" | "id[i]":"ratio[id[i]]"\t" \
	 	 	 	 	 	 	 	\
	 	 	 	 	 	 	 	class[id[i]]"\n" >> "principle_'$tolerance'_'$ratio_limit'.tsv"
	 	 	 	 	 	 	 	}
	 	 	 	 	 	 	}
	 	 	 	 	 	}
	 	 	 	 	}
	 	 	 	}
 	 	 	delete norm_delta[id[i]]
 	 	 	}
 	 	}' $data
 #######################################
 #######################################
 #######################################
 ## classification arrangement
 data1="canopus_summary.tsv"
 data2="canopus.tsv"
 data3="results/canopus_pp.tsv"
 definition_limit="0.9"
 awk -F $'\t' '
 	{
 	if(NR==FNR)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	p_file=FILENAME
 	 	 	\
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
 	 	 	\
 	 	 	id=a[n]
 	 	 	\
 	 	 	specific[id]=$col_specific
 	 	 	\
 	 	 	level[id]=$col_level
 	 	 	\
 	 	 	subclass[id]=$col_subclass
 	 	 	\
 	 	 	class[id]=$col_class
 	 	 	\
 	 	 	superclass[id]=$col_superclass
 	 	 	}
 	 	}
 	if(FILENAME~/canopus.tsv/)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	close(p_file)
 	 	 	\
 	 	 	p_file=FILENAME
 	 	 	\
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
 	 	 	\
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
 	 	 	printf "id\t"  "definition_source\t"  "definition\t"  "definition_pp\t"  "definition_description\t" \
 	 	 	\
 	 	 	"specific\t"  "specific_pp\t"  "level_5\t"  "level_5_pp\t" \
 	 	 	\
 	 	 	"subclass\t"  "subclass_pp\t"  "class\t"  "class_pp\t"  "superclass\t"  "superclass_pp\n" > "stat_classification.tsv"
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	for(i=2; i<=NF; i++)
 	 	 	 	{
 	 	 	 	c_pp[col_class_name[i]]=sprintf("%.4f",$i)
 	 	 	 	}
 	 	 	if(level[$col_id] != "" && c_pp[level[$col_id]]+0 >= "'$definition_limit'"+0)
 	 	 	 	{
 	 	 	 	definition_source="level_5"
 	 	 	 	\
 	 	 	 	definition=level[$col_id]
 	 	 	 	}
 	 	 	else if(subclass[$col_id] != "" && c_pp[subclass[$col_id]]+0 >= "'$definition_limit'"+0)
 	 	 	 	{
 	 	 	 	definition_source="subclass"
 	 	 	 	\
 	 	 	 	definition=subclass[$col_id]
 	 	 	 	}
 	 	 	else if(class[$col_id] != "" && c_pp[class[$col_id]]+0 >= "'$definition_limit'"+0)
 	 	 	 	{
 	 	 	 	definition_source="class"
 	 	 	 	\
 	 	 	 	definition=class[$col_id]
 	 	 	 	}
 	 	 	else if(superclass[$col_id] != "" && c_pp[superclass[$col_id]]+0 >= "'$definition_limit'"+0)
 	 	 	 	{
 	 	 	 	definition_source="superclass"
 	 	 	 	\
 	 	 	 	definition=superclass[$col_id]
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	definition_source="null"
 	 	 	 	\
 	 	 	 	definition="null"
 	 	 	 	\
 	 	 	 	c_pp[definition]="null"
 	 	 	 	\
 	 	 	 	description[definition]="null"
 	 	 	 	}
 	 	 	printf $col_id"\t"  definition_source"\t"  definition"\t"  c_pp[definition]"\t"  description[definition]"\t" \
 	 	 	\
 	 	 	specific[$col_id]"\t"  c_pp[specific[$col_id]]"\t"  level[$col_id]"\t"  c_pp[level[$col_id]]"\t" \
 	 	 	\
 	 	 	subclass[$col_id]"\t"  c_pp[subclass[$col_id]]"\t"  class[$col_id]"\t"  c_pp[class[$col_id]]"\t" \
 	 	 	\
 	 	 	superclass[$col_id]"\t"  c_pp[superclass[$col_id]]"\n" >> "stat_classification.tsv"
 	 	 	}
 	 	}
 	}' $data1 $data2 $data3
 ##################
 ## translate
 data=$(awk -F $'\t' '
 	{
 	if(NR==1)
 	 	{
 	 	for(i=1; i<=NF; i++)
 	 	 	{
 	 	 	if($i~/^definition$/)
 	 	 	 	{
 	 	 	 	col_definition=i
 	 	 	 	}
 	 	 	}
 	 	}
 	if(NR>=2 && $col_definition!="null")
 	 	{
 	 	print $col_definition
 	 	}
 	}' stat_classification.tsv | sort | uniq)
 trans :zh -b <<< echo "$data"
 ################################################################################
 ################################################################################
 ################################################################################
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
 	 	 	\
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
 data3="results/neg_RT.tsv"
 compar1="Raw"
 compar2="Pro"
 log=10
 log_to=2
 awk -F ["\t"@] -v OFS=$'\t' '
 	BEGIN{
 	 	maxNF=0
 	 	\
 	 	rows=0
 	 	\
 	 	compar1n=0
 	 	\
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
 	 	\
 	 	num[FNR]=$2
 	 	\
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
 	 	if(FNR>=2)
 	 	 	{
 	 	 	if(sum[1,$1]=="" && sum[2,$1]=="") ######## revise data 
 	 	 	 	{
	 	 	 	sum[1,$1]=0
	 	 	 	\
	 	 	 	sum[2,$1]=0
	 	 	 	\
	 	 	 	rt[$1]=$rtcolum
	 	 	 	\
	 	 	 	mz[$1]=$mzcolum
	 	 	 	\
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
 	 	\
 	 	"log'$log'_delta_area\t"  "pro_to_raw\t"  "log'$log_to'_pro_to_raw\t"  "variety\t"  "number\t" > "algorithm.tsv"
 	 	\
 	 	printf "rt\t"  "m/z\n" >> "algorithm.tsv"
 	 	\
 	 	for(i=1; i<=rows; i++)
 	 	 	{
 	 	 	for(j=3; j<=maxNF; j++)
 	 	 	 	{
 	 	 	 	if(id[i,j]!="")
 	 	 	 	 	{
 	 	 	 	 	raw=(sum[1,id[i,j]]/compar1n)
 	 	 	 	 	\
 	 	 	 	 	pro=(sum[2,id[i,j]]/compar2n)
 	 	 	 	 	\
 	 	 	 	 	log_raw=log(raw)/log('$log')
 	 	 	 	 	\
 	 	 	 	 	log_pro=log(pro)/log('$log')
 	 	 	 	 	\
 	 	 	 	 	delta_area=pro-raw
 	 	 	 	 	\
 	 	 	 	 	if(raw!=0)
 	 	 	 	 	 	{
 	 	 	 	 	 	to_raw=pro/raw
 	 	 	 	 	 	\
 	 	 	 	 	 	norm_to_raw=log(to_raw)/log('$log_to')
 	 	 	 	 	 	}
 	 	 	 	 	else
 	 	 	 	 	 	{
 	 	 	 	 	 	to_raw="infinity"
 	 	 	 	 	 	\
 	 	 	 	 	 	norm_to_raw="infinity"
 	 	 	 	 	 	}
 	 	 	 	 	if(delta_area+0>0)
 	 	 	 	 	 	{
 	 	 	 	 	 	norm_delta=log(delta_area)/log('$log')
 	 	 	 	 	 	\
 	 	 	 	 	 	variety="increase"
 	 	 	 	 	 	}
 	 	 	 	 	else if(delta_area<0)
 	 	 	 	 	 	{
 	 	 	 	 	 	norm_delta=-1*(log(-1*delta_area)/log('$log'))
 	 	 	 	 	 	\
 	 	 	 	 	 	variety="decrease"
 	 	 	 	 	 	}
 	 	 	 	 	else
 	 	 	 	 	 	{
 	 	 	 	 	 	norm_delta=0
 	 	 	 	 	 	}
 	 	 	 	 	printf class[i]"\t"  id[i,j]"\t"  log_raw"\t"  log_pro"\t" \
 	 	 	 	 	\
 	 	 	 	 	norm_delta"\t"  to_raw"\t"  norm_to_raw"\t"  variety"\t"  num[i]"\t" \
 	 	 	 	 	\
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
 	 	\
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
 	 	\
 	 	$col_delta \
 	 	\
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
 	 	 	\
 	 	 	name[$col_id]=$col_name
 	 	 	\
 	 	 	formula[$col_id]=$col_formula
 	 	 	\
 	 	 	smiles[$col_id]=$col_smiles
 	 	 	\
 	 	 	inchi[$col_id]=$col_inchi
 	 	 	\
 	 	 	inchikey2D[$col_id]=$col_inchikey2D
 	 	 	\
 	 	 	score[$col_id]=$col_score
 	 	 	\
 	 	 	xlogp[$col_id]=$col_xlogp
 	 	 	\
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
 	 	 	\
 	 	 	"inchi\t"  "inchikey2D\t"  "links\n" >> "com_compound.tsv"
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	printf $i"\t" >> "com_compound.tsv"
 	 	 	 	}
 	 	 	printf simi[$col_list_id]"\t"  name[$col_list_id]"\t"  formula[$col_list_id]"\t"  xlogp[$col_list_id]"\t" \
 	 	 	\
 	 	 	smiles[$col_list_id]"\t" inchi[$col_list_id]"\t"  inchikey2D[$col_list_id]"\t" \
 	 	 	\
 	 	 	links[$col_list_id]"\n" >> "com_compound.tsv"
 	 	 	}
 	 	}
 	}' $data $list
 ##################
 #sort 
 data_o="com_compound.tsv"
 data_r=$(cat <(head -n 1 $data_o) <(sort -t $'\t' -k 6 -n -r <(sed '1d' $data_o)))
 echo "$data_r" > com_compound.tsv
 ##################
 #############################
 #############################
  ## transform principle
 data1="com_compound.tsv"
 data2="results/stat_classification.tsv"
 tolerance="0.01"
 ratio_limit="1.5"
 #similarity_limit="0.5"
 awk -F $'\t' '
 	{
 	if(NR==FNR)
 	 	{
	 	if(NR==1)
	 	 	{
	 	 	for(i=1; i<=NF; i++)
	 	 	 	{
	 	 	 	if($i~/^id/)
	 	 	 	 	{
	 	 	 	 	col_id=i
	 	 	 	 	}
	 	 	 	if($i~/norm_delta/)
	 	 	 	 	{
	 	 	 	 	col_delta=i
	 	 	 	 	}
	 	 	 	if($i~/class/)
	 	 	 	 	{
	 	 	 	 	col_class=i
	 	 	 	 	}
	 	 	 	if($i~/pro\/raw/)
	 	 	 	 	{
	 	 	 	 	col_ratio=i
	 	 	 	 	}
	 	 	 	}
	 	 	}
	 	if(NR>=2)
	 	 	{
	 	 	id[FNR]=$col_id
	 	 	\
	 	 	norm_delta[$col_id]=$col_delta
	 	 	\
	 	 	class[$col_id]=$col_class
	 	 	\
	 	 	ratio[$col_id]=$col_ratio
	 	 	}
	 	}
 	if(NR!=FNR)
 	 	{
 	 	if(FNR==1)
	 	 	{
	 	 	for(i=1; i<=NF; i++)
	 	 	 	{
	 	 	 	if($i~/^id$/)
	 	 	 	 	{
	 	 	 	 	col_id=i
	 	 	 	 	}
	 	 	 	if($i~/definition_source/)
	 	 	 	 	{
	 	 	 	 	col_source=i
	 	 	 	 	}
	 	 	 	if($i~/^definition$/)
	 	 	 	 	{
	 	 	 	 	col_definition=i
	 	 	 	 	}
	 	 	 	if($i~/^level_5$/)
	 	 	 	 	{
	 	 	 	 	col_level=i
	 	 	 	 	}
	 	 	 	if($i~/^subclass$/)
	 	 	 	 	{
	 	 	 	 	col_subclass=i
	 	 	 	 	}
	 	 	 	if($i~/^class$/)
	 	 	 	 	{
	 	 	 	 	col_class=i
	 	 	 	 	}
	 	 	 	if($i~/^superclass$/)
	 	 	 	 	{
	 	 	 	 	col_superclass=i
	 	 	 	 	}
	 	 	 	}
	 	 	}
	 	if(FNR>=2)
	 	 	{
	 	 	if($col_class!="")
	 	 	 	{
	 	 	 	comparison[$col_id]=$col_class
	 	 	 	}
	 	 	else
	 	 	 	{
	 	 	 	comparison[$col_id]=$col_superclass
	 	 	 	}
	 	 	}
 	 	}
 	}
 	END{
 	 	printf "source\t" "target\t" "source_delta\t" "target_delta\t" "trans_probability\t" \
 	 	\
 	 	"ratio\t" "co_classification\n" > "all_principle_'$tolerance'_'$ratio_limit'.tsv"
 	 	\
 	 	for(i in id)
 	 	 	{
 	 	 	if(ratio[id[i]]+0>="'$ratio_limit'"+0)
	 	 	 	{
	 	 	 	for(j in norm_delta)
	 	 	 	 	{
	 	 	 	 	if(norm_delta[id[i]]+0>0 && norm_delta[j]+0<0 || norom_delta[id[i]]+0<0 && norm_delta[j]+0>0)
	 	 	 	 	 	{
	 	 	 	 	 	l=(norm_delta[j])*(-1)-"'$tolerance'"
	 	 	 	 	 	\
	 	 	 	 	 	m=(norm_delta[j])*(-1)+"'$tolerance'"
	 	 	 	 	 	\
	 	 	 	 	 	if(norm_delta[id[i]]+0>=l && norm_delta[id[i]]+0<=m && comparison[id[i]]+0==comparison[j])
	 	 	 	 	 	 	{
	 	 	 	 	 	 	trans=norm_delta[id[i]]+norm_delta[j]
	 	 	 	 	 	 	\
	 	 	 	 	 	 	if(trans+0<0)
	 	 	 	 	 	 	 	{
	 	 	 	 	 	 	 	trans*=(-1)
	 	 	 	 	 	 	 	}
	 	 	 	 	 	 	if(id[i]+0>j)
	 	 	 	 	 	 	 	{
	 	 	 	 	 	 	 	printf id[i]"\t" j"\t" norm_delta[id[i]]"\t" norm_delta[j]"\t" \
	 	 	 	 	 	 	 	\
	 	 	 	 	 	 	 	trans"\t" id[i]":"ratio[id[i]]" | "j":"ratio[j]"\t" \
	 	 	 	 	 	 	 	\
	 	 	 	 	 	 	 	comparison[id[i]]"\n" >> "all_principle_'$tolerance'_'$ratio_limit'.tsv"
	 	 	 	 	 	 	 	}
	 	 	 	 	 	 	if(id[i]+0<j)
	 	 	 	 	 	 	 	{
	 	 	 	 	 	 	 	printf j"\t" id[i]"\t" norm_delta[j]"\t" norm_delta[id[i]]"\t" \
	 	 	 	 	 	 	 	\
	 	 	 	 	 	 	 	trans"\t" j":"ratio[j]" | "id[i]":"ratio[id[i]]"\t" \
	 	 	 	 	 	 	 	\
	 	 	 	 	 	 	 	comparison[id[i]]"\n" >> "all_principle_'$tolerance'_'$ratio_limit'.tsv"
	 	 	 	 	 	 	 	}
	 	 	 	 	 	 	}
	 	 	 	 	 	}
	 	 	 	 	}
	 	 	 	}
 	 	 	delete norm_delta[id[i]]
 	 	 	}
 	 	}' $data1 $data2
 ######################
 ######################
 #data(alignment)
 data1=neg.csv
 data2=pos.csv
 awk -F "," '
  	{
  	if(NR==FNR)
  	 	{
  	 	if(FNR==1)
  	 	 	{
  	 	 	for(i=1; i<=NF; i++)
  	 	 	 	{
  	 	 	 	if($i~/row ID/)
  	 	 	 	 	{
  	 	 	 	 	col_id=i
  	 	 	 	 	}
  	 	 	 	else if($i~/row m\/z/)
  	 	 	 	 	{
  	 	 	 	 	col_mz=i
  	 	 	 	 	}
  	 	 	 	else if($i~/row retention time/)
  	 	 	 	 	{
  	 	 	 	 	col_rt=i
  	 	 	 	 	}
  	 	 	 	else
  	 	 	 	 	{
  	 	 	 	 	col_name[$i]=i
  	 	 	 	 	}
  	 	 	 	}
  	 	 	n=asorti(col_name,sort_name)
  	 	 	}
  	 	if(FNR>=2)
  	 	 	{
  	 	 	for(i=1; i<=n; i++)
  	 	 	 	{
  	 	 	 	}
  	 	 	}
  	 	}
  	}' $data1 $data2
 #######################
 #######################
 #######################
 #search glycosides
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
 ##database idenfication (mgf vs msp)
 data1="com_lignans_and_iridoids.tsv"
 datasource="../initial_8_neg.mgf"
 data2="../datasource.mgf"
 database="../MSMS-Public-Neg-VS15.msp"
 data3="../database.msp"
 tolerance=0.01
 ms2_tolerance=0.2
 in_tolerance=100
 weight1=80
 weight2=20
 point_limit=80
 sed 's/\r//g' $database > ../database.msp
 sed 's/\r//g' $datasource > ../datasource.mgf
 awk -F $'\t' '
 	{
 	if(NR==FNR)
 	 	{
 	 	if(NR==1)
 	 	 	{
 	 	 	p_file=FILENAME
 	 	 	\
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if($i~/^id/)
 	 	 	 	 	{
 	 	 	 	 	col_id=i
 	 	 	 	 	}
 	 	 	 	if($i~/m\/z/)
 	 	 	 	 	{
 	 	 	 	 	col_mz=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(NR>=2)
 	 	 	{
 	 	 	pre_mz[$col_id]=$col_mz
 	 	 	\
 	 	 	printf "info: datafile1 is ready. The data is " $col_id " to " $col_mz "\n"  #info
 	 	 	}
 	 	}
 	if(FILENAME~/datasource.mgf/)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	close(p_file)
 	 	 	\
 	 	 	p_file=FILENAME
 	 	 	\
 	 	 	FS="[=]||[ ]"
 	 	 	\
 	 	 	printf "info: datafile2 is ready." #info
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	if($0~/FEATURE_ID/)
 	 	 	 	{
 	 	 	 	id=$2
 	 	 	 	}
 	 	 	if($0~/MSLEVEL/)
 	 	 	 	{
 	 	 	 	ms_level=$2
 	 	 	 	}
 	 	 	if(($0~/^[0-9]/) && ms_level=2)
 	 	 	 	{
 	 	 	 	count[id]+=1
 	 	 	 	\
 	 	 	 	data2_mz[id,count[id]]=$1
 	 	 	 	\
 	 	 	 	data2_in[id,count[id]]=$2
 	 	 	 	\
 	 	 	 	if(data2_max_intensity[id]<$2)
 	 	 	 	 	{
 	 	 	 	 	data2_max_intensity[id]=$2
 	 	 	 	 	}
 	 	 	 	printf "Info: "$1"_"$2" >>> "id"\n"
 	 	 	 	}
 	 	 	}
 	 	}
 	if(FILENAME~/.msp/)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	for(i in data2_in)
 	 	 	 	{
 	 	 	 	split(i,p,"[\034]")
 	 	 	 	\
 	 	 	 	data2_in[i]*=(100/data2_max_intensity[p[1]])
 	 	 	 	}
 	 	 	msp_num+=1
 	 	 	\
 	 	 	FS="[:][ ]||[\t]"
 	 	 	\
 	 	 	close(p_file)
 	 	 	\
 	 	 	p_file=FILENAME
 	 	 	\
 	 	 	n=split($0,c,"[:][ ]")
 	 	 	\
 	 	 	name=c[n]
 	 	 	\
 	 	 	printf "Info: datafile3 is ready. The first name is " name "\n"  #info
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	if($1~/NAME/)
 	 	 	 	{
 	 	 	 	if(belong_id[msp_num]+0>=1)
 	 	 	 	 	{
 	 	 	 	 	for(i=1; i<=database_count[msp_num]; i++)
 	 	 	 	 	 	{
 	 	 	 	 	 	mz=data3_mz[msp_num,i]
 	 	 	 	 	 	\
 	 	 	 	 	 	intensity=(data3_in[msp_num,i]/max_in[msp_num])*100
 	 	 	 	 	 	\
 	 	 	 	 	 	if(intensity+0>=5)
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	point='$weight1'
 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	full_marks[msp_num]+=point
 	 	 	 	 	 	 	}
 	 	 	 	 	 	else
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	point='$weight2'
 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	full_marks[msp_num]+=point
 	 	 	 	 	 	 	}
 	 	 	 	 	 	for(j=1; j<=belong_id[msp_num]; j++)
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	for(k=1; k<=count[subdirectory[msp_num,j]]; k++)
 	 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	 	l1=data2_mz[subdirectory[msp_num,j],k]-'$ms2_tolerance'
 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	m1=data2_mz[subdirectory[msp_num,j],k]+'$ms2_tolerance'
 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	if(mz+0>=l1 && mz+0<=m1+0)
 	 	 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	 	 	l2=data2_in[subdirectory[msp_num,j],k]-'$in_tolerance'
 	 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	 	m2=data2_in[subdirectory[msp_num,j],k]+'$in_tolerance'
 	 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	 	if(intensity+0>=l2+0 && intensity+0<=m2+0)
 	 	 	 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	 	 	 	assign_point[subdirectory[msp_num,j],msp_num]+=point;
 	 	 	 	 	 	 	 	 	 	\
 	 	 	 	 	 	 	 	 	 	break;
 	 	 	 	 	 	 	 	 	 	}
 	 	 	 	 	 	 	 	 	}
 	 	 	 	 	 	 	 	}
 	 	 	 	 	 	 	}
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	msp_num+=1
 	 	 	 	\
 	 	 	 	name=$2
 	 	 	 	}
 	 	 	if($1~/PRECURSORMZ/)
 	 	 	 	{
 	 	 	 	for(i in pre_mz)
 	 	 	 	 	{
 	 	 	 	 	if($2+0>=pre_mz[i]-'$tolerance' && $2<=pre_mz[i]+'$tolerance')
 	 	 	 	 	 	{
 	 	 	 	 	 	printf "Info: the MS1 are " $2 " vs " pre_mz[i] "\n"
 	 	 	 	 	 	\
 	 	 	 	 	 	data3_mz[msp_num]=$2
 	 	 	 	 	 	\
 	 	 	 	 	 	belong_id[msp_num]+=1
 	 	 	 	 	 	\
 	 	 	 	 	 	subdirectory[msp_num,belong_id[msp_num]]=i
 	 	 	 	 	 	\
 	 	 	 	 	 	assign[i,msp_num]=name
 	 	 	 	 	 	\
 	 	 	 	 	 	sep[i,msp_num]=msp_num
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	if(belong_id[msp_num]+0>=1)
 	 	 	 	{
 	 	 	 	if($1~/PRECURSORTYPE/)
 	 	 	 	 	{
 	 	 	 	 	data3_type[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_formula[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_ontology[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_inchikey[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_smiles[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_rt[msp_num]=$2; getline;
 	 	 	 	 	\
 	 	 	 	 	data3_ccs[msp_num]=$2;
 	 	 	 	 	}
 	 	 	 	if($1~/Num Peaks/)
 	 	 	 	 	{
 	 	 	 	 	data3_peaks[msp_num]=$2
 	 	 	 	 	}
 	 	 	 	if($1~/^[0-9]/)
 	 	 	 	 	{
 	 	 	 	 	database_count[msp_num]+=1
 	 	 	 	 	\
 	 	 	 	 	data3_mz[msp_num,database_count[msp_num]]=$1
 	 	 	 	 	\
 	 	 	 	 	data3_in[msp_num,database_count[msp_num]]=$2
 	 	 	 	 	\
 	 	 	 	 	if(max_in[msp_num]<$2)
 	 	 	 	 	 	{
 	 	 	 	 	 	max_in[msp_num]=$2
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}
 	if(msp_num!="" && (FILENAME~"'$data1'"))
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	FS="\t"
 	 	 	\
 	 	 	close(p_file)
 	 	 	\
 	 	 	printf $0"\t" "custom_idenfication\n" > "../idenfication_'$data1'.tsv"
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	printf $0 >> "../idenfication_'$data1'.tsv"
 	 	 	\
 	 	 	for(i in assign_point)
 	 	 	 	{
 	 	 	 	if(i~"^"$1"\034" && assign_point[i]+0>='$point_limit')
 	 	 	 	 	{
 	 	 	 	 	score[i]=assign_point[i]
 	 	 	 	 	}
 	 	 	 	}
 	 	 	for(n=asort(score,sort_score); n>=1; n--)
 	 	 	 	{
 	 	 	 	for(i in score)
 	 	 	 	 	{
 	 	 	 	 	if(score[i]==sort_score[n])
 	 	 	 	 	 	{
 	 	 	 	 	 	printf "\t" score[i] "("  full_marks[sep[i]]  ")"  " | " \
 	 	 	 	 	 	\
 	 	 	 	 	 	data3_type[sep[i]] " | "  \
 	 	 	 	 	 	\
 	 	 	 	 	 	data3_formula[sep[i]] " | " \
 	 	 	 	 	 	\
 	 	 	 	 	 	data3_ontology[sep[i]] " | " \
 	 	 	 	 	 	\
 	 	 	 	 	 	assign[i] " | " \
 	 	 	 	 	 	\
 	 	 	 	 	 	>> "../idenfication_'$data1'.tsv"
 	 	 	 	 	 	\
 	 	 	 	 	 	delete score[i]; break;
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	delete score
 	 	 	\
 	 	 	printf "\n" >> "../idenfication_'$data1'.tsv"
 	 	 	}
 	 	}
 	}' $data1 $data2 $data3 $data1
 #############################
 #test xlogp parameter
 data="com_lignans_and_iridoids.tsv"
 awk -F $'\t' '
 	{
 	if(NR==1)
 	 	{
	 	for(i=1; i<=NF; i++)
	 	 	{
	 	 	if($i~/similarity/)
	 	 	 	{
	 	 	 	col_similarity=i
	 	 	 	}
	 	 	}
	 	printf $0"\n" > "0.4_'$data'"
	 	}
	if(NR>=2)
	 	{
	 	if($col_similarity+0>=0.4)
	 	 	{
	 	 	printf $0"\n" >> "0.4_'$data'"
	 	 	}
	 	}
 	}' $data
 data="0.5_com_compound.tsv"
 awk -F $'\t' '
  	{
  	if(NR==1)
  	 	{
  	 	for(i=1; i<=NF; i++)
  	 	 	{
  	 	 	if($i~/similarity/)
	 	 	 	{
	 	 	 	col_similarity=i
	 	 	 	}
	 	 	if($i~/^rt$/)
	 	 	 	{
	 	 	 	col_rt=i
	 	 	 	}
	 	 	if($i~/xlogp/)
	 	 	 	{
	 	 	 	col_xlogp=i
	 	 	 	}
	 	 	if($i~/m\/z/)
	 	 	 	{
	 	 	 	col_mz=i
	 	 	 	}
  	 	 	}
  	 	}
  	if(NR>=2)
  	 	{
  	 	if($col_xlogp!="" && $col_mz+0>=400 && $col_mz+0<=800)
  	 	 	{
  	 	 	n+=1
  	 	 	\
	  	 	sum_x+=$col_rt
	  	 	\
	  	 	figure_x[n]=$col_rt
	  	 	\
	  	 	sum_y+=$col_xlogp
	  	 	\
	  	 	figure_y[n]=$col_xlogp
	  	 	}
  	 	}
  	}
  	END{
  	 	average_x=sum_x/n
  	 	\
  	 	average_y=sum_y/n
  	 	\
  	 	for(i=1; i<=n; i++)
  	 	 	{
  	 	 	b_up+=((figure_x[i]-average_x)*(figure_y[i]+(-1)*average_y))
  	 	 	\
  	 	 	b_down+=((figure_x[i]-average_x)^2)
  	 	 	\
  	 	 	r_up_square_root=b_up
  	 	 	\
  	 	 	r_down_1+=((figure_x[i]-average_x)^2)
  	 	 	\
  	 	 	r_down_2+=((figure_y[i]+(-1)*average_y)^2)
  	 	 	}
  	 	b=b_up/b_down
  	 	\
  	 	a=average_y-b*average_x
  	 	\
  	 	r_square=(r_up_square_root^2)/(r_down_1*r_down_2)
  	 	\
  	 	printf "b="b  "\t"  "a="a  "\t"  "y=" b "x+" a "\t"  "r²="r_square  "\n"  average_x "\t"  average_y "\n"
  	 	}' $data
 ###########################
 ###B_ratio
 ################# need revise
 data1="../mobile_phase_system.tsv"
 data2="0.5_com_compound.tsv"
 awk -F $'\t' '
  	{
  	if(NR==FNR)
  	 	{
  	 	if(FNR==1)
  	 	 	{
	  	 	for(i=1; i<=NF; i++)
	  	 	 	{
	  	 	 	if($i~/^time/)
	  	 	 	 	{
	  	 	 	 	col_time=i
	  	 	 	 	}
	  	 	 	if($i~/^B/)
	  	 	 	 	{
	  	 	 	 	col_B=i
	  	 	 	 	}
	  	 	 	}
	  	 	}
  	 	if(FNR>=2)
  	 	 	{
  	 	 	x_rt[FNR-1]=$col_time
  	 	 	max_rows=FNR-1
  	 	 	x1=$col_time
  	 	 	y1=$col_B
  	 	 	getline; 
  	 	 	x2=$col_time
  	 	 	y2=$col_B
  	 	 	if(x2!="")
  	 	 	 	{
  	 	 	 	b[FNR-1]=(y2-y1)/(x2-x1)
  	 	 	 	a[FNR-1]=y1-b[FNR-1]*x1
  	 	 	 	}
  	 	 	}
  	 	}
  	if(FILENAME~/'$data2'/)
  	 	{
  	 	if(FNR==1)
  	 	 	{
  	 	 	printf "" > "b_'$data2'"
  	 	 	for(i=1; i<=NF; i++)
	  	 	 	{
	 	 	 	if($i~/^rt$/)
	 	 	 	 	{
	 	 	 	 	col_rt=i
	 	 	 	 	}
	 	 	 	printf $i >> "b_'$data2'"
	 	 	 	if(i!=NF)
	 	 	 	 	{
	 	 	 	 	printf "\t" >> "b_'$data2'"
	 	 	 	 	}
	 	 	 	else
	 	 	 	 	{
	 	 	 	 	printf "\n" >> "b_'$data2'"
	 	 	 	 	}
	 	 	 	if(col_rt==i)
	 	 	 	 	{
	 	 	 	 	printf "B_ratio\t" >> "b_'$data2'"
	 	 	 	 	}
	  	 	 	}
  	 	 	}
  	 	if(FNR>=2)
  	 	 	{
  	 	 	for(i=1; i<=NF; i++)
  	 	 	 	{
  	 	 	 	printf $i >> "b_'$data2'"
	 	 	 	if(i!=NF)
	 	 	 	 	{
	 	 	 	 	printf "\t" >> "b_'$data2'"
	 	 	 	 	}
	 	 	 	else
	 	 	 	 	{
	 	 	 	 	printf "\n" >> "b_'$data2'"
	 	 	 	 	}
	 	 	 	if(col_rt==i)
	 	 	 	 	{
	 	 	 	 	for(j=1; j<=max_rows; j++)
	 	 	 	 	 	{
	 	 	 	 	 	if(x_rt[j]+0<=$col_rt+0 && x_rt[j+1]+0>=$col_rt+0)
	 	 	 	 	 	 	{
	 	 	 	 	 	 	B_ratio=b[j]*$col_rt+a[j]
	 	 	 	 	 	 	printf B_ratio >> "b_'$data2'"
	 	 	 	 	 	 	break;
	 	 	 	 	 	 	}
	 	 	 	 	 	}
	 	 	 	 	printf "\t" >> "b_'$data2'"
	 	 	 	 	}
  	 	 	 	}
  	 	 	}
  	 	}
  	}' $data1 $data2
 ###########################
 ### filter_class
 data="results/stat_classification.tsv"
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
 	 	 	printf "" > "filter_0_class.tsv"
 	 	 	printf class[i]"\n" >> "filter_0_class.tsv"
 	 	 	}
 	 	}' $data
 cat filter_0_class.tsv | trans :zh -b > filter_0_trans_class.tsv
 ###########################
 data1="filter_0_class.tsv"
 data2="results/canopus_pp.tsv"
 savepath="results/canopus_pp_filter.tsv"
 awk -F $'\t' '
  	{
  	if(NR==FNR)
  	 	{
	  	if(FNR==1)
	  	 	{
	  	 	}
	  	if(FNR>=2)
	  	 	{
	  	 	}
	  	}
	if(NR>FNR)
	 	{
	 	}
  	}' $data
 ###########################
 data1="filter_class.tsv"
 data2="com_compound.tsv"
 savename="less_com_compound.tsv"
 awk -F $'\t' '
  	{
  	if(NR==FNR)
  	 	{
  	 	filter_class[$2]=$1
  	 	}
  	if(NR>FNR)
  	 	{
  	 	if(FNR==1)
  	 	 	{
  	 	 	for(i=1; i<=NF; i++)
  	 	 	 	{
  	 	 	 	if($i~/classification/)
  	 	 	 	 	{
  	 	 	 	 	col_class=i
  	 	 	 	 	}
  	 	 	 	if($i~/pro\/raw/)
  	 	 	 	 	{
  	 	 	 	 	col_ratio=i
  	 	 	 	 	}
  	 	 	 	}
  	 	 	printf $0"\n" > "'$savename'"
  	 	 	print col_class,col_ratio
  	 	 	}
  	 	if(FNR>=2)
  	 	 	{
  	 	 	for(i in filter_class)
  	 	 	 	{
  	 	 	 	if($col_class==filter_class[i] && $col_ratio+0>=0)
  	 	 	 	 	{
  	 	 	 	 	printf $0"\n" >> "'$savename'"
  	 	 	 	 	\
  	 	 	 	 	break;
  	 	 	 	 	}
  	 	 	 	}
  	 	 	}
  	 	}
  	}' $data1 $data2
 ###################
 #boxplot
 Rscript=$(echo '
library(ggplot2)
source<-read.table(file="boxplot.tsv",header=T,sep="\t")
boxplot<-ggplot(source,aes(x=classification,y=pro.raw,fill=classification))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot()+
  geom_point(position="jitter",shape=21,size=3)+
  theme(
       legend.position="none",
       axis.text.x = element_blank(),
       text = element_text(size=30),
       )+
  geom_text(data=source, aes(x=classification, y=30, label=classification),
	    color="black", fontface="bold", size=8, angle= 90, 
	    inherit.aes = FALSE )+
  stat_summary(fun="mean",geom="point",shape=23,size=2.5,fill="grey")+
  facet_zoom(ylim = c(0, 10))
pdf("boxplot.pdf",width=25,height=15)
boxplot
dev.off()
')
 ##################
 ###########################
 ###for sunplot
 #data1="filter_class.tsv"
 data1="filter_class.csv" #lignans and iridoids
 data2="results/canopus_pp.tsv"
 data3="com_compound.tsv"
 savename="for_sun.tsv"
 #ex_export="com_lignans_and_iridoids.tsv"
 ex_export="com_carboxylic_acids.tsv"
 similarity_limit="0.4"
 class_pp_limit="0.9"
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
  	 	 	\
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
	  	 	 	 	 	printf $col_id"\t"  class_set[i]"\t"  $col_log_raw"\t"  $col_log_pro"\n" \
	  	 	 	 	 	\
	  	 	 	 	 	>> "'$savename'"
	  	 	 	 	 	\
	  	 	 	 	 	printf $0"\n" >> "'$ex_export'"
	  	 	 	 	 	}
	  	 	 	 	}
	  	 	 	}
  	 	 	}
  	 	}
  	}' $data1 $data2 $data3
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
 	 	 	 	 	\
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
  	 	 	 	printf sprintf("%.2f",rt_min)"\t"  $col_intensity"\t"  $col_sample"\t"  group_name[$col_sample]"\t" \
  	 	 	 	\
  	 	 	 	label_sig[id,$col_sample]"\t"  color[id,$col_sample] >> "'$savepath'" id ".tsv"
  	 	 	 	\
  	 	 	 	start_FNR[id]+=1
  	 	 	 	\
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
 ##use python rdkit 
 conda activate my-rdkit-env
 savepath="results/structure_2d"
 cd $savepath
 data="../../com_compound.tsv"
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
  	 	printf $col_id "_" $col_smiles "|||"
  	 	}
  	}' $data > list
 python_script="/home/wizard/Downloads/codes/python_files/draw_structure.py"
 python $python_script list
 #### build instance data
 data1="0703_all/lignans_and_iridoids.tsv"
 data2="0703_all/ftalign.tsv"
 savepath="instance_data.tsv"
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
  	 	 	 	if($i~/^classification/)
  	 	 	 	 	{
  	 	 	 	 	col_class=i
  	 	 	 	 	}
  	 	 	 	}
  	 	 	}
  	 	if(FNR>=2)
  	 	 	{
  	 	 	if($col_class~/Lignans/)
  	 	 	 	{
  	 	 	 	the_id[$col_id]=$col_id
  	 	 	 	}
  	 	 	}
  	 	}
  	if(NR>FNR)
  	 	{
  	 	if(FNR==1)
  	 	 	{
  	 	 	printf $1 > "'$savepath'"
  	 	 	for(i=2; i<=NF; i++)
  	 	 	 	{
  	 	 	 	n=split($i, a, "[_]")
  	 	 	 	for(j in the_id)
  	 	 	 	 	{
  	 	 	 	 	if(a[n]~"^"j"$")
  	 	 	 	 	 	{
  	 	 	 	 	 	col_the_id[a[n]]=i
  	 	 	 	 	 	printf "\t"a[n] >> "'$savepath'"
  	 	 	 	 	 	break;
  	 	 	 	 	 	}
  	 	 	 	 	}
  	 	 	 	}
  	 	 	printf "\n" >> "'$savepath'"
  	 	 	}
  	 	if(FNR>=2)
  	 	 	{
  	 	 	n=split($1, b, "[_]")
  	 	 	printf "Get row of "  FNR  " >>> " b[n] "\n"
  	 	 	for(k in the_id)
  	 	 	 	{
  	 	 	 	if(b[n]==k)
  	 	 	 	 	{
	 	  	 	 	printf b[n] >> "'$savepath'"
	 	  	 	 	for(i=2; i<=NF; i++)
	 	  	 	 	 	{
	 	  	 	 	 	for(j in col_the_id)
	 	  	 	 	 	 	{
	 	  	 	 	 	 	if(i==col_the_id[j])
	 	  	 	 	 	 	 	{
	 	  	 	 	 	 	 	printf "\t"$i >> "'$savepath'"
	 	  	 	 	 	 	 	break
	 	  	 	 	 	 	 	}
	 	  	 	 	 	 	}
	 	  	 	 	 	}
	 	  	 	 	printf "\n" >> "'$savepath'"
	 	  	 	 	break;
	  	 	 	 	}
	  	 	 	}
  	 	 	}
  	 	}
  	}' $data1 $data2 
 #### norm the instance data
 data="instance_data.tsv"
 savepath="norm_instance_data.tsv"
 awk -F $'\t' '
  	{
 	for(i=1; i<=NF; i++)
 	 	{
 	 	raw[NR,i]=$i
 	 	}
 	}
	END{ 	
	 	printf raw[1,1] > "'$savepath'"
	 	for(i=2; i<=NF; i++)
	 	 	{
	 	 	printf "\t"raw[1,i] >> "'$savepath'"
	 	 	}
	 	printf"\n" >> "'$savepath'"
 	 	for(a=2; a<=NR; a++)
 	 	 	{
 	 	 	printf raw[a,1] >> "'$savepath'"
	 	 	for(b=2; b<=NF; b++)
	 	 	 	{
	 	 	 	norm1[a,b]=raw[a,b]/raw[a,a];
	 	 	 	norm2[a,b]=raw[a,b]/raw[b,b];
	 	 	 	norms[a,b]=(norm1[a,b]+norm2[a,b])/2;
	 	 	 	printf "\t"norms[a,b] >> "'$savepath'"
	 	  	 	}
	 	  	printf "\n" >> "'$savepath'"
	  	 	}
 	 	}' $data
 #### cut the data
 data1="instance_data.tsv"
 data2="norm_instance_data.tsv"
 cut="51"
 awk -F $'\t' '
  	{
  	if(NR==FNR)
  	 	{
  	 	if(FNR<='$cut')
  	 	 	{
  	 	 	printf $1 > "cut_'$data1'"
  	 	 	for(i=2; i<='$cut'; i++)
  	 	 	 	{
  	 	 	 	printf "\t"$i >> "cut_'$data1'"
  	 	 	 	}
  	 	 	printf "\n" >> "cut_'$data1'"
  	 	 	}
  	 	}
  	if(NR>FNR)
  	 	{
  	 	if(FNR<='$cut')
  	 	 	{
  	 	 	printf $1 > "cut_'$data2'"
  	 	 	for(i=2; i<='$cut'; i++)
  	 	 	 	{
  	 	 	 	printf "\t"$i >> "cut_'$data2'"
  	 	 	 	}
  	 	 	printf "\n" >> "cut_'$data2'"
  	 	 	}
  	 	}
  	}' $data1 $data2
 ####
 Rscript ~/Downloads/codes/ggplot2_heatmap.R
 ################################### Add to Methodology Content ################################################################################
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
  	 	 	\
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
	  	 	 	 	 	printf $col_id"\t"  class_set[i]"\t"  $col_log_raw"\t"  $col_log_pro"\n" \
	  	 	 	 	 	\
	  	 	 	 	 	>> "'$savename'"
	  	 	 	 	 	\
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
  	 	 	printf stat_id[i]"\t" stat_id[i]"\t"  "1\t"  "0\t"  "null\t"  "null\t"  belong[i]"\n" \
  	 	 	\
  	 	 	>> "'$savepath'" belong[i] ".tsv"
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
 data1="filter_0_class.tsv" #lignans and iridoids
 data2="results/canopus_pp.tsv"
 data3="com_compound.tsv"
 ex_export="com_1011.tsv"
 similarity_limit="0.4"
 ###################################
 class_pp_limit="0.9"
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
  	 	 	\
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
	  	 	 	 	 	printf $col_id"\t"  class_set[i]"\t"  $col_log_raw"\t"  $col_log_pro"\n" \
	  	 	 	 	 	\
	  	 	 	 	 	>> "'$savename'"
	  	 	 	 	 	\
	  	 	 	 	 	printf $0"\n" >> "'$ex_export'"
	  	 	 	 	 	}
	  	 	 	 	}
	  	 	 	}
  	 	 	}
  	 	}
  	}' $data1 $data2 $data3
 ###################################
 #### stat num
 class_pp_limit=0.9
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
 #################
 #################
