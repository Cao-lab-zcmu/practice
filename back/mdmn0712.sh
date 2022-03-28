##########################
##########################
echo "We are all in the gutter, 
but some of us are looking at the stars.";

PS3='Please select the workflow to be executed. >>> '

select command in "default" "structure_extract" "classification_extract_sum" "classification_extract_filter" "fragment_tree_network" "fragment_tree_network_delta" "double_ion_network" "exit"
do
	default="structure_extract classification_extract_sum classification_extract_filter fragment_tree_network fragment_tree_network_delta"
	
	if [[ $command == "default" ]]
	then
	confirm=0
	
	 	until [[ $confirm == "yes" ]] || [[ $confirm == "no" ]]
	 	do
	 	read -p "Running all proccesses sequentially? [yes/no] >>> " confirm
	 	done;
	 	if [[ $confirm == "no" ]]
	 	then exit
	 	fi;
	 	
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

fingerid=$(awk -F $'\t' '
 	BEGIN{
 	 	all_id=0
 	 	max_id=-1
 	 	p_id="null"
 	 	}
 	{
 	if(FNR==1)
 	 	{
 	 	if(NR>FNR)
 	 	 	{
 	 	 	close(pfile)
 	 	 	}
 	 	split(FILENAME,a,"[/]");
 	
 	 	n=split(a[1],b,"[_]");
 	
 	 	id=b[n];
 	
 	 	if(id!=p_id)
 	 	 	{
 	 	 	all_id+=1
 	 	 	\
 	 	 	if(id>max_id)
 	 	 	 	{
 	 	 	 	max_id=id
 	 	 	 	}
 	 	 	}
 	 	file[id]=a[1]
 	 	
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
 	 	 	 	}
 	 	 	}
 	 	pfile=FILENAME
 	 	}
 	if(FNR>=2)
 	 	{
 	 	row_score=$col_score
 	 	
 	 	row_simi=$col_simi
 	 	
 	 	if(max["score",id]=="" || max["score",id]<row_score)
	 	 	{
 	 	 	max["score",id]=row_score
 	 	 	
 	 	 	name["score",id]=$col_name
 	 	
 	 	 	formula["score",id]=$col_formu
 	 	 	
 	 	 	simi["score",id]=$col_simi
 	 	
 	 	 	smiles["score",id]=$col_smiles
 	 	
 	 	 	inchi["score",id]=$col_inchi
 	 	
 	 	 	in2D["score",id]=$col_2D
 	 	 	
 	 	 	score["score",id]=$col_score
 	 	
 	 	 	xlogp["score",id]=$col_x
 	 	 	}
 	 	if(max["simi",id]=="" || max["simi",id]<row_simi)
	 	 	{
 	 	 	max["simi",id]=row_simi
 	 	 	
 	 	 	name["simi",id]=$col_name
 	 	
 	 	 	formula["simi",id]=$col_formu
 	 	 	
 	 	 	simi["simi",id]=$col_simi
 	 	
 	 	 	smiles["simi",id]=$col_smiles
 	 	
 	 	 	inchi["simi",id]=$col_inchi

 	 	 	in2D["simi",id]=$col_2D
 	 	 	
 	 	 	score["simi",id]=$col_simi
 	 	 
 	 	 	xlogp["simi",id]=$col_x
 	 	 	}
 	 	}
 	}
 	END{
 	 	printf "id\t"  "name\t"  "formula\t"  "similarity\t"  "smiles\t"  "inchi\t"  "inchikey2D\t" \
 	 	\
 	 	"score\t"  "xlogp\n";
 	 	\
 	 	for(l=1; l<=max_id; l++)
 	 	 	{
 	 	 	if(max["score",l]!="")
 	 	 	 	{
 	 	 	 	printf l"\t" \
 	 	 	 	\
 	 	 	 	name["score",l]"\t"  formula["score",l]"\t"  simi["score",l]"\t" \
 	 	 	 	\
 	 	 	 	smiles["score",l]"\t"  inchi["score",l]"\t"  in2D["score",l]"\t" \
 	 	 	 	\
 	 	 	 	score["score",l]"\t"  xlogp["score",l]"\n"
 	 	 	 	
 	 	 	 	printf formula["score",l]"\t"  file[l]"\n" >> "temp/Mo_filename"
 	 	 	 	}
 	 	 	if(max["simi",l]!="" && max["simi",l]!=max["score",l])
 	 	 	 	{
 	 	 	 	printf l"\t" \
 	 	 	 	\
 	 	 	 	name["simi",l]"\t"  formula["simi",l]"\t"  simi["simi",l]"\t" \
 	 	 	 	\
 	 	 	 	smiles["simi",l]"\t"  inchi["simi",l]"\t"  in2D["simi",l]"\t" \
 	 	 	 	\
 	 	 	 	score["simi",l]"\t"  xlogp["simi",l]"\n"
 	 	 	 	}
 	 	 	}
 	 	}' *_*/fingerid/*.tsv)

echo "$fingerid" > results/fingerid_sum.tsv

echo "structure_extract results have been successfully assembled into <results/fingerid_sum.tsv>";
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
	
datas="*_*/canopus/*.fpt"
	
canopus=$(awk -F $'\t' '
 	{
 	if(NR==FNR)
 	 	{
 	 	file[FNR]=$2
 	 	
 	 	mo[$2]=$1
 	 	
 	 	i=split($2,s,"[_]")
 	 	
 	 	the_id[FNR]=s[i]
 	 	
 	 	n=FNR
 	 	}
 	if((FILENAME ~ /'$data2'/))
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	close("'$data1'")
 	 	 	
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if(($i ~ /name/))
 	 	 	 	 	{
 	 	 	 	 	col_class=i
 	 	 	 	 	}
 	 	 	 	if(($i ~ /absolute/))
 	 	 	 	 	{
 	 	 	 	 	col_abindex=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	abindex[1,FNR]=$col_abindex
 	 	 	
 	 	 	class[1,$col_abindex]=$col_class
 	 	 	
 	 	 	rows_data2=FNR
 	 	 	}
 	 	}
 	if((FILENAME ~ /'$data3'/))
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	close("'$data2'")
 	 	 	
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if(($i ~ /name/))
 	 	 	 	 	{
 	 	 	 	 	col_class=i
 	 	 	 	 	}
 	 	 	 	if(($i ~ /absolute/))
 	 	 	 	 	{
 	 	 	 	 	col_abindex=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	abindex[2,FNR]=$col_abindex
 	 	 	
 	 	 	class[2,$col_abindex]=$col_class
 	 	 	
 	 	 	rows_data3=FNR
 	 	 	}
 	 	}
 	if((FILENAME ~ /fpt/))
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	if(p_filename=="")
 	 	 	 	{
 	 	 	 	close("'$data3'")
 	 	 	 	}
 	 	 	else if(FILENAME!=p_filename)
 	 	 	 	{
 	 	 	 	close(p_filename)
 	 	 	 	}
 	 	 	p_filename=FILENAME
 	 	 	\
 	 	 	for(i=1; i<=n; i++)
 	 	 	 	{
 	 	 	 	if((FILENAME ~ file[i]) && (FILENAME ~ mo[file[i]]))
 	 	 	 	 	{
 	 	 	 	 	break
 	 	 	 	 	}
 	 	 	 	else
 	 	 	 	 	{
 	 	 	 	 	x=i+1
 	 	 	 	 	}
 	 	 	 	}
 	 	 	if(x==n+1)
 	 	 	 	{
 	 	 	 	nextfile
 	 	 	 	}
 	 	 	split(FILENAME,a,"[/]");

 	 	 	m=split(a[1],b,"[_]");

 	 	 	id=b[m];

 	 	 	if((FILENAME ~ /\+.fpt/))
 	 	 	 	{
 	 	 	 	ion=1
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	ion=2
 	 	 	 	}
 	 	 	pp[id,abindex[ion,2]]=$1
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	pp[id,abindex[ion,FNR+1]]=$1
 	 	 	}
 	 	}
 	}
 	END{
 	 	if(abindex[1,rows_data2] > abindex[2,rows_data3])
 	 	 	{
 	 	 	maxindex=abindex[1,rows_data2]
 	 	 	}
 	 	else
 	 	 	{
 	 	 	maxindex=abindex[2,rows_data3]
 	 	 	}
 	 	printf "id\t"
 	 	\
 	 	for(i=0; i<=maxindex; i++)
 	 	 	{
 	 	 	if(class[1,i]!="" && i!=maxindex)
 	 	 	 	{
 	 	 	 	printf class[1,i]"\t"
 	 	 	 	}
 	 	 	else if(class[2,i]!="" && i!=maxindex)
 	 	 	 	{
 	 	 	 	printf class[2,i]"\t"
 	 	 	 	}
 	 	 	else if(i==maxindex && class[1,i]!="")
 	 	 	 	{
 	 	 	 	printf class[1,i]"\n"
 	 	 	 	}
 	 	 	else if(i==maxindex && class[2,i]!="")
 	 	 	 	{
 	 	 	 	printf class[2,i]"\n"
 	 	 	 	}
 	 	 	}
 	 	for(i=1; i<=n; i++)
 	 	 	{
 	 	 	printf the_id[i]"\t"
 	 	 	
 	 	 	for(j=0; j<=maxindex; j++)
 	 	 	 	{
 	 	 	 	if(class[1,j]!="" && j!=maxindex || class[2,j]!="" && j!=maxindex)
 	 	 	 	 	{
 	 	 	 	 	if(pp[the_id[i],j]=="")
 	 	 	 	 	 	{
 	 	 	 	 	 	printf 0"\t"
 	 	 	 	 	 	}
 	 	 	 	 	else
 	 	 	 	 	 	{
 	 	 	 	 	 	printf pp[the_id[i],j]"\t"
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	else if(class[1,j]!="" && j==maxindex || class[2,j]!="" && j==maxindex)
 	 	 	 	 	{
 	 	 	 	 	if(pp[the_id[i],j]=="")
 	 	 	 	 	 	{
 	 	 	 	 	 	printf 0"\n"
 	 	 	 	 	 	}
 	 	 	 	 	else
 	 	 	 	 	 	{
 	 	 	 	 	 	printf pp[the_id[i],j]"\n"
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}' $data1 $data2 $data3 $datas)

echo "$canopus" > results/canopus_pp.tsv

echo "

classification_extract_sum have been successfully written into <results/canopus_pp.tsv>"
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
	 	until [ -d $projectpath ] && [ -f $projectpath/canopus.tsv ] && [ -f $projectpath/.format ]
	 	do
	 	read -p "Please input the path of the sirius project >>> " projectpath;
	 	done;
	 	cd $projectpath;
	fi;	
#########################
echo "

In order to filter out unnecessary categories, it is necessary to filter by the given list.
You can also customize this list by simply including keywords which separated by space.";

inputlist=0; 
	until [[ $inputlist == "most specific class" ]] || [[ $inputlist == "superclass" ]] || [[ $inputlist == "customlist" ]]
	do
	read -p "Please select the desired list.[most specific class/superclass/customlist] >>> " inputlist;
	done;
	
	if [[ $inputlist == "most specific class" ]];
	then index="results/spe_class.tsv"
	
	elif [[ $inputlist == "superclass" ]];
	then index="results/super_class.tsv"
	
	elif [[ $inputlist == "customlist" ]];
	then
	 	until [[ $the_list != "" ]]
	 	do
	 	read -p "Please enter the classification specified. For example, you can input [lignan iridoid steriod] >>> " the_list;
	 	done;
	else echo "error"
	fi;
	
	
	
	
	
	
	
check=0
	until [[ "$check" == "yes" ]] || [[ "$check" == "no" ]]
	do
	read -p "Filter the classification througn min Standard Deviation limition ?" check;
	done;
	
	if [[ $check == "yes" ]]
	then
	
	 	until [[ "$sigma" > "0.01" ]] && [[ "$sigma" < "0.5" ]]
	 	do
	 	read -p "Please enter the sigma limition. 0.1~0.3 may work well. >>> " l_limit;
	 	done;
	
	data="results/canopus_pp.tsv"
	
	sigma=0.2
	
	filter=$(awk -F $'\t' -v OFS=$'\t' '
 	 	{
 	 	if(NR==1)
 	 	 	{
 	 	 	for(i=2; i<=NF; i++)
 	 	 	 	{
 	 	 	 	class[i]=$i
 	 	 	 	}
 	 	 	}
 	 	if(NR>=2)
 	 	 	{
 	 	 	for(i=2; i<=NF; i++)
 	 	 	 	{
 	 	 	 	id[NR]=$1
 	 	 	 	\
 	 	 	 	pp[NR,i]=$i
 	 	 	 	\
 	 	 	 	sum[i]+=$i
 	 	 	 	}
 	 	 	}
 	 	}
 	 	END{
 	 	 	for(i=2; i<=NF; i++)
 	 	 	 	{
 	 	 	 	u[i]=sum[i]/(NR-1)
 	 	 	 	\
 	 	 	 	for(j=2; j<=NR; j++)
 	 	 	 	 	{
 	 	 	 	 	S[i]+=(pp[j,i]-u[i])^2
 	 	 	 	 	}
 	 	 	 	sigma[i]=sqrt(S[i]/(NR-1))
 	 	 	 	\
 	 	 	 	if(sigma[i]<'$sigma')
 	 	 	 	 	{
 	 	 	 	 	delete class[i]
 	 	 	 	 	}
 	 	 	 	}
 	 	 	printf "id\t"
 	 	 	
 	 	 	for(i=2; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if(class[i]!="")
 	 	 	 	 	{
 	 	 	 	 	printf class[i]"\t"
 	 	 	 	 	}
 	 	 	 	}
 	 	 	printf "\n"
 	 	 	
 	 	 	for(i=2; i<=NR; i++)
 	 	 	 	{
 	 	 	 	printf id[i]"\t"
 	 	 	 	
 	 	 	 	for(j=2; j<=NF; j++)
 	 	 	 	 	{
 	 	 	 	 	if(class[j]!="")
 	 	 	 	 	 	{
 	 	 	 	 	 	printf pp[i,j]"\t"
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	printf "\n"
 	 	 	 	}
 	 	 	}' $data)
 	 	 	
 	echo "$filter" > results/canopus_pp_filter.tsv
 	
 	echo "classification_extract_filter have been successfully written into <results/canopus_pp_filter.tsv>"
 	
 	fi;
;;
###################################
###################################
###################################
###################################
fragment_tree_network)
echo "Run fragment_tree_network.";
echo "Fragment_tree_network attempts to build molecular networks between different features with fragment tree alignment theory."
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
#if [ -d $projectpath ] && [ -f $projectpath/canopus_summary.tsv ]
check=0
	until [[ $check == "yes" ]]
	do
	read -p "Make sure you have moved the fragment tree alignment results to the directory where the project is located. [yes/no] >>> " check;
	done;
echo "Please enter the minimum alignment similarity that you want to filter."
tlimit=0
	until [[ "$tlimit" > "0.3" ]] && [[ "$tlimit" < "0.9" ]]
	do
	read -p "0.4-0.7 is recommended >>> " tlimit;
	done;
mkdir temp/ftaligntemp;
id=$(awk 'BEGIN{NF==1}{print $1}' ftalign.tsv | awk -F _ '{print $NF}')

ftalign=$(paste <(echo "$id") <(cut -f 2- ftalign.tsv) | sed '1d')

re_ftalign=$(cat <(paste -s <(echo "$id")) <(echo "$ftalign"));

echo "Reformat the data..."

fta=$(
awk -F $'\t' -v OFS=$'\t' '
{for(i=1; i<=NF; i++){raw[NR,i]=$i}}
	END{
	for(a=2; a<=NR; a++)
	{
	 	for(b=2; b<=NF; b++)
	 	{
	 	norm1[a,b]=raw[a,b]/raw[a,a];
	 	\
	 	norm2[a,b]=raw[a,b]/raw[b,b];
	 	\
	 	norms[a,b]=(norm1[a,b]+norm2[a,b])/2;
	 	\
	 	 	if(norms[a,b]>='$tlimit' && norms[a,b]<=0.999)
	 	 	{
	 	 	 	if(raw[a,1]>=raw[1,b])
	 	 	 	{
	 	 	 	print raw[a,1],raw[1,b],norms[a,b]
	 	 	 	}
	 	 	 	else
	 	 	 	{
	 	 	 	print raw[1,b],raw[a,1],norms[a,b]
	 	 	 	}
	 	 	}
	 	 }
	 }
 	}' <(echo "$re_ftalign"))
sort <(echo "$fta") | uniq > temp/ftaligntemp/filter_net_$tlimit
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
plimit=0
	until [[ "$plimit" > "0.5" ]] && [[ "$plimit" < "0.99" ]]
	do
	read -p "Please enter the minimum posterior probability of the molecular fingerprint to be filtered. 0.9-0.99 is recommended. >>> " plimit;
	done;
###################
###################
plimit2=0
	until [[ "$plimit2" > "0.1" ]] && [[ "$plimit2" < "$plimit" ]]
	do
	read -p "The minimum posterior probability of the molecular fingerprint to be controled. 0.1-0.5 may work well. >>> " plimit2;
	done;
	
tlimit=$(ls temp/ftaligntemp/filter_net_* | awk -F _ '{print $NF}')

check_rep=$(awk 'END{print NR}' <(echo "$tlimit"))
	if [[ "$check_rep" > "1" ]]
	then
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
	
awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $3,$14}}' results/finfirst_set.tsv > temp/Mo_filename;

echo "Aquiring data from sirius index..."

database=$(cat *_*/compound.info)

echo "Run fragment_tree_network_delta."

data=$(awk -F ["\t":_"\]"] -v OFS=$'\t' '
 	BEGIN{
 	i=0;
 	}
 	{
 	if(NR==FNR)
 	 	{
 	 	if($1=="name")
 	 	 	{
 	 	 	i+=1;
 	 	 	id[i]=$NF;
 	 	 	}
 	 	if($1=="ionMass")
 	 	 	{
 	 	 	mz[id[i]]=$NF
 	 	 	}
 	 	if($1=="ionType")
 	 	 	{
 	 	 	type[id[i]]=$NF
 	 	 	}
 	 	if($1=="rt")
 	 	 	{
 	 	 	rt[id[i]]=$2;
 	 	 	}
 	 	}
 	if(NR!=FNR)
 	 	{
 	 	#source  target  ftalign  delta_mz  delta_rt  source_iontype  target_iontype;
 	 	\
 	 	print $1,$2,$3,sprintf("%.3f",mz[$2]-mz[$1]),sprintf("%.2f",rt[$2]-rt[$1]),type[$1],type[$2]
 	 	}
 	}' <(echo "$database") <(cat temp/ftaligntemp/filter_net_$tlimit))
 	
data_path=$(awk -F ["\t"_] -v OFS=$'\t' '
 	{
 	if(NR==FNR)
 	 	{
 	 	mo[$NF]=$1
 	 	}
 	if(NR!=FNR)
 	 	{
 	 	if(mo[$1]!="" && mo[$2]!="")
 	 	 	{
 	 	 	\
 	 	 	#<path>sourceFormula  <path>targetFormula
 	 	 	\
 	 	 	printf "*_" $1 "/fingerprints/" mo[$1] "*\n"   "*_" $2 "/fingerprints/" mo[$2] "*\n"
 	 	 	}
 	 	}
 	}' <(cat temp/Mo_filename) <(echo "$data") | sed '$d' | sort -u | awk '
 	{
 	if(system("test -f " $1))
 	 	{
 	 	printf ""
 	 	}
 	else
 	 	{
 	 	print $1
 	 	}
 	}' | paste -d " " -s)
 	
data_allfp=$(awk -F $'\t' -v OFS=$'\t' '
 	BEGIN{
 	 	n=0
 	 	}
 	{
 	if(FNR==1)
 	 	{
 	 	if(n>=1)
 	 	 	{
 	 	 	close(file)
 	 	 	}
 	 	file=FILENAME;
 	 	\
 	 	n+=1;
 	 	\
 	 	printf FILENAME"\n"$0"\n"
 	 	}
 	else
 	 	{
 	 	print $0
 	 	}
 	}' $data_path)
awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $2,$3}}' csi_fingerid.tsv > temp/ftaligntemp/pos_fingerprint_index;

awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $2,$3}}' csi_fingerid_neg.tsv > temp/ftaligntemp/neg_fingerprint_index;

posindex=$(awk -F $'\t' '{print $1}' temp/ftaligntemp/pos_fingerprint_index)

negindex=$(awk -F $'\t' '{print $1}' temp/ftaligntemp/neg_fingerprint_index)
#
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
 	 	 	if(x>f)
 	 	 	 	{
 	 	 	 	f=x  # calculate the max index.
 	 	 	 	}
 	 	 	count+=1;  # calculate the all fingerprints file number.
 	 	 	\
 	 	 	split($1,a,"[/]"); e=split(a[1],b,"[_]"); id=b[e]; #catch the id.
 	 	 	\
 	 	 	x=0;
 	 	 	}
 	 	else
 	 	 	{
 	 	 	x+=1;
 	 	 	\
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
 	 	if(fp[$1,1]!="" && fp[$2,1]!="")
 	 	 	{
 	 	 	if($6=="+" && $7=="+")
 	 	 	 	{
 	 	 	 	for(x=1; x<=posnum; x++)
 	 	 	 	 	{
 	 	 	 	 	if(fp[$1,x]>='$plimit' && fp[$2,x]<='$plimit2')
	 	 	 	 	 	{
	 	 	 	 	 	data_s[FNR,x]=pos[x]
	 	 	 	 	 	}
	 	 	 	 	else if(fp[$2,x]>='$plimit' && fp[$1,x]<='$plimit2')
	 	 	 	 	 	{
	 	 	 	 	 	data_t[FNR,x]=pos[x]
	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	if($6=="-" && $7=="-")
 	 	 	 	{
 	 	 	 	for(x=1; x<=negnum; x++)
 	 	 	 	 	{
 	 	 	 	 	if(fp[$1,x]>='$plimit' && fp[$2,x]<='$plimit2')
	 	 	 	 	 	{
	 	 	 	 	 	data_s[FNR,x]=neg[x]
	 	 	 	 	 	}
	 	 	 	 	else if(fp[$2,x]>='$plimit' && fp[$1,x]<='$plimit2')
	 	 	 	 	 	{
	 	 	 	 	 	data_t[FNR,x]=neg[x]
	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	if($6=="-" && $7=="+" || $7=="-" && $6=="+")
 	 	 	 	{
 	 	 	 	for(i=1; i<=posnum; i++)
 	 	 	 	 	{
 	 	 	 	 	for(j=1; j<=negnum; j++)
 	 	 	 	 	 	{
 	 	 	 	 	 	if(pos[i]==neg[j])
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	mirror[i]=j 	#if pos=i, the identical index of neg is mirror[i]
 	 	 	 	 	 	 	}
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	if($6=="-" && $7=="+")	
 	 	 	 	{
 	 	 	 	for(x=1; x<=f; x++)
 	 	 	 	 	{
 	 	 	 	 	if(fp[$1,mirror[x]]>='$plimit' && fp[$2,x]<='$plimit2')
 	 	 	 	 	 	{
 	 	 	 	 	 	data_s[FNR,x]=neg[mirror[x]]
 	 	 	 	 	 	}
 	 	 	 	 	else if(fp[$2,x]>='$plimit' && fp[$1,mirror[x]]<='$plimit2')
 	 	 	 	 	 	{
 	 	 	 	 	 	data_t[FNR,x]=neg[mirror[x]]
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	if($7=="-" && $6=="+")
 	 	 	 	{
 	 	 	 	for(x=1; x<=f; x++)
 	 	 	 	 	{
 	 	 	 	 	if(fp[$1,x] >= '$plimit' && fp[$2,mirror[x]]<='$plimit2')
 	 	 	 	 	 	{
 	 	 	 	 	 	data_s[FNR,x]=neg[mirror[x]]
 	 	 	 	 	 	}
 	 	 	 	 	else if(fp[$2,mirror[x]]>='$plimit' && fp[$1,x]<='$plimit2')
 	 	 	 	 	 	{
 	 	 	 	 	 	data_t[FNR,x]=neg[mirror[x]]
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
 	 	 	\
 	 	 	if(fp[$1,1]=="" || fp[$2,1]=="")
 	 	 	 	{
 	 	 	 	printf "NA@NA\n"
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	printf "source:"  #the source fingerprint uniq.
 	 	 	 	\
 	 	 	 	for(x=1; x<=f; x++)
 	 	 	 	 	{
 	 	 	 	 	if(data_s[FNR,x]!="")
 	 	 	 	 	 	{
 	 	 	 	 	 	printf data_s[FNR,x]","
 	 	 	 	 	 	}
 	 	 	 	 	};
 	 	 	 	printf "@target:"  #the target fingerprint uniq.
 	 	 	 	\
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
 	 	}
 	}' <(echo "$data_allfp") <(echo "$posindex") <(echo "$negindex") <(echo "$data") | sed -e 's/,@/\t/g; s/,$//g; s/@/\t/g' \
 	 \
 	 > results/source_target_tree_$tlimit.tsv;
 	 
	echo "All instances have written into <results/source_target_tree_$tlimit.tsv>."
;;
###################################
###################################
###################################
###################################
double_ion_network)
echo "Run double_ion_network."

echo "It attempts to establish a link between features in different ionic modes based on retention times and molecular fingerprints. Please make sure that you have completed both ion modes using SIRIUS calculations and that you have completed steps such as fragment_tree_network above."

pos_projectpath=0
	until [ -d $pos_projectpath ] && [ -f $pos_projectpath/canopus_summary.tsv ]
	do
	read -p "Please input the path of the pos_ion project. >>> " pos_projectpath;
	done;
echo "pos_projectpath has successfully entered."
neg_projectpath=0
	until [ -d $neg_projectpath ] && [ -f $neg_projectpath/canopus_summary.tsv ] 
	do
	read -p "Please input the path of the neg_ion project. >>> " neg_projectpath;
	done;
echo "neg_projectpath has successfully entered."
save_path=0
	until [ -d $save_path ] 
	do
	read -p "Please input the storage path. >>> " save_path;
	done;
mkdir $save_path/double_ion_network

cd $save_path/double_ion_network;
mkdir temp
mkdir results
plimit="0"
	until [[ "$plimit" > "0.5" ]] && [[ "$plimit" < "0.99" ]]
	do
	read -p "Please enter the minimum posterior probability of the molecular fingerprint to be filtered. 0.9-0.99 is recommended. >>> "
	done;
flimit="0"
	until [[ "$flimit" > "0" ]] && [[ "$flimit" < "1000" ]]
	do
	read -p "Please enter the minimum identical fingerprints to be filtered. >>> " flimit;
	done;
dtime="0"
	until [[ "$dtime" > "0.01" ]] && [[ "$dtime" < "0.30" ]]
	do
	read -p "Please enter the rententiontime tolerance. 0.1-0.15 is recommended (C18). >>> " dtime;
	done;
dmass="0"
	until [[ "$dmass" > "0" ]] && [[ "$dmass" < "100" ]]
	do
	read -p "Please enter the Max delta_ionMass. 5-50 may work well. >>> " flimit;
	done;
pos_source=$(cat $pos_projectpath/temp/Mo_filename)
a=$(awk 'END{print NR}' <(echo "$pos_source"))
neg_target=$(cat $neg_projectpath/temp/Mo_filename)
b=$(awk 'END{print NR}' <(echo "$neg_source"))
awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $2,$3}}' $pos_projectpath/csi_fingerid.tsv > temp/pos_fingerprint_index;

awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $2,$3}}' $neg_projectpath/csi_fingerid_neg.tsv > temp/neg_fingerprint_index;
	if [ -f temp/doubleion_net_$plimit_$flimit ]
	then rm temp/doubleion_net_$plimit_$flimit
	fi;
	for (( i=1; i <= a; i++))
	do
	sumname=$(awk '{if(NR==$i){print}}' <(echo "$pos_source"))
	moname=$(echo "$sumname" | awk -F $'\t' '{print $1}');
	filename=$(echo "$sumname" | awk -F $'\t' '{print $2}');
	x=$(echo $pos_projectpath/$filename/fingerprints/$moname* )
	 	if [[ x != "$pos_projectpath/$filename/fingerprints/$moname*" ]]
	 	then
	id=$(echo "$filename" | awk -F _ '{print $NF}');
	source_mass=$(awk '/ionMass/{if(NR==FNR){print $2}}' $pos_projectpath/$filename/compound.info);
	source_rt=$(awk -F ['\t':] '/^rt/{if(NR==FNR){print $2/60}}' $pos_projectpath/$filename/compound.info)
	sourcefp=$(paste temp/pos_fingerprint_index $pos_projectpath/$filename/fingerprints/$moname* | awk -F $'\t' -v OFS=$'\t' '$3>='$plimit'{print $1}')
	##############
	#(parallel jobs)
	parallel_lists=$(seq 1 $b | paste -d " " - - - - - -)
	p=$(awk 'END{print NR}' <(echo "$parallel_lists"))
	 	for ((e=1; e <= p; e++))
	 	do
	 	f=$(awk '{if(NR=='$e'){print}}' <(echo "$parallel_lists"))
	 	 	for i in $(echo "$f")
	 	 	do
	 	 	{ t_sumname=$(awk '{if(NR==$i){print}}' <(echo "$neg_source"))
	 	 	t_filename=$(echo "$t_sumname" | awk -F $'\t' '{print $2}');
	 	 	t_moname=$(echo "$t_sumname" | awk -F $'\t' '{print $1}');
	 	 	y=$(echo $neg_projectpath/$t_filename/fingerprints/$t_moname*)
	 	 	 	if [[ $y != "$neg_projectpath/$t_filename/fingerprints/$t_moname*" ]]
	 	 	 	then
	 	 	 	target_mass=$(awk '/ionMass/{if(NR==FNR){print $2}}' $neg_projectpath/$t_filename/compound.info);
	 	 	 	target_rt=$(awk -F ['\t':] '/^rt/{if(NR==FNR){print $2/60}}' $neg_projectpath/$t_filename/compound.info);
	 	 	 	deltart=$(echo | awk '{print "'$source_rt'"-"'$target_rt'"}')
	 	 	 	deltamass=$(echo | awk '{print "'$source_mass'"-"'$target_mass'"}')
	 	 	 	 	if [[ "$deltart" < "$dtime" ]] && [[ "$deltart" > "0" ]] || [[ "$deltart" > "-$dtime" ]] && [[ "$deltart" < "0" ]]
	 	 	 	 	then
	 	 	 	 	 	if [[ "$deltamass" < "dmass" ]] && [[ "$deltamass" > "0" ]] || [[ "$deltamass" > "-$dmass" ]] && [[ "$deltamass" < "0" ]]
	 	 	 	 	 	then
	 	 	 	 	 	targetfp=$(paste temp/neg_fingerprint_index $neg_projectpath/$t_filename/fingerprints/$t_moname* | awk -F $'\t' -v OFS=$'\t' '$3>='$plimit'{print $1}')
	 	 	 	 	 	assemble=$(sort <(sort -m <(echo "$sourcefp") <(echo "$targetfp")) | uniq -d)
	 	 	 	 	 	fpnum=$(awk 'END{print NR}' <(echo "$assemble"))
	 	 	 	 	 	 	if [[ "$fpnum" > "$flimit" ]]
	 	 	 	 	 	 	then
	 	 	 	 	 	 	t_id=$(echo "$t_filename" | awk -F _ '{print $NF}');
	 	 	 	 	 	 	identicalfp=$(paste -s -d $'\t' <(echo "$assemble") | awk -F $'\t' -v OFS=, '{print $0}') 
	 	 	 	 	 	 	awk -F $'\t' -v OFS=$'\t' '{print $0,"'$t_id'","null","'$deltamass'","null","null","'$identicalfp'"}' <(echo "$id") > temp/single_doubleion_net_${id}_${t_id}_$plimit_$flimit;
	 	 	 	 	 	 	fi;
	 	 	 	 	 	fi;
	 	 	 	 	fi;
	 	 	 	fi;
	 	 	}&
	 	 	done;
	 	 	wait;
	 	done;
	 	fi;
	echo "Parallel jobs have finished $i($a)."
	done;
cat temp/single_doubleion_net_*_$plimit_$flimit > temp/doubleion_net_$plimit_$flimit
cat <(echo source$'\t'target$'\t'tfalign$'\t'delta_mass$'\t'fp_source_uniq$'\t'fp_target_uniq$'\t'identicalfp) temp/doubleion_net_$plimit_$flimit > results/doubleion_net_$plimit_$flimit.tsv;
echo "All instances have written into <results/doubleion_net_$plimit_$flimit.tsv>."
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

