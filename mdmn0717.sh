##########################
##########################
echo "We are all in the gutter, 
but some of us are looking at the stars.";

PS3='Please select the workflow to be executed. >>> '

select command in \
  	\
 	"default" \
 	\
 	"structure_extract" \
 	\
 	"classification_extract_sum" \
 	\
 	"classification_extract_filter" \
 	\
 	"fragment_tree_network" \
 	\
 	"fragment_tree_network_delta" \
 	\
 	"compound_idenfication" \
 	\
 	"double_ion_network" \
 	\
 	"exit" 
do
	
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
 	 	\
 	 	max_id=-1
 	 	\
 	 	p_id="null"
 	 	\
 	 	printf "Info: loading the data..."
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
 	 	 	 	if($i~/links/)
 	 	 	 	 	{
 	 	 	 	 	col_links=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	pfile=FILENAME
 	 	
 	 	printf "Info: the data of " id " (ID) has input. The filename is " FILENAME ".\n"
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
 	 	 	
 	 	 	links["score",id]=$col_links
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
 	 	 	
 	 	 	score["simi",id]=$col_score
 	 	 
 	 	 	xlogp["simi",id]=$col_x
 	 	 	
 	 	 	links["simi",id]=$col_links
 	 	 	}
 	 	}
 	}
 	END{
 	 	printf "id\t"  "name\t"  "formula\t"  "similarity\t"  "smiles\t"  "inchi\t"  "inchikey2D\t" \
 	 	\
 	 	"score\t"  "xlogp\t"  "links\n" > "results/fingerid_sum.tsv"
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
 	 	 	 	score["score",l]"\t"  xlogp["score",l]"\t"  links["score",l]"\n" >> "results/fingerid_sum.tsv"
 	 	 	 	
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
 	 	 	 	score["simi",l]"\t"  xlogp["simi",l]"\t"  links["simi",l]"\n" >> "results/fingerid_sum.tsv"
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
 	 	\
 	 	if(max_id<id || max_id=="")
 	 	 	{
 	 	 	max_id=id
 	 	 	}
 	 	if(score[id]<$col_score || score[id]=="")
 	 	 	{
 	 	 	score[id]=$col_score
 	 	 	\
 	 	 	data[id]=$0
 	 	 	}
 	 	}
 	}
 	END{
 	 	printf "id\t"  "name\t"  "formula\t"  "similarity\t"  "smiles\t"  "inchi\t"  "inchikey2D\t" \
 	 	\
 	 	"score\t"  "xlogp\t"  "links\n" > "results/fingerid_first_score.tsv"
 	 	\
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
 	if(system("test -f " $2 "/canopus/" $1 "*.fpt"))
 	 	{
 	 	printf ""
 	 	}
 	else
 	 	{
 	 	printf $2 "/canopus/" $1 "*.fpt "
 	 	}
 	}' $data1)
	
awk -F $'\t' '
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
 	 	 	
 	 	 	printf "Info: data_file name of '$data1' has been input.\n"
 	 	 	
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
 	 	 	
 	 	 	printf "Info: data_file name of '$data2' has been input.\n"
 	 	 	
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
 	 	 	 	
 	 	 	 	printf "Info: data_file name of " p_filename " has been input.\n"
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
 	 	printf "id\t" > "results/canopus_pp.tsv"
 	 	\
 	 	for(i=0; i<=maxindex; i++)
 	 	 	{
 	 	 	if(class[1,i]!="" && i!=maxindex)
 	 	 	 	{
 	 	 	 	printf class[1,i]"\t" >> "results/canopus_pp.tsv"
 	 	 	 	}
 	 	 	else if(class[2,i]!="" && i!=maxindex)
 	 	 	 	{
 	 	 	 	printf class[2,i]"\t" >> "results/canopus_pp.tsv"
 	 	 	 	}
 	 	 	else if(i==maxindex && class[1,i]!="")
 	 	 	 	{
 	 	 	 	printf class[1,i]"\n" >> "results/canopus_pp.tsv"
 	 	 	 	}
 	 	 	else if(i==maxindex && class[2,i]!="")
 	 	 	 	{
 	 	 	 	printf class[2,i]"\n" >> "results/canopus_pp.tsv"
 	 	 	 	}
 	 	 	}
 	 	for(i=1; i<=n; i++)
 	 	 	{
 	 	 	printf the_id[i]"\t" >> "results/canopus_pp.tsv"
 	 	 	
 	 	 	for(j=0; j<=maxindex; j++)
 	 	 	 	{
 	 	 	 	if(class[1,j]!="" && j!=maxindex || class[2,j]!="" && j!=maxindex)
 	 	 	 	 	{
 	 	 	 	 	if(pp[the_id[i],j]=="")
 	 	 	 	 	 	{
 	 	 	 	 	 	printf 0"\t" >> "results/canopus_pp.tsv"
 	 	 	 	 	 	}
 	 	 	 	 	else
 	 	 	 	 	 	{
 	 	 	 	 	 	printf pp[the_id[i],j]"\t" >> "results/canopus_pp.tsv"
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	else if(class[1,j]!="" && j==maxindex || class[2,j]!="" && j==maxindex)
 	 	 	 	 	{
 	 	 	 	 	if(pp[the_id[i],j]=="")
 	 	 	 	 	 	{
 	 	 	 	 	 	printf 0"\n" >> "results/canopus_pp.tsv"
 	 	 	 	 	 	}
 	 	 	 	 	else
 	 	 	 	 	 	{
 	 	 	 	 	 	printf pp[the_id[i],j]"\n" >> "results/canopus_pp.tsv"
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}' $data1 $data2 $data3 $datas

echo "classification_extract_sum have been successfully written into <results/canopus_pp.tsv>"
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
	 	read -p "Please enter the posterior probabilities limition. 0.5~0.99 may work well. >>> " definition_limit;
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
	 	 	 	"subclass\t"  "subclass_pp\t"  "class\t"  "class_pp\t"  "superclass\t"  "superclass_pp\n" \
	 	 	 	\
	 	 	 	> "results/stat_classification.tsv"
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
	 	 	 	 	\
	 	 	 	 	definition=level[$col_id]
	 	 	 	 	}
	 	 	 	else if(subclass[$col_id] != "" && c_pp[subclass[$col_id]] >= "'$definition_limit'")
	 	 	 	 	{
	 	 	 	 	definition_source="subclass"
	 	 	 	 	\
	 	 	 	 	definition=subclass[$col_id]
	 	 	 	 	}
	 	 	 	else if(class[$col_id] != "" && c_pp[class[$col_id]] >= "'$definition_limit'")
	 	 	 	 	{
	 	 	 	 	definition_source="class"
	 	 	 	 	\
	 	 	 	 	definition=class[$col_id]
	 	 	 	 	}
	 	 	 	else if(superclass[$col_id] != "" && c_pp[superclass[$col_id]] >= "'$definition_limit'")
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
	 	 	 	superclass[$col_id]"\t"  c_pp[superclass[$col_id]]"\n" >> "results/stat_classification.tsv"
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

echo "Please enter the minimum alignment similarity that you want to filter."

tlimit=0
	until [[ "$tlimit" > "0.1" ]] && [[ "$tlimit" < "0.9" ]]
	do
	read -p "0.4-0.7 is recommended >>> " tlimit;
	done;

mkdir temp/ftaligntemp

data=ftalign.tsv

ftalign=$(awk -F $'\t' -v OFS=$'\t' '
 	{
 	for(i=1; i<=NF; i++)
 	 	{
 	 	if(NR==1 && i!=1 || NR!=1 && i==1)
 	 	 	{
 	 	 	n=split($i,x,"[_]")
 	 	 	\
 	 	 	raw[NR,i]=x[n]
 	 	 	}
 	 	else
 	 	 	{
 	 	 	raw[NR,i]=$i
 	 	 	}
 	 	}
 	}
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
 	 	}' $data)
 	
sort <(echo "$ftalign") | uniq > temp/ftaligntemp/filter_net_$tlimit

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
	read -p "Please enter the minimum posterior probability of the molecular fingerprint to be filtered. 0.9-0.99 is recommended. >>> " plimit;
	done;

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

echo "Aquiring data from sirius index..."

data1="*_*/compound.info"

data2="temp/ftaligntemp/filter_net_$tlimit"

echo "Run fragment_tree_network_delta."

data=$(awk -F ["\t":_"\]"] -v OFS=$'\t' '
 	BEGIN{
 	 	printf "id\t"  "m/z\t"  "rt\n" > "results/mz_and_rt.tsv"
 	 	}
 	{
 	if(FILENAME~/compound.info/)
 	 	{
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
 	 	 	
 	 	 	printf sprintf("%.2f",$NF)"\n" >> "results/mz_and_rt.tsv"
 	 	 	}
 	 	}
 	if(FILENAME~/ftaligntemp/)
 	 	{
 	 	#source  target  ftalign  delta_mz  delta_rt  source_iontype  target_iontype;
 	 	
 	 	print $1,$2,$3,sprintf("%.3f",mz[$2]-mz[$1]),sprintf("%.2f",rt[$2]-rt[$1]),type[$1],type[$2]
 	 	}
 	}' $data1 $data2)
 	
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
 	}' temp/Mo_filename <(echo "$data") | sort -u | awk '
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
compound_idenfication)

echo "Run compound_idenfication."















echo "compound_idenfication "

exit

;;
###################################
###################################
double_ion_network)

echo "Run double_ion_network."

echo "Double_ion_network attempts to establish a link between features in different ionic modes based on retention times and molecular fingerprints. Please make sure that you have completed both ion modes using SIRIUS calculations and that you have completed steps such as fragment_tree_network above."

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

