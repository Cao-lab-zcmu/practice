##########################
##########################
echo "We are all in the gutter, 
but some of us are looking at the stars."; \
PS3='Please select the workflow to be executed. >>> '
select command in "default" "structure_extract" "classification_extract_sum" "classification_extract_select" "fragment_tree_network" "fragment_tree_network_delta" "double_ion_network" "exit"
do
	default="structure_extract classification_extract_sum classification_extract_select fragment_tree_network fragment_tree_network_delta"
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
	else list=$( echo $command )
	fi;\
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
ls -l | awk '/^d/ {print $NF}' | awk '/_/ {print}' > temp/namelist;
echo "The temporary directory and the output directory have been successfully created.";
echo "Running structure_extract.";
a=$(awk 'END{print NR}' temp/namelist)
parallel_lists=$(seq 1 $a | paste -d " " - - - - - - -)
p=$(awk 'END{print NR}' <(echo "$parallel_lists"))
	for ((e=1; e <= p; e++))
	do
	f=$(awk '{if(NR=='$e'){print}}' <(echo "$parallel_lists"))
	 	for i in $(echo "$f")
	 	do
	 	{ filename=$(awk '{if(NR=='$i'){print}}' temp/namelist)
	 	x=$(echo $filename/fingerid/*.tsv)
	 	 	if [[ $x != "$filename/fingerid/*.tsv" ]]
	 	 	then
	 	 	id=$(echo $filename | awk -F _ '{print $NF}');
	 	 	cat $filename/fingerid/*.tsv | sort -t $'\t' -n -k 5 -r | uniq | sed 's/\r//g' | awk -F $'\t' -v OFS=$'\t' '{if(NR==1){print $0,"id","filename"};if(NR>=2){print $0,"'$id'","'$filename'"}}' > temp/fintemp/sum$id;
	 	 	sed -n '1p;2p' temp/fintemp/sum$id > temp/fintemp/first$id;
	 	 	fi;
	 	}&
	 	done;
	 	wait;
	echo "Parallel jobs have finished $e($p)."
	done;
cat temp/fintemp/first* | sort -t $'\t' -n -k 5 -r | uniq > results/finfirst_set.tsv;
echo "structure_extract results have been successfully assembled into <results/finfirst_set.tsv>";
;;
###################################
###################################
###################################
###################################
classification_extract_sum)
echo "Run classification_extract_sum.";
	if [ -d "temp" ] && [ -f canopus.tsv ] && [ -f .format ]
	then
	echo "Project path acknowledged."
	else
	 	until [ -d $projectpath ] && [ -f $projectpath/canopus.tsv ] && [ -f $projectpath/.format ]
	 	do
	 	read -p "Please input the path of the sirius project >>> " projectpath;
	 	done;
	 	cd $projectpath;
	fi;	
mkdir temp/canotemp;
awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $3,$14}}' results/finfirst_set.tsv > temp/Mo_filename;
awk -F $'\t' '{if(NR>=2){print $1}}' canopus.tsv > temp/canotemp/pos_canoindex;
#pos_cano_title=$(awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $4}}' canopus.tsv | paste -s);
awk -F $'\t' '{if(NR>=2){print $1}}' canopus_neg.tsv > temp/canotemp/neg_canoindex;
#neg_cano_title=$(awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $4}}' canopus_neg.tsv | paste -s);
list=$(cat temp/Mo_filename);
echo "Note that here the CANOPUS output is extracted by the highest ranking of the structure identification results. Based on <temp/Mo_filename>."
echo "Clear generated files.";
	if [ -f temp/canotemp/allpos ] || [ -f temp/canotemp/allneg ]
	then
	rm temp/canotemp/all*; 
	fi
echo "CANOPUS files being extracted..."
a=$(awk 'END{print NR}' <(echo "$list"))
parallel_lists=$(seq 1 $a | paste -d " " - - - - - - - - - - - -)
p=$(awk 'END{print NR}' <(echo "$parallel_lists"))
	for ((e=1; e <= p; e++))
	do
	f=$(awk '{if(NR=='$e'){print}}' <(echo "$parallel_lists"))
	 	for i in $(echo "$f")
	 	do
	 	{ sumname=$(awk '{if(NR=='$i'){print}}' <(echo "$list"))
	 	moname=$(echo "$sumname" | awk -F $'\t' '{print $1}');
	 	filename=$(echo "$sumname" | awk -F $'\t' '{print $2}');
	 	x=$(echo $filename/canopus/$moname*)
	 	 	if [[ $x != "$filename/canopus/$moname*" ]]
	 	 	then
	 	 	id=$(echo $filename | awk -F _ '{print $NF}');
	 	 	iontype=$( grep ionType $filename/compound.info | awk -F ] '{print $NF}')
	 	 	 	if [[ "$iontype" == "+" ]]
	 	 	 	then canoindex="pos_canoindex"; ion="pos"
	 	 	 	else canoindex="neg_canoindex"; ion="neg"
	 	 	 	fi;
	 	 	paste -s temp/canotemp/$canoindex $filename/canopus/$moname* | sed 's/\r//g' | awk -F $'\t' -v OFS=$'\t' '{if(NR==1){print $0,"id"};if(NR>=2){print $0,"'$id'"}}' > temp/canotemp/single_cano$id;
 	 	 	fi;
 	 	}&
 	 	done;
 	 	wait;
 	echo "Parallel jobs have finished $e($p)."
	done;
cat temp/canotemp/single_cano* | sort | uniq > temp/canotemp/all$ion;
echo "CANOPUS files have extracted successfully."
###############################################
	if [[ "$iontype" == "+" ]]
	then indexname=$(awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $4}}' canopus.tsv | paste -s)
	else indexname=$(awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $4}}' canopus_neg.tsv | paste -s)
	fi;\
cat <(echo "$indexname") <(sort temp/canotemp/all* | uniq | sed '1d') > results/cano_"$ion"_pp.tsv;
echo "classification_extract_sum have been successfully written into <results/cano_"$ion"_pp.tsv>"
;;
###################################
###################################
###################################
###################################
classification_extract_select)
echo "Run classification_extract_select.";
	if [ -d temp ] && [ -f canopus.tsv ] && [ -f .format ]
	then
	echo "Project path acknowledged."
	else
	 	until [ -d $projectpath ] && [ -f $projectpath/canopus.tsv ] && [ -f $projectpath/.format ]
	 	do
	 	read -p "Please input the path of the sirius project >>> " projectpath;
	 	done;
	 	cd $projectpath;
	fi;	
awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $4}}' canopus_summary.tsv | sort -k 1 | uniq > results/spe_class.tsv;
awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $8 }}' canopus_summary.tsv | sort -k 1 | uniq > results/super_class.tsv;
#########################
echo "In order to filter out unnecessary categories, it is necessary to filter by the given list. Two lists are provided by default (depend on your idenfication results), but you can also customize this list by simply including keywords.";
inputlist=0; \
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
	 	until [ -f $listpath ] && [[ $listpath != "" ]]
	 	do
	 	read -p "Please enter the path of the customlist >>> " listpath;
	 	done;
	 	index=$(echo "$listpath")
	else echo "error"
	fi;
awk '{if(NR==1){print}}' results/cano_*_pp.tsv | sed 's/\t/\n/g' > temp/canotemp/searchindex;
#########################
	if [ -f temp/canotemp/line_n ]
	then rm temp/canotemp/line_n
	fi;\
IFS_def=$IFS;\
IFS=$'\n';\
	for line_n in $(cat $index)
	do
	grep -F -i -n "$line_n" temp/canotemp/searchindex >> temp/canotemp/line_n
	done;\
IFS=$IFS_def;
awknum=$(awk -F : -v OFS=$'\t' '{print $1}' temp/canotemp/line_n | sort -n)
	for numline in $(echo "$awknum")
	do
	awk -F $'\t' -v OFS=$'\t' '{print $'$numline'}' results/cano_*_pp.tsv > temp/canotemp/awk$numline
	done;
paste <(awk -F $'\t' '{if(NR==1){print "id"}if(NR>=2){print $NF}}' results/cano_*_pp.tsv) temp/canotemp/awk* > results/select_cano_pp.tsv
echo "classification_extract_select have been successfully written into <results/select_cano_pp.tsv>"
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
	 	norm2[a,b]=raw[a,b]/raw[b,b];
	 	norms[a,b]=(norm1[a,b]+norm2[a,b])/2;
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
echo "Run fragment_tree_network_delta."
echo "If you have a large amount of data (maybe more than 1000 features or more), please enter the super administrator user and set "ulimit" before running fragment_ tree_ network_ delta. For example, enter >>> ulimt - N 65536 "
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
awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $2,$3}}' csi_fingerid.tsv > temp/ftaligntemp/pos_fingerprint_index;
awk -F $'\t' -v OFS=$'\t' '{if(NR>=2){print $2,$3}}' csi_fingerid_neg.tsv > temp/ftaligntemp/neg_fingerprint_index;
echo "The following module attempts to compute the differential fingerprints of connected clusters based on the molecular fingerprint results of sirius 4."
plimit=0
	until [[ "$plimit" > "0.5" ]] && [[ "$plimit" < "0.99" ]]
	do
	read -p "Please enter the minimum posterior probability of the molecular fingerprint to be filtered. 0.9-0.99 is recommended. >>> " plimit;
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
a=$(awk 'END{print NR}' temp/ftaligntemp/filter_net_$tlimit)
parallel_lists=$(seq 1 $a | paste -d " " - - - - - -)
p=$(awk 'END{print NR}' <(echo "$parallel_lists"))
	for ((e=1; e <= p; e++))
	do
	f=$(awk '{if(NR=='$e'){print}}' <(echo "$parallel_lists"))
	 	for i in $(echo "$f")
	 	do
	 	{ source=$(awk '{if(NR=='$i'){print $1}}' temp/ftaligntemp/filter_net_$tlimit)
	 	sourceion=$(grep ionMass *_$source/compound.info | awk '{print $2}' | sed 's/[[:space:]]//g' | sed -n '1p;1q')
	 	target=$(awk '{if(NR=='$i'){print $2}}' temp/ftaligntemp/filter_net_$tlimit)
	 	targetion=$(grep ionMass *_$target/compound.info | awk '{print $2}' | sed 's/[[:space:]]//g' | sed -n '1p;1q')
	 	deltaion=$(echo | awk '{print "'$sourceion'"-"'$targetion'"}')
	 	iontype=$( grep ionType *_$source/compound.info | awk -F ] '{print $NF}')
	 	 	if [[ "$iontype" == "+" ]]
	 	 	then
	 	 	printindex=pos_fingerprint_index
	 	 	else
	 	 	printindex=neg_fingerprint_index
	 	 	fi;
	 	echo "$iontype" >> temp/ftaligntemp/iontype
	 	source_moname=$(grep ''_$source'$' temp/Mo_filename | awk -F $'\t' -v OFS=$'\t' '{print $1}')
	 	target_moname=$(grep ''_$target'$' temp/Mo_filename | awk -F $'\t' -v OFS=$'\t' '{print $1}')
	 	source_filename=$(grep ''_$source'$' temp/Mo_filename | awk -F $'\t' -v OFS=$'\t' '{print $2}')
	 	target_filename=$(grep ''_$target'$' temp/Mo_filename | awk -F $'\t' -v OFS=$'\t' '{print $2}')
	 	x=$(echo $source_filename/fingerprints/$source_moname*)
	 	y=$(echo $target_filename/fingerprints/$target_moname*) 
	 	 	if [[ $x != "$source_filename/fingerprints/$source_moname*" ]] && [[ $y != "$target_filename/fingerprints/$target_moname*" ]]
	 	 	then
	 	 	sourcefp=$(paste temp/ftaligntemp/$printindex *_$source/fingerprints/$source_moname* | awk -F $'\t' -v OFS=$'\t' ' $3>='$plimit' {print $1}');
	 	 	targetfp=$(paste temp/ftaligntemp/$printindex *_$target/fingerprints/$target_moname* | awk -F $'\t' -v OFS=$'\t' ' $3>='$plimit' {print $1}');
	 	 	source_line=$(sort <(sort -m <(echo "$sourcefp" | sort) <(echo "$targetfp" | sort) <(echo "$targetfp" | sort)) | uniq -u);
	 	 	source_uniq=$(echo "$source_line" | paste -s -d , | awk -F : -v OFS=: '{print "source",$0}');
	 	 	target_line=$(sort <(sort -m <(echo "$sourcefp" | sort) <(echo "$sourcefp" | sort) <(echo "$targetfp" | sort)) | uniq -u);
	 	 	target_uniq=$(echo "$target_line" | paste -s -d , | awk -F : -v OFS=: '{print "target",$0}');
	 	 	else
	 	 	source_uniq="NA";
	 	 	target_uniq="NA";
	 	 	source_line="NA";
	 	 	target_line="NA";
	 	 	fi;
	 	awk -F $'\t' -v OFS=$'\t' '{if(NR=="'$i'"){print $0,"'$deltaion'","'$source_uniq'","'$target_uniq'"}}' temp/ftaligntemp/filter_net_$tlimit > temp/ftaligntemp/single_refilter_net_${source}_${target}_$tlimit;
	 	cat <(echo "$source_line") <(echo "$target_line") > temp/ftaligntemp/single_fpsample_${source}_${target}_$tlimit;
	 	}&
	 	done;
	 	wait;
	echo "Parallel jobs have finished $e($p)."
	done;
cat temp/ftaligntemp/single_refilter_net_*_$tlimit > temp/ftaligntemp/refilter_net_$tlimit
cat temp/ftaligntemp/single_fpsample_*_$tlimit > temp/ftaligntemp/fpsample
echo "Reformat the output files..."
############
sort temp/ftaligntemp/fpsample | uniq | awk -F _ -v OFS=_ '{print "Frag",$0}' > temp/ftaligntemp/refor_fpsample;
st_sample=$(paste <(sed '$d' temp/ftaligntemp/refor_fpsample) <(sed '1d' temp/ftaligntemp/refor_fpsample));
loop_end=$(paste <(awk 'END{print}' temp/ftaligntemp/refor_fpsample) <(awk '{if(NR==1){print}}' temp/ftaligntemp/refor_fpsample));
loop_st_sample=$(cat <(echo "$st_sample") <(echo "$loop_end"));
############
cat <(echo source$'\t'target$'\t'tfalign$'\t'delta_mass$'\t'fp_source_uniq$'\t'fp_target_uniq) temp/ftaligntemp/refilter_net_$tlimit <(echo "$loop_st_sample") > results/source_target_tree_$tlimit.tsv;
iontype=$(cat temp/ftaligntemp/iontype | sort -u)
uniform=$(awk 'END{print NR}' <(echo "$iontype"))
	if [[ "$uniform" > "1" ]]
	then
	echo "The presence of both positive and negative ion modes in the project was detected. This may cause errors."
	else
	 	if [[ "$iontype" == "+" ]]
	 	then
	 	printindex=pos_fingerprint_index
	 	else
	 	printindex=neg_fingerprint_index
	 	fi;	
	sharedname=$(awk -F $'\t' -v OFS=$'\t' '{print "Frag_"$1}' temp/ftaligntemp/$printindex);
	revise=$(awk -F $'\t' -v OFS=$'\t' '{print $2}' temp/ftaligntemp/$printindex | awk -F " " -v OFS=$'\t' '{print $1,$2$3$4$5}');
	cat <(echo sharedname$'\t'SMARTS$'\t'complementary) <(paste <(echo "$sharedname") <(echo "$revise")) > results/cyto_"$printindex".tsv;
	echo "All instances have written into <results/source_target_tree_$tlimit.tsv>, <results/cyto_$printindex.tsv>."
	fi;
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
	fi;\
	for (( i=1; i <= a; i++))
	do
	sumname=$(awk '{if(NR==$i){print}}' <(echo "$pos_source"))
	moname=$(echo "$sumname" | awk -F $'\t' '{print $1}');
	filename=$(echo "$sumname" | awk -F $'\t' '{print $2}');
	x=$(echo $pos_projectpath/$filename/fingerprints/$moname* )
	 	if [[ x != "$pos_projectpath/$filename/fingerprints/$moname*" ]]
	 	then
	id=$(echo "$filename" | awk -F _ '{print $NF}');
	source_mass=$(grep ionMass $pos_projectpath/$filename/compound.info | awk '{print $2}' | sed 's/[[:space:]]//g');
	source_rt=$(grep rt $pos_projectpath/$filename/compound.info | awk '{print $2}' | awk -F : '{min=$1/60;{print min}}');
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
	 	 	 	target_mass=$(grep ionMass $neg_projectpath/$t_filename/compound.info | awk '{print $2}' | sed 's/[[:space:]]//g');
	 	 	 	target_rt=$(grep rt $neg_projectpath/$t_filename/compound.info | awk '{print $2}' | awk -F : '{min=$1/60;{print min}}');\
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
cat <(echo source$'\t'target$'\t'tfalign$'\t'delta_mass$'\t'fp_source_uniq$'\t'fp_target_uniq$'\t'identicalfp) temp/doubleion_net_$plimit_$flimit > results/doubleion_net_$plimit_$flimit.tsv;\
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

