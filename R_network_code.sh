mkdir r_network
mkdir r_network/temp
mkdir r_network/ms2_figures
mkdir r_network/structure_svg
mkdir r_network/frags_figures
mkdir r_network/annonation
mkdir r_network/nodes_figures
#############################
#############################
#############################
#####
##############
######for Unix
#edges
file="results/source_target_tree*"
data=$(awk '{if(($1 ~ /^[0-9]/)){print $0}}' $file)
################
################
################
awk -F ["\t":,] '{
 	if(NR==1)
 	 	{
 	 	printf "from\t"  "to\t"  "length\t"  "width\t"  "label\t"  "arrows\t"  "dashes\t"  "smooth\t"  "shadow\t"  "title\n"
 	 	}
 	if(NR>=2)
 	 	{
 	 	norm=sprintf("%.2f",$3);
 	 	\
 	 	norm2=sprintf("%.2f",$4);
 	 	\
 	 	len=150/$3;
 	 	\
 	 	width=10*$3;
 	 	\
 	 	printf $1"\t"  $2"\t"  len"\t"  width"\t"  norm2"\t"  "to\t"  "FALSE\t"  "FALSE\t"  "FALSE\t" \
 	 	\
 	 	"<b>ftalign_similarity:"norm"</b>"  "</br>"  "<b>delta_m/z:"norm2"</b>"  "</br>";\
 	 	\
 	 	for(i=5; i<=NF; i++)
 	 	 	{
 	 	 	if(i>=6 && $i=="NA")
 	 	 	 	{
 	 	 	 	NA=1
 	 	 	 	}
	 	 	if($i=="source")
 	 	 	 	{
 	 	 	 	a=i+1;
 	 	 	 	}
 	 	 	if($i=="target")
 	 	 	 	{
 	 	 	 	b=i;
 	 	 	 	\
 	 	 	 	c=i+1
 	 	 	 	}
 	 	 	};
 	 	if(NA==1)
 	 	 	{
 	 	 	printf "fp:NA\n";
 	 	 	}
 	 	if(NA!=1)
 	 	 	{
 	 	 	printf "</br>source_uniq_fp:</br>";
 	 	 	\
 	 	 	for(x=a; x<b; x++)
 	 	 	 	{
 	 	 	 	printf "<a href=+frags_figures|" $x ".svg+target=+_blank+>" $x "</a> ";
 	 	 	 	}
 	 	 	printf "</br>target_uniq_fp:</br>";
 	 	 	\
 	 	 	for(y=c; y<=NF; y++)
 	 	 	 	{
 	 	 	 	if(y!=NF)
 	 	 	 	 	{
 	 	 	 	 	printf "<a href=+frags_figures|" $y ".svg+target=+_blank+>" $y "</a> ";
 	 	 	 	 	}
 	 	 	 	else
 	 	 	 	 	{
 	 	 	 	 	printf "<a href=+frags_figures|" $y ".svg+target=+_blank+>" $y "</a>\n";
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}
 	}' <(echo "$data") | sed -e 's/|/\//g;s/+/\"/g;$d' > r_network/edges.tsv
###<a href="test.txt"target="_blank">点击打开本地文件</a>
############################
############################
############################
############################
############################
############################
############################
###to create nodes(table)
#sort sample
mzmine_path="../neg.csv"
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
 	 	 	printf $1 "\t" $2 "\t" $3 "\t"
 	 	 	\
 	 	 	for(i=1; i<l; i++)
 	 	 	 	{
 	 	 	 	printf b[i] "\t"
 	 	 	 	}
 	 	 	printf b[l] "\n"
 	 	 	}
 	 	if(NR>=2)
 	 	 	{
 	 	 	printf $1"\t" $2"\t" $3"\t"
 	 	 	\
 	 	 	for(i=1; i<l; i++)
 	 	 	 	{
 	 	 	 	printf $col_sample[b[i]] "\t"
 	 	 	 	}
 	 	 	printf $col_sample[b[l]] "\n"
 	 	 	}
 	 	}' $mzmine_path)
 	echo "$mzmine_data" > r_network/mzmine_table.tsv
 	fi;
mzmine="r_network/mzmine_table.tsv"
###########################
###########################
awk -F $'\t' '
 	{
 	if(NR==1)
 	 	{
 	 	printf "id\tlabel\tgroup\tvalue\tshape\timage\tshadow\ttitle\n"
 	 	}
 	if(NR>=2)
 	 	{
 	 	norm=sprintf("%.4f",$2);
 	 	\
 	 	normrt=sprintf("%.2f",$3);
 	 	\
 	 	printf $1"\t"  norm  "\tNA"  "\t5"  "\timage"  "\tnodes_figures|"$1".svg\t"  "FALSE\t"\
 	 	\
 	 	"<b>ID:"$1"</b>" "</br>" "m/z:"norm "</br>" "RT(min):"normrt "</br>"\
 	 	\
 	 	"<a href=+ms2_figures|"$1".svg+target=+_blank+>See_MS2_spectra</a>" "</br>"\
 	 	\
 	 	"<a href=+structure_svg|"$1".svg+target=+_blank+>See_Structure_svg</a>" "</br>"\
 	 	\
 	 	"<a href=+annonation|"$1"+target=+_blank+>See_annonation</a>"\
 	 	\
 	 	"\n"
 	 	}
 	}' $mzmine | sed -e 's/|/\//g;s/+/\"/g;$d' > r_network/nodes.tsv
####tooltips#<a href="test.txt"target="_blank">点击打开本地文件</a>
###########################
###########################
############################
###to create ms2 figure(table)
data1="temp/Mo_filename"
data2="*_*/spectrum.ms"
data3=$(awk -F $'\t' '
 	{
 	if(system("test -f " $2 "/spectra/" $1 "*.tsv"))
 	 	{
 	 	printf ""
 	 	}
 	else
 	 	{
 	 	printf $2 "/spectra/" $1 "*.tsv "
 	 	}
 	}' $data1)
# to   r_network/ms2_figures/<id>_<formula>.tsv
mkdir results/ms2_figures
savepath="results/ms2_figures/"
awk -F ["\t"" "] '
 	BEGIN{
 	 	file_number=0
 	 	}
 	{
 	if(FNR==1)
 	 	{
 	 	if(FILENAME!=p_filename)
 	 	 	{
 	 	 	file_number+=1
 	 	 	}
 	 	if(file_number>=2)
 	 	 	{
 	 	 	close(p_filename)
 	 	 	}
 	 	p_filename=FILENAME
 	 	\
 	 	if((FILENAME~/spectrum.ms/))
 	 	 	{
 	 	 	step=2
 	 	 	}
 	 	else if((FILENAME~/spectra/))
 	 	 	{
 	 	 	step=3
 	 	 	}
 	 	}
 	if(NR==FNR)
 	 	{
 	 	formu[FNR]=$1
 	 	\
 	 	file[FNR]=$2
 	 	\
 	 	n=split($2,a,"[_]")
 	 	\
 	 	id[FNR]=a[n]
 	 	\
 	 	rows_data1=FNR
 	 	}
 	else if(step==2)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	the_id=$2
 	 	 	start=0
 	 	 	}
 	 	if(($1~/ms2peaks/))
 	 	 	{
 	 	 	start=FNR
 	 	 	}
 	 	if(start!=0)
 	 	 	{
 	 	 	n_ms2=FNR-start
 	 	 	\
 	 	 	msms[the_id,"mz",n_ms2]=$1
 	 	 	\
 	 	 	msms[the_id,"inten",n_ms2]=$2
 	 	 	\
 	 	 	if(max_n=="" || max_n<n_ms2)
 	 	 	 	{
 	 	 	 	max_n=n_ms2
 	 	 	 	}
 	 	 	if(max_inten[the_id]=="" || max_inten[the_id]<$2)
 	 	 	 	{
 	 	 	 	max_inten[the_id]=$2
 	 	 	 	}
 	 	 	}	
 	 	}
 	else if(step==3)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	if(step3_pro>=2)
 	 	 	 	{
 	 	 	 	close("'$savepath'"  the_id  ".tsv")
 	 	 	 	}
 	 	 	for(i=1; i<=rows_data1; i++)
 	 	 	 	{
 	 	 	 	if((FILENAME ~ file[i]) && (FILENAME ~ formu[i]))
 	 	 	 	 	{
 	 	 	 	 	the_id=id[i]
 	 	 	 	 	\
 	 	 	 	 	formula=formu[i]
 	 	 	 	 	\
 	 	 	 	 	break
 	 	 	 	 	}
 	 	 	 	else x=n+1
 	 	 	 	}
 	 	 	if(x==rows_data1+1)
 	 	 	 	{
 	 	 	 	nextfile
 	 	 	 	}
 	 	 	step3_pro+=1
 	 	 	\
 	 	 	if(step3_pro==1)
 	 	 	 	{
 	 	 	 	FS="[\t]"
 	 	 	 	if(FNR==1)
 	 	 	 	 	{
 	 	 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	 	 	{
 	 	 	 	 	 	if(($i~/mz/))
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	col_mz=i
 	 	 	 	 	 	 	}
 	 	 	 	 	 	if(($i~/rel/))
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	col_rel=i
 	 	 	 	 	 	 	}
 	 	 	 	 	 	if($i~/formula/)
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	col_formula=i
 	 	 	 	 	 	 	}
 	 	 	 	 	 	if($i~/ionization/)
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	col_ion=i
 	 	 	 	 	 	 	}
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	printf "mz\t"  "rel.intensity\t"  "formula\t"  "ionization\n" > "'$savepath'"  the_id  ".tsv"
 	 	 	\
 	 	 	for(j=1; j<=max_n; j++)
 	 	 	 	{
 	 	 	 	if(msms[the_id,"mz",j]!="")
 	 	 	 	 	{
 	 	 	 	 	printf msms[the_id,"mz",j]"\t" \
 	 	 	 	 	\
 	 	 	 	 	(msms[the_id,"inten",j]/max_inten[the_id])*100 "\n" >> "'$savepath'"  the_id  ".tsv"
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	norm=$col_rel*(-1)
 	 	 	\
 	 	 	printf $col_mz"\t"  norm"\t"  $col_formula"\t"  $col_ion"\n" >> "'$savepath'"  the_id  ".tsv"
 	 	 	}
 	 	}
 	}' $data1 $data2 $data3
###########################
#create match point
mkdir results/ms2_figures_match
data="results/ms2_figures/*.tsv"
savepath="results/ms2_figures_match/"
mz_tolerance="0.01"
in_tolerance="20"
awk -F $'\t' '
 	{
 	if(FNR==1)
 	 	{
 	 	file_num+=1
 	 	n=split(FILENAME,a,"[/]")
 	 	split(a[n],b,"[.]")
 	 	id=b[1]
 	 	if(file_num==1)
 	 	 	{
 	 	 	p_id=id
 	 	 	p_file=FILENAME
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	if($i~/mz/)
 	 	 	 	 	{
 	 	 	 	 	col_mz=i
 	 	 	 	 	}
 	 	 	 	if($i~/intensity/)
 	 	 	 	 	{
 	 	 	 	 	col_intensity=i
 	 	 	 	 	}
 	 	 	 	if($i~/formula/)
 	 	 	 	 	{
 	 	 	 	 	col_formula=i
 	 	 	 	 	}
 	 	 	 	if($i~/ionization/)
 	 	 	 	 	{
 	 	 	 	 	col_ion=i
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	if(file_num>=2)
 	 	 	{
 	 	 	close(p_file)
 	 	 	for(i in neg_mz)
 	 	 	 	{
 	 	 	 	for(j in pos_mz)
 	 	 	 	 	{
 	 	 	 	 	if(neg_mz[i]>=pos_mz[j]-"'$mz_tolerance'" && neg_mz[i]<=pos_mz[j]+"'$mz_tolerance'")
 	 	 	 	 	 	{
 	 	 	 	 	 	if((-1)*neg_in[i]>=pos_in[j]-"'$in_tolerance'" && (-1)*neg_in[i]<=pos_in[j]+"'$in_tolerance'")
 	 	 	 	 	 	 	{
 	 	 	 	 	 	 	match_p[i]+=1
 	 	 	 	 	 	 	match_p[j]+=1
 	 	 	 	 	 	 	}
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	delete neg_mz
 	 	 	delete pos_mz
 	 	 	printf "mz\t"  "rel.intensity\t"  "match\t"  "formula\t"  "ionization\n" > "'$savepath'" p_id ".tsv"
 	 	 	for(i=2; i<=max_FNR[p_file]; i++)
 	 	 	 	{
 	 	 	 	if(match_p[p_file,i]=="")
	 	 	 	 	{
	 	 	 	 	match_p[p_file,i]=0
	 	 	 	 	}
 	 	 	 	printf mz[p_file,i]"\t" intensity[p_file,i]"\t" match_p[p_file,i]"\t" formula[p_file,i]"\t" ion[p_file,i]"\n" \
 	 	 	 	\
 	 	 	 	>> "'$savepath'" p_id ".tsv"
 	 	 	 	}
 	 	 	if(p_id==2664)
 	 	 	 	{
 	 	 	 	printf file_num"\n"
 	 	 	 	}
 	 	 	close("'$savepath'" p_id ".tsv")
 	 	 	p_id=id
 	 	 	p_file=FILENAME
 	 	 	}
 	 	}
 	if(FNR>=2)
 	 	{
 	 	max_FNR[FILENAME]=FNR
 	 	mz[FILENAME,FNR]=$col_mz
 	 	intensity[FILENAME,FNR]=$col_intensity
 	 	formula[FILENAME,FNR]=$col_formula
 	 	ion[FILENAME,FNR]=$col_ion
 	 	if($col_intensity>0)
 	 	 	{
 	 	 	pos_mz[FILENAME,FNR]=$col_mz
 	 	 	pos_in[FILENAME,FNR]=$col_intensity
 	 	 	}
 	 	else
 	 	 	{
 	 	 	neg_mz[FILENAME,FNR]=$col_mz
 	 	 	neg_in[FILENAME,FNR]=$col_intensity
 	 	 	}
 	 	}
 	}
 	END{
 	 	for(i in neg_mz)
 	 	 	{
 	 	 	for(j in pos_mz)
 	 	 	 	{
 	 	 	 	if(neg_mz[i]>=pos_mz[j]-"'$mz_tolerance'" && neg_mz[i]<=pos_mz[j]+"'$mz_tolerance'")
 	 	 	 	 	{
 	 	 	 	 	if((-1)*neg_in[i]>=pos_in[j]-"'$in_tolerance'" && (-1)*neg_in[i]<=pos_in[j]+"'$in_tolerance'")
 	 	 	 	 	 	{
 	 	 	 	 	 	match_p[i]+=1
 	 	 	 	 	 	match_p[j]+=1
 	 	 	 	 	 	}
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	delete neg_mz
 	 	delete pos_mz
 	 	printf "mz\t"  "rel.intensity\t"  "match\t"  "formula\t"  "ionization\n" > "'$savepath'" p_id ".tsv"
 	 	for(i=2; i<=max_FNR[p_file]; i++)
 	 	 	{
 	 	 	if(match_p[p_file,i]=="")
 	 	 	 	{
 	 	 	 	match_p[p_file,i]=0
 	 	 	 	}
 	 	 	printf mz[p_file,i]"\t" intensity[p_file,i]"\t" match_p[p_file,i]"\t" formula[p_file,i]"\t" ion[p_file,i]"\n" \
 	 	 	\
 	 	 	>> "'$savepath'" p_id ".tsv"
 	 	 	}
 	 	}' $data
###########################
########Supplementary label and annonate
mkdir results/ms2_figures_label
data1="results/re_neg_RT.tsv"
data2="results/fingerid_first_score.tsv"
data3="results/ms2_figures_match/*.tsv"
savepath="results/ms2_figures_label/"
distance="10"
plus="5"
awk -F $'\t' '
 	BEGIN{
 	 	srand()
 	 	}
 	{
 	if(FILENAME~/re_neg_RT/)
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
 	 	 	p_file=FILENAME
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	precursor[$col_id]=$col_mz
 	 	 	rt[$col_id]=$col_rt
 	 	 	}
 	 	}
 	if(FILENAME~/fingerid/)
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
 	 	 	p_file=FILENAME
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	similarity[$col_id]=$col_similarity
 	 	 	}
 	 	}
 	if(FILENAME~/ms2_figures/)
 	 	{
	 	if(FNR==1)
	 	 	{
	 	 	close(p_file)
	 	 	file_num+=1
	 	 	p_file=FILENAME
	 	 	n=split(FILENAME,a,"[/]")
	 	 	split(a[n],b,"[.]")
	 	 	id=b[1]
	 	 	printf $0"\t"  "y_label\t" "precursor m/z\t" "RT (min)\t" "similarity\n" \
	 	 	\
	 	 	> "'$savepath'"a[n]
	 	 	if(file_num==1)
	 	 	 	{
	 	 	 	for(i=1; i<=NF; i++)
	 	 	 	 	{
	 	 	 	 	if($i~/mz/)
	 	 	 	 	 	{
	 	 	 	 	 	col_mz=i
	 	 	 	 	 	}
	 	 	 	 	if($i~/intensity/)
	 	 	 	 	 	{
	 	 	 	 	 	col_intensity=i
	 	 	 	 	 	}
	 	 	 	 	}
	 	 	 	}
	 	 	}
	 	if(FNR>=2)
	 	 	{
	 	 	printf $0"\t" >> "'$savepath'"a[n]
	 	 	mz[FNR]=$col_mz
	 	 	intensity[FNR]=$col_intensity
	 	 	if(intensity[FNR]>0)
	 	 	 	{
	 	 	 	y_label[FNR]=intensity[FNR]+"'$distance'"
	 	 	 	}
	 	 	else
	 	 	 	{
	 	 	 	y_label[FNR]=intensity[FNR]-"'$distance'"
	 	 	 	}
	 	 	for(i=2; i>=1; i--)
	 	 	 	{
	 	 	 	if(mz[FNR-i]!="")
	 	 	 	 	{
	 	 	 	 	distance=((mz[FNR]-mz[FNR-i])^2+(intensity[FNR]-intensity[FNR-i])^2)^(1/2)
	 	 	 	 	if(distance<"'$distance'"*(1/(2*i)))
	 	 	 	 	 	{
	 	 	 	 	 	if(y_label[FNR]>0)
	 	 	 	 	 	 	{
	 	 	 	 	 	 	# y_label[FNR]+=("'$plus'"+int(16*rand()+16*rand()))
	 	 	 	 	 	 	y_label[FNR]=0
	 	 	 	 	 	 	}
	 	 	 	 	 	else
	 	 	 	 	 	 	{
	 	 	 	 	 	 	# y_label[FNR]-=("'$plus'"+int(16*rand()+16*rand()))
	 	 	 	 	 	 	y_label[FNR]=0
	 	 	 	 	 	 	}
	 	 	 	 	 	}
	 	 	 	 	}
	 	 	 	}
	 	 	printf y_label[FNR] >> "'$savepath'"a[n]
	 	 	if(FNR==2)
	 	 	 	{
	 	 	 	printf "\t"sprintf("%.4f",precursor[id])"\t"  sprintf("%.2f",rt[id])"\t" sprintf("%.2f",similarity[id]) \
	 	 	 	\
	 	 	 	>> "'$savepath'"a[n]
	 	 	 	}
	 	 	printf "\n" >> "'$savepath'"a[n]
	 	 	}
 	 	}
 	}' $data1 $data2 $data3
#plot the map ms2_figures
###########################
###########################
###########################
###table_for_R_paint
	if [ -f results/select_canopus_filter_step3.tsv ]
	then 
	cp results/select_canopus_filter_step3.tsv r_network/temp/select
	elif [ -f results/select_canopus_filter_step2.tsv ]
	then 
	cp results/select_canopus_filter_step2.tsv r_network/temp/select
	else 
	cp results/select_canopus_filter_step1.tsv r_network/temp/select
	fi;
mkdir r_network/nodes_figures
data="r_network/temp/select"
awk -F $'\t' -v OFS=$'\t' '
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
	 	id[NR]=$1;
	 	\
	 	N=NR-1;
	 	\
	 	printf "class\tgroup\tposterior_probability\n" > "r_network/nodes_figures/class_"$1".tsv"
	 	\
	 	for(i=2; i<=NF; i++)
	 	 	{
	 	 	if(i!=NF)
	 	 	 	{
	 	 	 	printf class[i]"\tNA\t"$i"\n" >> "r_network/nodes_figures/class_"$1".tsv"
	 	 	 	}
	 	 	else
	 	 	 	{
	 	 	 	printf class[i]"\tNA\t"$i >> "r_network/nodes_figures/class_"$1".tsv"
	 	 	 	\
	 	 	 	if(NR>=3)
 	 	 	 	 	{
 	 	 	 	 	close("r_network/nodes_figures/class_"id[N]".tsv")
 	 	 	 	 	}
	 	 	 	}
	 	 	}
	 	}
 	}' $data
###########################
###########################
##########################
######features area to file.
mzmine="r_network/mzmine_table.tsv"
awk -F $'\t' '
 	{
 	if(NR==1)
 	 	{
 	 	for(i=1;i<=NF;i++)
 	 	 	{
 	 	 	sample[i]=$i
 	 	 	}
 	 	}
 	if(NR>=2)
 	 	{
 	 	id[NR]=$1;
 	 	\
 	 	printf "number\tsample\tarea\tlog10_area\n" > "r_network/nodes_figures/area_"$1".tsv"
 	 	\
 	 	N=NR-1
 	 	\
 	 	for(i=4;i<=NF;i++)
 	 	 	{
 	 	 	n=i-3;
 	 	 	if($i!=0)
 	 	 	 	{
 	 	 	 	norm=sprintf("%.2f",log($i)/log(10))
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	norm=0
 	 	 	 	};
 	 	 	if(i!=NF)
 	 	 	 	{
 	 	 	 	printf n"\t"sample[i]"\t"$i"\t"norm"\n" >> "r_network/nodes_figures/area_"$1".tsv"
 	 	 	 	}
 	 	 	else
 	 	 	 	{
 	 	 	 	printf n"\t"sample[i]"\t"$i"\t"norm >> "r_network/nodes_figures/area_"$1".tsv"
 	 	 	 	if(NR>=3)
 	 	 	 	 	{
 	 	 	 	 	close("r_network/nodes_figures/area_"id[N]".tsv")
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}
 	}' $mzmine
#########################
#########################
pr=$(awk -F $'\t' -v OFS=$'\t' '
 	{
 	if(NR>=2)
 	 	{
 	 	class="r_network/nodes_figures/class_"$1".tsv"
 	 	\
 	 	area="r_network/nodes_figures/area_"$1".tsv"
 	 	\
 	 	if(system("test -f " class))
 	 	 	{
 	 	 	c=0
 	 	 	}
 	 	else
 	 	 	{
 	 	 	c=1
 	 	 	}
 	 	if(system("test -f " area))
 	 	 	{
 	 	 	a=0
 	 	 	}
 	 	else
 	 	 	{
 	 	 	a=1
 	 	 	}
 	 	printf $1  "_"  c  "_"  a  "@"
 	 	}
 	}' $mzmine)
####R nodes image class;area
script=$(echo '
	setwd("r_network/nodes_figures")
	args<-commandArgs(TRUE)
	#####################
 	if ( "tidyverse" %in% rownames(installed.packages())==FALSE)
 	{
	install.packages("tidyverse")
	}
	library(tidyverse)
	library(grid)
	n=args[1]
	seq <- as.character(unlist(strsplit(n, split="@")))
        #############
	for(sum in seq)
	 	{
	 	sp=as.character(unlist(strsplit(sum, split="_")))
	 	id=sp[1]
	 	checkfile=sp[2]
	 	checkfile2=sp[3]
	 	file=paste0("class_",id,".tsv")
	 	file2=paste0("area_",id,".tsv")
	 	savename=paste0(id,".svg")
	 	##########
	 	if(checkfile==1)
	 	 	{
	 	 	data <- read.csv(file=file,header=T,sep="\t")
	 	 	data <- data %>% gather(key = "observation", value="value", -c(1,2)) 
	 	 	nObsType <- nlevels(as.factor(data$observation))
	 	 	data <- data %>% arrange(group, class)
	 	 	data$id <- rep( seq(1, nrow(data)/nObsType), each=nObsType)
	 	 	label_data <- data %>% group_by(id, class) %>% summarize(tot=sum(value))
	 	 	number_of_bar <- nrow(label_data)
	 	 	angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
	 	 	label_data$hjust <- ifelse( angle < -90, 1, 0)
	 	 	label_data$angle <- ifelse(angle < -90, angle+180, angle)
	 	 	p <- 
	 	 	 ggplot(data) +
	  	 	 geom_bar(aes(x=as.factor(id), y=value, fill=observation), 
	  	 	          stat="identity", alpha=0.5
	  	 	          ) +
	  	 	 scale_fill_brewer(palette = "Paired") +
	  	 	 ylim(-1,1) +
	  	 	 theme_minimal() +
	  	 	 theme(
	    	 	       axis.text = element_blank(),
	    	 	       axis.title.y = element_blank(),
	    	 	       panel.grid = element_blank(),
	    	 	       legend.position = "none",
	  	 	      ) +
	  	 	labs(x = "", y = "")+
	  	 	guides(fill=guide_legend(title=NULL)) +
	  	 	coord_polar() +
	  	 	geom_text(data=label_data, aes(x=id, y=0, label=class, hjust=hjust),
	  	 	          color="black", alpha=0.6, fontface="bold", size=1, angle= label_data$angle, 
	  	 	          inherit.aes = FALSE )
	 	 	}
	 	if(checkfile2==1)
	 	 	{
	 	 	data <- read.csv(file=file2,header=T,sep="\t")
	 	 	num <- nrow(data)
	 	 	max <- 10
	 	 	l <- ggplot(data, aes(x = number+2*num, y = log10_area, fill = sample)) + 
	  	 	     geom_bar(stat = "identity", position=position_dodge(0.7), width = 0.5, alpha=0.5) +
	  	 	              scale_fill_brewer(palette = "Paired") +
  	   	 	     theme_minimal() +
	   	 	     theme(
	    	 	           axis.text = element_blank(),
	    	 	           axis.title = element_blank(),
	    	 	           panel.grid = element_blank(),
	    	 	           legend.position = "none",
	   	 	          ) +
  	   	 	     geom_text(data=data, aes(x=number+2*num, fontface="bold", y=0, label=sample),
  	   	 	               color="black", size=1, alpha=0.5
  	   	 	               )+
  	   	 	     geom_text(data=data, aes(x=number+2*num, y=max/2, fontface="bold",
  	   	 	               label=log10_area), color="black", size=1, alpha=0.5
  	   	 	               )+
  	   	 	     labs(x = "", y = "")+
  	   	 	     xlim(0,num*3.43) +
  	   	 	     ylim(0,max) +
  	   	 	     guides(fill=guide_legend(title=NULL)) +
  	   	 	     coord_polar(theta = "y")
  	 	 	}
######################
	 	svg(savename)
	 	grid.newpage()
	 	if(checkfile==1)
	 	 	{
	 	 	print(p)
	 	 	}
	 	vp <- viewport(x = 0.5, y = 0.51, width = 0.72, height = 0.72)
	 	pushViewport(vp)
	 	if(checkfile2==1)
	 	 	{
	 	 	print(l, vp = vp)
	 	 	}
	 	dev.off()
	 	}
	')
###################### bg = "transparent"
#################
Rscript <(echo "$script") $pr
#################
#################
#################
###a rough classification define;
data="results/select_canopus_filter_step3.tsv"
awk -F $'\t' -v OFS=$'\t' '{
	if(NR==1)
	 	{
	 	m=NR+1;
	 	max[m]=0;
	 	for(i=2; i<=NF; i++)
	 	 	{
	 	 	group[i]=$i;
	 	 	}
	 	}
	if(NR>=2)
	 	{
	 	id[NR]=$1;
	 	for(j=2; j<=NF; j++)
	 	 	{
	 	 	if(max[NR]<$j)
	 	 	 	{
	 	 	 	max[NR]=$j;
	 	 	 	class_d[NR]=group[j];
	 	 	 	}
	 	 	}
	 	m=NR+1;
	 	max[m]=0;
	 	}
	}
	END{
	for(a=1;a<=NR;a++)
	 	{
	 	if(a==1)
	 	 	{
	 	 	print "id","class","pp"
	 	 	};
	 	if(a>=2)
	 	 	{
	 	 	print id[a],class_d[a],max[a]
	 	 	}
	 	}
	}' $data > test.tsv
###############
###############
###############
climit="0.95"
numlimit="0.2"
mzmine="r_network/mzmine_table.tsv"
data1=results/canopus_pp_filter.tsv
data2=$(echo "$mzmine")
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
forboxplot=$(awk -F ["\t"@] -v OFS=$'\t' '
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
 	 	printf "class\t"  "id\t"  "log'$log'_delta_area\t"  "pro_to_raw\t"  "log'$log_to'_pro_to_raw\t"  "variety\t"  "number\t";
 	 	\
 	 	printf "rt\t"  "m/z\n"
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
 	 	 	 	 	delta_area=pro-raw
 	 	 	 	 	\
 	 	 	 	 	to_raw=pro/raw
 	 	 	 	 	\
 	 	 	 	 	norm_to_raw=log(to_raw)/log('$log_to')
 	 	 	 	 	\
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
 	 	 	 	 	printf class[i]"\t"  id[i,j]"\t"  norm_delta"\t"  to_raw"\t"  norm_to_raw"\t"  variety"\t"  num[i]"\t" \
 	 	 	 	 	\
 	 	 	 	 	rt[id[i,j]]"\t"  mz[id[i,j]]"\n"
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}' <(echo "$stat") $(echo "$data2"))
echo "$forboxplot" > boxplot.tsv
#Rscript	boxplot
Rscript=$(echo '
library(ggplot2)
source<-read.table(file="boxplot.tsv",header=T,sep="\t")
boxplot<-ggplot(source,aes(x=class,y=log2_pro_to_raw,fill=class))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot()+
  geom_dotplot(binaxis ="y",
               stackdir="center",
               position="jitter",
               stackratio=0.2,
               dotsize = 0.25)+
  theme_minimal() +
  theme(
       legend.position="none",
       axis.text.x = element_blank(),
       text = element_text(size=30),
       )+
  geom_text(data=source, aes(x=class, y=-4, label=class),
	    color="black", fontface="bold", size=5, angle= 90, 
	    inherit.aes = FALSE )+
  geom_text(data=source, aes(x=class, y=5, label=number),
	    color="black", fontface="bold", size=5, 
	    inherit.aes = FALSE ) + 
  stat_summary(fun="mean",geom="point",shape=23,size=2.5,fill="grey")
pdf("boxplot.pdf",width=25,height=15)
boxplot
dev.off()
')
##############
lignans_iridoids=$(awk -F $'\t' -v OFS=$'\t' '
 	{
 	if(NR==1)
 	 	{
 	 	print "id","rt","m/z","classification","variety","pro/raw"
 	 	}
 	if(($1~/lignan/) || ($1~/ridoid/))
 	 	{
 	 	print $2,$8,$9,$1,$6,$4
 	 	}
  	}' <(echo "$forboxplot")
)
data=results/fingerid_first_score.tsv
list=lignans_and_iridoids.tsv
com=$(awk -F $'\t' '
 	BEGIN{
 	 	rows=0
 	 	}	
 	{
 	if(NR==FNR)
 	 	{
 	 	if(FNR>=2)
 	 	 	{
 	 	 	rows+=1
 	 	 	\
 	 	 	id[FNR]=$13;
 	 	 	\
 	 	 	simi[FNR]=$11
 	 	 	\
 	 	 	name[FNR]=$6
 	 	 	\
 	 	 	formula[FNR]=$3
 	 	 	\
 	 	 	smiles[FNR]=$7
 	 	 	\
 	 	 	pubchem[FNR]=$10
 	 	 	}
 	 	}
 	if(NR!=FNR)
 	 	{
 	 	if(FNR==1)
 	 	 	{
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	printf $i"\t"
 	 	 	 	}
 	 	 	printf "similarity\t"  "Name\t"  "Formula\t"  "Smiles\n";
 	 	 	}
 	 	if(FNR>=2)
 	 	 	{
 	 	 	for(i=1; i<=NF; i++)
 	 	 	 	{
 	 	 	 	printf $i"\t"
 	 	 	 	}
 	 	 	for(j=2;j<=rows;j++)
 	 	 	 	{
 	 	 	 	if(id[j]==$1)
 	 	 	 	 	{
 	 	 	 	 	printf simi[j]"\t"  name[j]"\t"  formula[j]"\t"  smiles[j]"\n";
 	 	 	 	 	}
 	 	 	 	}
 	 	 	}
 	 	}
 	}' $data $list)
 echo "$com" > com_lignans_and_iridoids.tsv
