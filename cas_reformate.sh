###############
data="cas"
savepath="cas_arrange.tsv"
awk -F "[:]||['][,][ ][']||['][,][ ][\"]||[\"][,][ ][\"]||[\"][,][ ][']" '
 	{
 	if($1~/^\[/)
 	 	{
 	 	name=$1
 	 	id_num+=1
 	 	id[id_num]=name
 	 	print id_num, id[id_num]
 	 	}
 	if($1~/BEGIN_compound/)
 	 	{
 	 	num[id_num]+=1
 	 	getline;
 	 	for(i=1; i<=NF; i++)
 	 	 	{
 	 	 	if($i ~ /[0-9](.*)-[0-9][0-9]-[0-9](.*)$/ || $i ~ /^CAS/)
	 	 	 	{
	 	 	 	data[id_num,num[id_num]]=$i
	 	 	 	break;
	 	 	 	}
 	 	 	}
 	 	}
 	}
 	END{
 	 	printf "number\t"  "name\t"  "cas\n" > "'$savepath'"
 	 	for(i=1; i<=id_num; i++)
 	 	 	{
 	 	 	printf i"\t"id[i] >> "'$savepath'"
 	 	 	print i"\t"id[i]
 	 	 	for(j=1; j<=num[i]; j++)
 	 	 	 	{
 	 	 	 	printf "\t"data[i,j] >> "'$savepath'"
 	 	 	 	}
 	 	 	printf "\n" >> "'$savepath'"
 	 	 	}
 	 	}' $data
##############
data="cas_arrange.tsv"
sed -i -e 's/\[//g; s/\]//g; s/{//g; s/}//g' $data
