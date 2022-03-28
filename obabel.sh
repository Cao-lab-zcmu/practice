cd ~/operation/back/0703_all/results/structure_2d/
data="../../com_compound.tsv"
data_repel="pubchem.tsv"
awk -F $'\t' '
 	{
 	if(NR==FNR)
 	 	{
	 	if(FNR==1)
	 	 	{
	 	 	for(i=1; i<=NF; i++)
	 	 	 	{
	 	 	 	if($i~/^id$/)
	 	 	 	 	{
	 	 	 	 	col_id=i
	 	 	 	 	}
	 	 	 	if($i~/smiles/)
	 	 	 	 	{
	 	 	 	 	col_smilse=i
	 	 	 	 	}
	 	 	 	if($i~/links/)
	 	 	 	 	{
	 	 	 	 	col_links=i
	 	 	 	 	}
	 	 	 	}
	 	 	}
	 	if(FNR>=2)
	 	 	{
	 	 	# print FNR," >>> ",$col_id
	 	 	if($col_links~/^PubChem/)
	 	 	 	{
	 	 	 	split($col_links, a, "[:][(]||[,][ ]||[)][;]||[)]")
	 	 	 	cid[$col_id]=a[2]
	 	 	 	# print a[2]
	 	 	 	}
	 	 	}
	 	}
	if(NR>FNR)
	 	{
	 	if(FNR==1)
	 	 	{
	 	 	file=FILENAME
	 	 	}
	 	if(FNR>=2)
	 	 	{
	 	 	repel[$1]=$1
	 	 	}
	 	}
	}
 	END{ 	
 	 	close(file)
 	 	RS="|||"
 	 	printf "id\t"  "cid\n" > "cid_metadata.tsv"
 	 	for(i in cid)
 	 	 	{
 	 	 	n+=1
 	 	 	printf i"\t"  cid[i]"\n" >> "cid_metadata.tsv"
 	 	 	if(n==1)
 	 	 	 	{
 	 	 	 	cid_set=cid[i]
 	 	 	 	}
 	 	 	if(n>=2)
 	 	 	 	{
 	 	 	 	cid_set=cid_set","cid[i]
 	 	 	 	}
 	 	 	}
 	 	# print cid_set
 	 	printf "id\t"  "cid\t" "IUPACName\t" "canonical_smiles\t"  "isomeric_smiles\n" >> "pubchem.tsv"
 	 	for(i in cid)
 	 	 	{
 	 	 	if(repel[i]=="")
 	 	 	 	{
	 	 	 	"curl https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" \
	 	 	 	\
	 	 	 	cid[i] \
	 	 	 	\
	 	 	 	"/property/IUPACName,CanonicalSMILES,IsomericSMILES/CSV" \
	 	 	 	\
	 	 	 	|getline data;
	 	 	 	gsub(/(.*)\n/,"",data)
	 	 	 	gsub(/,/,"\t",data)
	 	 	 	printf i"\t" data"\n" >> "pubchem.tsv"
	 	 	 	close("curl https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" \
	 	 	 	\
	 	 	 	cid[i] \
	 	 	 	\
	 	 	 	"/property/CanonicalSMILES,IsomericSMILES/CSV")
	 	 	 	}
	 	 	}
 	 	}' $data $data_repel
 ######################
 data="pubchem.tsv"
 sort -r $data | uniq > re_pubchem.tsv
 ######################
 data="re_pubchem.tsv"
 awk -F $'\t' '
  	{
  	if(FNR==1)
  	 	{
  	 	for(i=1; i<=NF; i++)
  	 	 	{
  	 	 	if($i~/^id$/)
  	 	 	 	{
  	 	 	 	col_id=i
  	 	 	 	}
  	 	 	if($i~/isomeric/)
  	 	 	 	{
  	 	 	 	col_isomeric=i
  	 	 	 	}
  	 	 	}
  	 	}
  	if(FNR>=2)
  	 	{
  	 	print $col_isomeric
  	 	system("obabel -:"  $col_isomeric  " -osvg -O "  $col_id  ".svg")
  	 	close("obabel -:"  $col_isomeric  " -osvg -O "  $col_id  ".svg")
  	 	}
  	}' $data
 #####################
 sed -i 's/white/transparent/g' [0-9]*.svg
 #####################
 ls *.svg | awk '
  	{
  	if($0~/cairo/)
  	 	{
  	 	next
  	 	}
  	else
  	 	{
  	 	system("cairosvg "  $0  " -o cairo_"  $0)
  	 	close("cairosvg "  $0  " -o cairo_"  $0)
  	 	}
  	}'
 #####################
 #####################
 #####################
 #####################
 # select <- c(279, 458, 574, 1107, 1445, 2227, 2529, 2664, 2824, 3380, 3918, 3938)
 idset=279,458,574,1107,1445,2227,2529,2664,2824,3380,3918,3938
 mkdir smiles_draw
 #path=~/operation/back/0703_all
 data="../../com_lignans_and_iridoids.tsv"
 savepath="smiles_draw"
 awk -F $'\t' '
  	BEGIN{
  	 	n=split("'$idset'", a, "[,]")
  	 	print n, a[1], a[2]
  	 	}
  	{
  	if(FNR==1)
  	 	{
  	 	print FILENAME
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
  	 	for(i in a)
  	 	 	{
  	 	 	if($col_id==a[i])
  	 	 	 	{
  	 	 	 	print $col_id, $col_smiles
  	 	 	 	system("obabel -:\""  $col_smiles  "\" -osvg -O '$savepath'/"  $col_id  ".svg")
  	 	 	 	close("obabel -:\""  $col_smiles  "\" -osvg -O '$savepath'/"  $col_id  ".svg")
  	 	 	 	}
  	 	 	}
  	 	}
  	}' $data
 ###########################
 ###########################
 idset=279,458,574,1107,1445,2227,2529,2664,2824,3380,3918,3938
 data1="re_pubchem.tsv"
 data2="../../com_lignans_and_iridoids.tsv"
 savepath="smiles_draw"
 awk -F $'\t' '
  	BEGIN{
  	 	n=split("'$idset'", a, "[,]")
  	 	print n, a[1], a[2]
  	 	}
  	{
  	if(FNR==1)
  	 	{
  	 	for(i=1; i<=NF; i++)
  	 	 	{
  	 	 	if($i~/^id$/)
  	 	 	 	{
  	 	 	 	col_id=i
  	 	 	 	}
  	 	 	if($i~/isomeric/)
  	 	 	 	{
  	 	 	 	col_isomeric=i
  	 	 	 	}
  	 	 	if($i~/^smiles/)
  	 	 	 	{
  	 	 	 	col_smiles=i
  	 	 	 	}
  	 	 	}
  	 	}
  	if(FNR>=2 && FILENAME ~/pubchem/)
  	 	{
  	 	for(i in a)
  	 	 	{
  	 	 	if($col_id==a[i])
  	 	 	 	{
	 	  	 	iso[$col_id]=$col_isomeric
	 	  	 	}
	 	  	}
  	 	}
  	if(FNR>=2 && FILENAME ~/com_/)
  	 	{
  	 	for(i in a)
  	 	 	{
  	 	 	if($col_id==a[i])
  	 	 	 	{
	 	  	 	smiles[$col_id]=$col_smiles
	 	  	 	}
	 	  	}
  	 	}
  	}
  	END{
  	 	printf "id\t"  "smiles\t"  "isomeric\n" > "'$savepath'/structure_chemdraw.tsv"
  	 	for(i in a)
  	 	 	{
  	 	 	printf a[i]"\t"  smiles[a[i]]"\t"  iso[a[i]]"\n" > "'$savepath'/structure_chemdraw.tsv"
  	 	 	}
  	 	}' $data1 $data2
 #####################
 #####################
 #####################
 #####################
 #####################
 ## complementation get SDF data
 mkdir sdf
 data="../../com_compound.tsv"
 savepath="sdf/"
 awk -F $'\t' '
 	{
 	if(NR==FNR)
 	 	{
	 	if(FNR==1)
	 	 	{
	 	 	for(i=1; i<=NF; i++)
	 	 	 	{
	 	 	 	if($i~/^id$/)
	 	 	 	 	{
	 	 	 	 	col_id=i
	 	 	 	 	}
	 	 	 	if($i~/smiles/)
	 	 	 	 	{
	 	 	 	 	col_smilse=i
	 	 	 	 	}
	 	 	 	if($i~/links/)
	 	 	 	 	{
	 	 	 	 	col_links=i
	 	 	 	 	}
	 	 	 	}
	 	 	}
	 	if(FNR>=2)
	 	 	{
	 	 	# print FNR," >>> ",$col_id
	 	 	if($col_links~/^PubChem/)
	 	 	 	{
	 	 	 	split($col_links, a, "[:][(]||[,][ ]||[)][;]||[)]")
	 	 	 	cid[$col_id]=a[2]
	 	 	 	# print a[2]
	 	 	 	}
	 	 	}
	 	}
	}
 	END{	
 	 	for(i in cid)
 	 	 	{
 	 	 	if(repel[i]=="")
 	 	 	 	{
	 	 	 	system("curl https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" \
	 	 	 	\
	 	 	 	cid[i] \
	 	 	 	\
	 	 	 	"/SDF > '$savepath'" i ".SDF")
	 	 	 	close("curl https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" \
	 	 	 	\
	 	 	 	cid[i] \
	 	 	 	\
	 	 	 	"/SDF > '$savepath'" i ".SDF")
	 	 	 	}
	 	 	}
 	 	}' $data
 #####################
 sed -i '1c\ ' sdf/[0-9]*.SDF
 #####################
 ls sdf/[0-9]*.SDF | awk -F $'\t' '
  	{
  	system("obabel "  $0  " -iSDF -osvg -O "  $col_id  ".svg")
  	close("obabel "  $0  " -iSDF -osvg -O "  $col_id  ".svg")
  	}' 
 #####################
 sed -i -e 's/white/transparent/g; s/stroke-width="2.0"/stroke-width="4.0"/g;' sdf/[0-9]*SDF.svg
 #####################
 ls sdf/[0-9]*.svg | awk '
  	{
  	if($0~/cairo/)
  	 	{
  	 	next
  	 	}
  	else
  	 	{
  	 	system("cairosvg "  $0  " -o "  $0 ".cairo.svg")
  	 	close("cairosvg "  $0  " -o "  $0 ".cairo.svg")
  	 	}
  	}'
 #####################
 ##################### 
 #####################
 ## mavin molconvert
 #####################
  ls smiles_draw/[0-9]*_1.sdf | awk -F $'\t' '
  	{
  	split($0, a, "[_][1][.][s][d][f]")
  	system("obabel "  $0  " -iSDF -osvg -O "  a[1]  ".svg")
  	close("obabel "  $0  " -iSDF -osvg -O "  a[1]  ".svg")
  	print a[1]".svg"
  	}'
 #####################
 ### create name
 molconvert "name:common,all" -s "COC(=O)C1=COC(C2C1CC=C2CO)OC3C(C(C(C(O3)CO)O)O)O"

