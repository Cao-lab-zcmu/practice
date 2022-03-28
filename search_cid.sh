
filename=$(cat list | awk '{name="/media/wizard/back/0703_all/*_"$0"/fingerid/*.tsv"; printf name" "}')

grep "445858" $filename
