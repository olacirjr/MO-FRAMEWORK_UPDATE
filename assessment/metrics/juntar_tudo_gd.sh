#entradas="novo/novo"\("0,0"\)" smpso/smpso"\("2,2"\)" smpso/smpso"\("2,3"\)" smpso/smpso"\("2,4"\)" smpso/smpso"\("3,2"\)" smpso/smpso"\("3,3"\)" smpso/smpso"\("3,4"\)" smpso/smpso"\("4,2"\)" smpso/smpso"\("4,3"\)" smpso/smpso"\("4,4"\)""

entradas=$(cat results/titles.txt | tail -1)

comando=""
EFS="
"
#for problem in 1 2 3 4 5 6 7; do
problem=$1
	rm -f "results/all-$problem-gd.txt"	
	for objectives in 3 5 8 10 15 20; do
	#for objectives in 2 3 5 8 9; do
		for nomeSaida in $entradas; do
			comando=$comando"results/$nomeSaida""$problem-$objectives"_"fronts.txt "
		done
		java -Xmx1G -cp assessment/metrics/ gdp $comando >> "results/all-$problem-gd.txt"

		unset comando;
		echo $EFS >> "results/all-$problem-gd.txt"
		echo $EFS >> "results/all-$problem-gd.txt"
		echo $EFS >> "results/all-$problem-gd.txt"
		echo $EFS >> "results/all-$problem-gd.txt"
		echo $EFS >> "results/all-$problem-gd.txt"
		echo $EFS >> "results/all-$problem-gd.txt"
		echo $EFS >> "results/all-$problem-gd.txt"

		echo "$problem-$objectives - GD"
	done
#done
