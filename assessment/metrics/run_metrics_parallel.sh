problems="dtlz1 dtlz2 dtlz3 dtlz4 dtlz5 dtlz6 dtlz7 wfg1 wfg2 wfg3 wfg4 wfg5 wfg6 wfg7 wfg8 wfg9"
#problems="dtlz1"
objectives="3 5 8 10 15 20"
#objectives="3"
# metrics="gdp igdp r2 hv"
metrics="hv"

entradas=$(cat results/titles.txt | tail -1)
comando=""
EFS="
"

for metric in $metrics; do
	for problem in $problems; do
		for objective in $objectives; do
			for nomeSaida in $entradas; do
				comando=$comando"results/$nomeSaida""$problem-$objective"_"fronts.txt "
			done
			front=`echo $problem | tr 'a-z' 'A-Z'`
			
			if [ $front = "DTLZ3" ] || [ $front = "DTLZ4" ]; then
				front="DTLZ2"
			fi
			if [ $front = "DTLZ6" ]; then
				front="DTLZ5"
			fi
			front=$front'_'$objective

			
	# 		java -Xmx1G -cp assessment/metrics/ r2 $comando >> "results_bkp/all-$problem-r2.$objective.txt"
	# 		echo files: $comando
	# 		echo front: "assessment/pareto/$front assessment/pareto/REF_$objective"
	# 		echo saida: "results_bkp/all-$problem-r2.$objective.txt"

			qsub assessment/metrics/qsubMetrics.sh "$metric" "$comando" "assessment/metrics/pareto/$front assessment/metrics/pareto/REF_$objective" "results/all-$problem-$metric.$objective.txt"

			unset comando;
		done
	done
done
