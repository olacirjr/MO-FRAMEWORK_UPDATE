#/bin/bash
#########################################
## Parametros que le pasamos al script ##
#########################################
#$ -S /bin/bash
#######################################
# Usar el directorio de trabajo actual
#######################################
#$ -cwd
# Tiempo de trabajo
#$ -l h_rt=2400:00:00
# juntar la salida estandar y de error en un solo fichero
#$ -j y
###########################
# usar colas indicadas
###########################
#$ -q pegasus.q,gemini.q,loki.q,libra.q,others.q

##$ -t 1-30:1
##$ -o /dev/null

# scrt=/home/olacir/metrics/tmp.$JOB_ID
# mkdir -p $scrt

# cp --parents $2 $3 $scrt
# cp --parents assessment/metrics/$1.class assessment/metrics/Front.class assessment/metrics/hv/wfg  $scrt

# cd $scrt

echo "Host --> $HOSTNAME"

echo Init Time: `date`
init=`date +'%s'`

echo "java -Xmx1G -cp assessment/metrics/ $1 \"$(echo $2 | sed -e "s/ /\" \"/g")\" > $4"
eval "java -Xmx1G -cp assessment/metrics/ $1 \"$(echo $2 | sed -e "s/ /\" \"/g")\" > $4"

final=`date +'%s'`
echo End  Time: `date` --- run in $(( (($final-$init)/3600) )):$(( (($final-$init)/60)%60 )):$(( (($final-$init))%60 )) ---

# cp $scrt/results/* $SGE_O_WORKDIR/results/
# rm -rf $scrt #/home/olacir/metrics/

mkdir -p "$SGE_O_WORKDIR/logs/metrics/"
mv "$SGE_O_WORKDIR/qsubMetrics.sh.o$JOB_ID" "$SGE_O_WORKDIR/logs/metrics/"

echo "Pronto!" --`date`--

