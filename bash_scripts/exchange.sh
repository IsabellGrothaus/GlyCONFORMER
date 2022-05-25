# Bash script for calculating round-trip times of each replica

# Input 1: number of replica
# Input 2: timestep between replica exchange attempts in ps. 

declare -i n

nrep=$1
dt=$2
n=($nrep-1)

# demux trajectories with built of Giovanni Bussi to allow processing of long trajectories
demux.pl md.log

# Extract exchange probability between replica from md.log file
grep "Replica exchange statistics" -A11 md.log > REST_analysis.dat

# Calculation of runtrip time for each replica (to go from lowest to highest replica and backwards) + number of total exchanges

for i in `seq 0 ${n}`

do
    awk -v z=$i -v nrep=$nrep -v dt=$dt 'BEGIN{a=0; i=0; avg=0}{if($2 == z && a == 0){a=1; t=(NR-1)*dt; next}; if($(nrep+1) == z && a == 1){a=2; next}; if($2 == z && a == 2){a=1; r[i]=$1-t; i++; t=$1; next}}END{print "# Total number of round-trip times:",i; for(b=0;b<i;b++) avg+=r[b]/1000; print " # Average round-trip time (ns):", avg/i; print "# Values of round-trip times (ns):"; for(b=0;b<i;b++) print r[b]/1000}' replica_index.xvg  > rtt_rep${i}.dat
done

for i in `seq 0 ${n}`
do
 printf "$i "
 cat rtt_rep${i}.dat | grep Average
 cat rtt_rep${i}.dat | grep Average >> REST_analysis.dat
done

