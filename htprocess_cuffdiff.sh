#! /bin/bash

# Created by Roger Barthelson to run in the iPlant DE
# A step of the HTProcess pipeline
# Runs cuffdiff on all BAM files in input HTProcess_BAM directories
# Uses custom or std annotation gtf or gff files 



export PATH=$PATH:/usr/local2/cufflinks-2.2.0.Linux_x86_64:/usr/local2/samtools-0.1.18/:/usr/local2/bowtie2-2.1.0/
#export PATH=$PATH:/usr/local2/tophat-2.0.11.Linux_x86_64:/usr/local2/samtools-0.1.19:/usr/local2/bowtie2-2.1.0

while getopts a:b:c:d:e:f:g:h:k:l:m:n:o:p:q:r:s:t:u:v:x:y:z option
do
        case "${option}"
        in
        
                a) user_config=${OPTARG};;
                b) single_cnt=${OPTARG};;
#                c) library_type=${OPTARG};;
#                 d) compatible_hits=${OPTARG};;
#                 e) no_eff_length=${OPTARG};;
#                 f) no_length_corr=${OPTARG};;
                g) trnscrpt_asmbly=${OPTARG};;
#                 h) pre_mrna=${OPTARG};;
#                 k) option1=${OPTARG};;
#                 l) max_intron=${OPTARG};;
#                 m) assembly_mode=${assembly_mode};;
#                 r) reference=${OPTARG};;
#                 s) min_isoform=${OPTARG};;
#                 t) max_overhang=${OPTARG};;
        esac
done

Conditon=''
count=''
array=''
trnscrpt_asmbly="${trnscrpt_asmbly}"

#dirCnt=`ls | grep 'HTProcess_BAM' | wc -l`

for x in HTProcess_BAM*/manifest_file.txt
do

Condition=`grep condition $x | sed 's/condition=//'` 
if [ $single_cnt = individual ];
then
count=`grep -e '!!OSBAM-file' -e '!!BAM_file' $x | grep -c -e 'TrmPr1' -e TrmSos -`
array=(`grep -e '!!OSBAM-file' -e '!!BAM_file' $x |  grep -e 'TrmPr1' -e TrmSos - | sed 's/!!BAM_file//' | sed 's/!!OSBAM-file//'`);
fi
if [ $single_cnt = combined ];
then
count=`grep -e '!!BAM_file' $x | grep -c -e 'TrmPr1' -e TrmS_ -`
array=(`grep -e '!!BAM_file' $x | grep -e 'TrmPr1' -e TrmS_ - | sed 's/!!BAM_file//'`);
fi
#done

# for ((j = 0; j < "$dirCnt"; j += 1));
# do
# condCnt=`"${count[$j]}"`
# for ((k = 0; k < "$condCnt"; k += 1));
# Kondition+=`

# tot=0
# for m in ${count[@]}; do
#   let tot+=$m
# done
# echo "Total: $tot"
for ((i = 0; i < "$count"; i += 1));
do
# let j=$i+1
bamfile=`echo "${array[$i]}" | sed 's/,.*//'`
Cond=`echo "$Condition" | sed 's/,.*//'`
echo "$Cond	$bamfile" >> new_config.txt
done
wait
done

######################################
#PART  TWO
######################################

if [[ -r $user_config ]];
then cat $user_config > cuffdiff_config.txt
fi
if [[ -z $user_config ]];
then cat new_config.txt > cuffdiff_config.txt
cat cuffdiff_config.txt > config.txt
fi
echo "HTPROCESS CuffDiff is starting `date`" > HTProcess.log
cat config.txt >> HTProcess.log

switch=`cut -f 1 config.txt`
printf "$switch" >> LABELS

#run cuffdiff


names=(`cut -f 1 config.txt | uniq`)
nameList=(`cut -f 1 config.txt`)
bamList=(`cut -f 2 config.txt`)

ttl=`echo ${names[@]} | wc -w`
tot=`echo ${nameList[@]} | wc -w`
echo "TTL: $ttl" >> HTProcess.log
echo "Total: $tot" >> HTProcess.log
Labels=''

for ((m = 0; m < "$ttl"; m += 1)); do
y="${names[$m]}"
Labels+=','"$y"
printf ' ' >> inputfiles

for ((j = 0; j < "$tot"; j += 1)); do
if [ ${nameList[$j]} == $y ];
then
printf ','"${bamList[$j]}" >> inputfiles
fi
inPuts=`less inputfiles | sed 's/ ,/ /g'`
done
done
#printf "$Labels" > LABELS
labelz=`echo "$Labels" | sed 's/^,//'`
cp HTProcess_BAM*/*.bam ./
echo "....................................." >> HTProcess.log
echo " Files present in directory. " >> HTProcess.log
ls >> HTProcess.log
/usr/local2/cufflinks-2.2.0.Linux_x86_64/cuffdiff -p 4 -o CUFFDIFFOUTPUT -L $labelz $trnscrpt_asmbly $inPuts
echo "cuffdiff -p 4 -o CUFFDIFFOUTPUT -L $labelz ${trnscrpt_asmbly} $inPuts" >> HTProcess.log
rm *.bam

echo "Sorting output data and graphing graphs with CummeRbund" >> HTProcess.log
perl /usr/local2/rogerab/HTProcess/perlcleanup.pl

rm CleanupOut
rm cuffdiff_config.txt
rm new_config.txt


echo "HTPROCESS CuffDiff is finished `date`" >> HTProcess.log
