#! /bin/bash
#export PATH=$PATH:/usr/local2/cufflinks-2.2.0.Linux_x86_64:/usr/local2/samtools-0.1.18/:/usr/local2/bowtie2-2.1.0/
#export PATH=$PATH:/usr/local2/tophat-2.0.11.Linux_x86_64:/usr/local2/samtools-0.1.19:/usr/local2/bowtie2-2.1.0
#Requires samtools and bowtie2
#Written by Roger Barthelson

while getopts a:b:c:d:e:f:g:h:k:l:m:n:o:p:q:r:s:t:u:v:x:y:z option
do
        case "${option}"
        in
        
                a) name_meth=${OPTARG};;
                b) single_cnt=${OPTARG};;
#                 c) library_type=${OPTARG};;
#                 d) compatible_hits=${OPTARG};;
#                 e) no_eff_length=${OPTARG};;
#                 f) no_length_corr=${OPTARG};;
#                 g) label=${OPTARG};;
#                 h) pre_mrna=${OPTARG};;
#                 k) option1=${OPTARG};;
#                 l) max_intron=${OPTARG};;
#                 m) assembly_mode=${assembly_mode};;
#                 r) reference=${OPTARG};;
#                 s) min_isoform=${OPTARG};;
#                 t) max_overhang=${OPTARG};;
        esac
done

set -x
command1=''
if [[ -n "${name_meth}" ]] ;
then
naming="${name_meth}"
fi


# fi
# if [[ -n "${library_type}" ]] ;
# then
# command1="$command1 --library-type $library_type"
# fi
# if [[ -n "${compatible_hits}" ]] ;
# then
# command1="$command1 --compatible-hits-norm"
# fi

mv HTProcess_BAM_* HTProcess_BAM
cp HTProcess_BAM/manifest_file.txt ./
mkdir HTProcess_BAMstats
exp_condition=`grep 'condition' manifest_file.txt | sed 's/condition//' | sed 's/\=//'`
mkdir intermediate_files
cp manifest_file.txt HTProcess_BAMstats/
cp HTProcess_BAM*/HTProcess.log HTProcess_BAMstats/
echo "......................................................" >>HTProcess_BAMstats/HTProcess.log
echo "......................................................" >> HTProcess_BAMstats/manifest_file.txt
echo "" >> HTProcess_BAMstats/manifest_file.txt
echo "" >> HTProcess_BAMstats/manifest_file.txt
echo "......................................................" >> HTProcess_BAMstats/manifest_file.txt
echo "HTPROCESS BAMstats is starting `date`" >> HTProcess_BAMstats/HTProcess.log

if [ $single_cnt = individual ];
then
count=`grep -e '!!OSBAM-file' -e '!!BAM_file' manifest_file.txt | grep -c -e 'TrmPr1' -e TrmSos -`
# Library_name=`grep 'Library_name' manifest_file.txt | sed 's/Library_name//' | sed 's/\=//'`
# condition=`grep 'condition' manifest_file.txt | sed 's/condition//' | sed 's/\=//'`
array=(`grep -e '!!OSBAM-file' -e '!!BAM_file' manifest_file.txt |  grep -e 'TrmPr1' -e TrmSos - | sed 's/!!BAM_file//' | sed 's/!!OSBAM-file//'`);
fi
if [ $single_cnt = combined ];
then
count=`grep -e '!!BAM_file' manifest_file.txt | grep -c -e 'TrmPr1' -e TrmS_ -`
# Library_name=`grep 'Library_name' manifest_file.txt | sed 's/Library_name//' | sed 's/\=//'`
# condition=`grep 'condition' manifest_file.txt | sed 's/condition//' | sed 's/\=//'`
array=(`grep -e '!!BAM_file' manifest_file.txt | grep -e 'TrmPr1' -e TrmS_ - | sed 's/!!BAM_file//'`);
fi

echo "........................................" >> HTProcess_BAMstats/manifest_file.txt
echo "!!!BAMstats Output files!!!" >> HTProcess_BAMstats/manifest_file.txt

cd HTProcess_BAM

for ((i = 0; i < "$count"; i += 1));
do
bamfile=`echo "${array[$i]}" | sed 's/,.*//'`
#singletest=`echo $bamfile | grep 'TrmS_'`
#if [[ -z "${singletest}" ]] ;

BAMstatsout=`echo $bamfile`
newname="$BAMstatsout"'.bamstats'
echo "Bamfile equals $bamfile" >> ../HTProcess_BAMstats/HTProcess.log
samtools index $bamfile
samtools idxstats $bamfile > $newname
echo "command: samtools index $bamfile" >> ../HTProcess_BAMstats/HTProcess.log
echo "command: samtools idxstats $bamfile" >> ../HTProcess_BAMstats/HTProcess.log
wait
if [ $i = 0 ];
then
cut -f 1,3 $newname > combo.tsv
if [ "$naming" = 'bam' ] ;
then
printf 'id'"\t""$BAMstatsout" > header.tsv
fi
if [ "$naming" = 'condition' ] ;
then
printf 'condition'"\t""$exp_condition" > header.tsv
fi
cp $newname 'group.tsv'
fi
if [[ $i -gt 0 ]];
then
if [ "$naming" = 'bam' ] ;
then
printf "\t""$BAMstatsout" >> header.tsv
fi
if [ "$naming" = 'condition' ] ;
then
printf "\t""$exp_condition" >> header.tsv
fi

cut -f 1,3 $newname > single.tsv

cp combo.tsv group.tsv

join -1 1 -2 1 group.tsv single.tsv > combo.tsv
fi

done
printf "\n" >> header.tsv

# cut -f 1 combo.tsv >> counts.tsv
# for ((j = 0; j < "$count"; j += 1)); 
# do
# let y=($j*4)+2
# cut -f $y combo.tsv >> counts.tsv
# done
cat combo.tsv | sed 's/ /\t/g' > counts.tsv
cat header.tsv counts.tsv >> ../BAM_counts.tsv
cp ../BAM_counts.tsv ../HTProcess_BAMstats/BAM_counts.tsv


echo "......................................................" >> ../HTProcess_BAMstats/manifest_file.txt
echo "BAM_counts.tsv" >> ../HTProcess_BAMstats/manifest_file.txt
echo "" >> ../HTProcess_BAMstats/manifest_file.txt
echo "......................................................" >> ../HTProcess_BAMstats/manifest_file.txt
#cat temp >> ../HTProcess_BAMstats/manifest_file.txt
#condition=`grep 'condition' manifest_file.txt | sed 's/condition//' | sed 's/\=//'`.
cp *.bai ../intermediate_files/
cp *.bamstats ../intermediate_files/
cd ..

rm manifest_file.txt
echo "Finished Running BAMstats. `date`" >> HTProcess_BAMstats/HTProcess.log

mv HTProcess_BAMstats HTProcess_BAMstats'_'"$exp_condition"
rm -r HTProcess_BAM
set +x

