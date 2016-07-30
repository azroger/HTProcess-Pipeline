#! /bin/bash


#export PATH=/usr/local2/bin/:/usr/local2/rogerab/kmergenie-1.6476/:$PATH

while getopts a:b:c:d:e:f:g:h:k:l:m:n:o:p:q:r:s:t:u:v:x:y:z:i:j option
do
        case "${option}"
        in
        
                a) diploid=${OPTARG};;
                b) lowKmer=${OPTARG};;
                c) highKmer=${OPTARG};;
                d) interval=${OPTARG};;
                e) onePass=${OPTARG};;
        esac
done
#Program1
command1=''
if [[ -n "${diploid}" ]] ;
then
command1="$command1 --diploid"
fi
if [[ -n "${onePass}" ]] ;
then
command1="$command1 --one-pass"
fi
if [[ -n "${lowKmer}" ]] ;
then
command1="$command1 -l $lowKmer"
fi
if [[ -n "${highKmer}" ]] ;
then
command1="$command1 -k $highKmer"
fi
if [[ -n "${interval}" ]] ;
then
command1="$command1 -s $interval"
fi
if [[ -n "${sample}" ]] ;
then
command1="$command1 -e $sample"
fi

Trim=0
if [ -e HTProcess_Reads_T1 ];
then 
mv HTProcess_Reads_T1 HTProcess_Reads
Trim=1
fi
if [ -e HTProcess_Reads_T2 ];
then 
mv HTProcess_Reads_T2 HTProcess_Reads
Trim=1
fi
cp HTProcess_Reads/manifest_file.txt ./

Library_name=`grep 'Library_name' manifest_file.txt | sed 's/Library_name//' | sed 's/\=//'`
touch allreads.txt
if [[ "$Trim" == 1 ]];
then
count=`grep -c '!TRIMMED_Pr' manifest_file.txt`
array=(`grep '!TRIMMED_Pr' manifest_file.txt | sed 's/!TRIMMED_Pr//'`)
for ((i = 0; i < "$count"; i += 1));
do
read1=`echo "${array[$i]}" | sed 's/,.*//'`
read2=`echo "${array[$i]}" | sed 's/.*,//'`
echo "$read1" >> allreads.txt
echo "$read2" >> allreads.txt
done
singlename='TrmS_'"$Library_name"'.fastq'
if [[ -e HTProcess_Reads/$singlename ]];
then
echo "$singlename" >> allreads.txt
fi
fi

if [[ "$Trim" == 0 ]];
then
count=`grep -c '!PPP' manifest_file.txt`
array=(`grep '!PPP' manifest_file.txt | sed 's/!PPP//'`)
for ((i = 0; i < "$count"; i += 1));
do
read1=`echo "${array[$i]}" | sed 's/,.*//'`
read2=`echo "${array[$i]}" | sed 's/.*,//'`
echo "$read1" >> allreads.txt
echo "$read2" >> allreads.txt
done
countz=`grep -c '!ZZZ' manifest_file.txt`
arrayS=(`grep '!ZZZ' manifest_file.txt | sed 's/!ZZZ//'`)
for ((i = 0; i < "$countz"; i += 1));
do
readS=`echo "${arrayS[$i]}" | sed 's/,.*//'`
echo "$readS" >> allreads.txt
done
fi
cp allreads.txt HTProcess_Reads/
cd HTProcess_Reads
/usr/local2/bin/python2.7 /usr/local2/rogerab/kmergenie-1.6476/kmergenie allreads.txt -t 3 $command1  > log.txt
mkdir ../intermediate_files
mv histograms* ../intermediate_files/
cd ..
if [[ $Trim == 1 ]];
then
cat "intermediate_files/histograms_report.html" | sed "s/KmerGenie report/KmerGenie report for trimmed reads in the library named $Library_name/" > temp1.html
fi
if [[ $Trim == 0 ]];
then
cat "intermediate_files/histograms_report.html" | sed "s/KmerGenie report/KmerGenie report for untrimmed reads in the library named $Library_name/" > temp1.html
fi
cat temp1.html | sed 's/<\/body>//' | sed 's/<\/html>/<p>/' > temp2.html
echo "READS ANALYZED:" >> temp2.html
cat allreads.txt >> temp2.html
echo "<p>" >> temp2.html
echo "STANDARD OUT for Kmergenie:" >> temp2.html
cat HTProcess_Reads/log.txt >> temp2.html
echo '</body>' >> temp2.html
echo '</html>' >> temp2.html
cat temp2.html > Kmergenie_report.html
rm temp*
