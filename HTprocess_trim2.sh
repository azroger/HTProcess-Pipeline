#! /bin/bash

# Created by Roger Barthelson to run in the iPlant DE
# Second main step of HTProcess pipeline
# Takes output (HTProcess_Reads directory) from HTProcess_fastqc
# or from HTProcess-prepare_directories_and_run_fastqc (workflow)
# Runs trimmomatic on all fastq files input
# Runs FastQC on trimmed reads
# Creates combined summary file in html that includes embedded images

#export PATH=/usr/local2/rogerab/HTProcess:${PATH}
#export PATH=/usr/local2/FastQC:${PATH}
#export PATH=/usr/local2/Trimmomatic-0.32:${PATH}

while getopts a:b:c:d:e:f:g:h:k:l:m:n:o:p:q:r:s:t:u:v:x:y:z:i:j option
do
        case "${option}"
        in
        
                a) headcrop1=${OPTARG};;
                b) crop1=${OPTARG};;
                c) ICmismatch1=${OPTARG};;
                d) ICpalindrome1=${OPTARG};;
                e) ICmatch1=${OPTARG};;
                f) leading1=${OPTARG};;
                g) trailing1=${OPTARG};;
                h) SWwindowSize1=${OPTARG};;
                k) SWquality1=${OPTARG};;
                l) MItargetLength1=${OPTARG};;
                m) MIstrictness1=${OPTARG};;
                n) headcrop2=${OPTARG};;
                o) crop2=${OPTARG};;
                p) ICmismatch2=${OPTARG};;
                q) ICpalindrome2=${OPTARG};;
                r) ICmatch2=${OPTARG};;
                s) leading2=${OPTARG};;
                t) trailing2=${OPTARG};;
                u) SWwindowSize2=${OPTARG};;
                v) SWquality2=${OPTARG};;
                x) minLength=${OPTARG};;                
                y) ADAPTERS=${OPTARG};;
                z) MItargetLength2=${OPTARG};;
                i) MIstrictness2=${OPTARG};;
        esac
done
#Program1
command1=''
if [[ -n "${headcrop1}" ]] ;
then
command1="$command1 HEADCROP:$headcrop1"
fi
if [[ -n "${ICmismatch1}" ]] ;
then
command1="$command1 ILLUMINACLIP:$ADAPTERS:$ICmismatch1:$ICpalindrome1:$ICmatch1"
fi
if [[ -n "${leading1}" ]] ;
then
command1="$command1 LEADING:$leading1"
fi
if [[ -n "${trailing1}" ]] ;
then
command1="$command1 TRAILING:$trailing1"
fi
if [[ -n "${SWwindowSize1}" ]] ;
then
command1="$command1 SLIDINGWINDOW:$SWwindowSize1:$SWquality1"
fi
if [[ -n "${MItargetLength1}" ]] ;
then
command1="$command1 MAXINFO:$MItargetLength1:$MIstrictness1"
fi
if [[ -n "${crop1}" ]] ;
then
command1="$command1 CROP:$crop1"
fi
if [[ -n "${minLength}" ]] ;
then
command1="$command1 MINLEN:$minLength"
fi
#Program2
command2=''
if [[ -n "${headcrop2}" ]] ;
then
command2="$command2 HEADCROP:$headcrop2"
fi
if [[ -n "${ICmismatch2}" ]] ;
then
command2="$command2 ILLUMINACLIP:$ADAPTERS:$ICmismatch2:$ICpalindrome2:$ICmatch2"
fi
if [[ -n "${leading2}" ]] ;
then
command2="$command2 LEADING:$leading2"
fi
if [[ -n "${trailing2}" ]] ;
then
command2="$command2 TRAILING:$trailing2"
fi
if [[ -n "${SWwindowSize2}" ]] ;
then
command2="$command2 SLIDINGWINDOW:$SWwindowSize2:$SWquality2"
fi
if [[ -n "${MItargetLength2}" ]] ;
then
command2="$command2 MAXINFO:$MItargetLength2:$MIstrictness2"
fi
if [[ -n "${crop2}" ]] ;
then
command2="$command2 CROP:$crop2"
fi
if [[ -n "${minLength}" ]] ;
then
command2="$command2 MINLEN:$minLength"
fi
#ILLUMINACLIP:barcodes.fa:2:40:15 \
#LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25

if [ -e HTProcess_ReadsX ];
then 
mv HTProcess_ReadsX HTProcess_Reads
fi
cp HTProcess_Reads/manifest_file.txt ./
mv HTProcess_Reads HTProcess_Reads1
mkdir HTProcess_Reads
cp manifest_file.txt HTProcess_Reads/
cp HTProcess_Reads1/HTProcess.log HTProcess_Reads/
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "!!!TRIMMED READS!!!" >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "HTPROCESS Trimmomatic is starting `date`" >> HTProcess_Reads/HTProcess.log

pairing=`grep 'pairing' manifest_file.txt | sed 's/pairing=//'`
fred=`grep -m 1 'encoding' manifest_file.txt | sed s/encoding.*\=//`
if [[ $fred == 1.5 ]];
then
Fred='-phred64' 
fi
if [[ $fred == 1.9 ]];
then
Fred='-phred33'
fi
mkdir singletemp
touch "orphans"
count=`grep -c '!PairTrim' manifest_file.txt`
Library_name=`grep 'Library_name' manifest_file.txt | sed 's/Library_name//' | sed 's/\=//'`
let twocount="2 * $count"
array=(`grep '!PairTrim' manifest_file.txt | sed 's/!PairTrim//'`);

for ((i = 0; i < "$twocount"; i += 1));
do
if [[ "${array[$i]}" == [1,2] ]];
then
echo "trim program is ${array[$i]}" >> HTProcess_Reads/HTProcess.log
program="${array[$i]}"
j=`expr $i + 1`
read1=`echo "${array[$j]}" | sed 's/,.*//'`
read2=`echo "${array[$j]}" | sed 's/.*,//'`
echo "read1 equals $read1" >> HTProcess_Reads/HTProcess.log
echo "read2 equals $read2" >> HTProcess_Reads/HTProcess.log
if [[ $program == 1 ]];
then
java -jar trimmomatic-0.32.jar PE -threads 4 "$Fred" -trimlog trimlogFilep1.txt HTProcess_Reads1/"$read1" HTProcess_Reads1/"$read2" HTProcess_Reads/"TrmPr1_""$read1" singletemp/"TrmS1_""$read1" HTProcess_Reads/"TrmPr2_""$read2" singletemp/"TrmS2_""$read2" $command1
echo "!TRIMMED_Pr TrmPr1_$read1,TrmPr2_$read2" >> HTProcess_Reads/manifest_file.txt
cat singletemp/"TrmS2_$read2" >> singletemp/"TrmS1_$read1"
mv singletemp/"TrmS1_$read1" singletemp/"TrmSos_$read1"
rm singletemp/"TrmS2_""$read2"
echo "!TRIMMED_OS TrmSos_$read1" >> orphans
fi
if [[ $program == 2 ]];
then
java -jar trimmomatic-0.32.jar PE -threads 4 "$Fred" -trimlog trimlogFilep1.txt HTProcess_Reads1/"$read1" HTProcess_Reads1/"$read2" HTProcess_Reads/"TrmPr1_""$read1" singletemp/"TrmS1_""$read1" HTProcess_Reads/"TrmPr2_""$read2" singletemp/"TrmS2_""$read2" $command2
echo "!TRIMMED_Pr TrmPr1_$read1,TrmPr2_$read2" >> HTProcess_Reads/manifest_file.txt
cat singletemp/"TrmS2_$read2" >> singletemp/"TrmS1_$read1"
mv singletemp/"TrmS1_$read1" singletemp/"TrmSos_$read1"
rm singletemp/"TrmS2_$read2"
echo "!TRIMMED_OS TrmSos_$read1" >> orphans
fi
fi
done
echo "Finished trimming paired reads files." >> HTProcess_Reads/HTProcess.log
#SINGLES
cnt=`grep -c '!SingleTrim' manifest_file.txt`
let twocnt="2 * $cnt"
arrayz=(`grep '!SingleTrim' manifest_file.txt | sed 's/!SingleTrim//'`);
for ((i = 0; i < "$twocnt"; i += 1));
do
if [[ "${arrayz[$i]}" == [1,2] ]];
then
echo "trim program is ${arrayz[$i]}" >> HTProcess_Reads/HTProcess.log
program="${arrayz[$i]}"
j=`expr $i + 1`
readS=`echo "${arrayz[$j]}" | sed 's/,.*//'`
echo "readS equals $readS" >> HTProcess_Reads/HTProcess.log
if [[ $program == 1 ]];
then
java -jar trimmomatic-0.32.jar SE -threads 4 "$Fred" -trimlog trimlogFileS.txt HTProcess_Reads1/"$readS" singletemp/"TrmSos_""$readS" $command1
echo "!TRIMMED_OS TrmSos_""$readS" >> orphans
fi
if [[ $program == 2 ]];
then
java -jar trimmomatic-0.32.jar SE -threads 4 "$Fred" -trimlog trimlogFileS.txt HTProcess_Reads1/"$readS" singletemp/"TrmSos_""$readS" $command2
echo "!TRIMMED_OS TrmSos_""$readS" >> orphans
fi
fi
done
newname='TrmS_'"$Library_name"'.fastq'
cat singletemp/TrmS* > HTProcess_Reads/"$newname"
echo "!TRIMMED_S $newname" >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "!!!TRIMMED ORPHAN AND INDIVIDUAL SINGLES!!!" >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "Not used for normal analysis with a completely uniform library" >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
cat "orphans" >> HTProcess_Reads/manifest_file.txt
rm -r HTProcess_Reads1
rm orphans
mv singletemp/* HTProcess_Reads/
rm -r singletemp


#START FASTQC SECTION
echo "Starting to run fastqc on trimmed reads." >> HTProcess_Reads/HTProcess.log
counter=`grep -c '!TRIMMED_Pr' HTProcess_Reads/manifest_file.txt`
array=(`grep '!TRIMMED_Pr' HTProcess_Reads/manifest_file.txt | sed 's/!TRIMMED_Pr//'`);
for ((i = 0; i < "$counter"; i += 1));
do
#j=`expr $i + 1`
read1=`echo "${array[$i]}" | sed 's/,.*//'`
read2=`echo "${array[$i]}" | sed 's/.*,//'`

		fastqc -t 4 HTProcess_Reads/$read1 -o ./
		wait

		fastqc -t 4 HTProcess_Reads/$read2 -o ./
		wait
done
		fastqc -t 4 HTProcess_Reads/$newname -o ./
#Orphan single section
counter=`grep -c '!TRIMMED_OS' HTProcess_Reads/manifest_file.txt`
array=(`grep '!TRIMMED_OS' HTProcess_Reads/manifest_file.txt | sed 's/!TRIMMED_OS//'`);
#array=`echo "${array[$i]}" | sed 's/,//'`
for ((i = 0; i < "$counter"; i += 1));
do		
readOS=`echo "${array[$i]}"`
		fastqc -t 4 HTProcess_Reads/$readOS -o ./
		wait
done		
		echo "Finished running fastqc on trimmed reads." >> HTProcess_Reads/HTProcess.log
#mv HTProcess_Reads/*fastqc ./
#mv HTProcess_Reads/*zip ./		
		
#Start HTML Summary Section

echo "Starting Summary of fastqc outputs for trimmed reads." >> HTProcess_Reads/HTProcess.log
echo "Starting creation of summary file for FASTQC reports" >> HTProcess_Reads/HTProcess.log
echo "First Phase of HTPROCESS1 FINISHED `date`" >> HTProcess_Reads/HTProcess.log

	echo '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Strict//EN">' > fastqc_summary.html
	echo '<html' >> fastqc_summary.html
	echo '<head><title>Summary of Fastqc Reports</title>' >> fastqc_summary.html
	echo '<style type="text/css">' >> fastqc_summary.html
	echo '	body { font-family: sans-serif; color: #0098aa; background-color: #FFF; font-size: 100%; border: 0; margin: 0; padding: 0; }' >> fastqc_summary.html
	echo '	h1 { font-family: sans-serif; color: #0098aa; background-color: #FFF; font-size: 300%; font-weight: bold; border: 0; margin: 0; padding: 0; }' >> fastqc_summary.html
	echo '	h2 { font-family: sans-serif; color: #0098aa; background-color: #FFF; font-size: 200%; font-weight: bold; border: 0; margin: 0; padding: 0; }' >> fastqc_summary.html	
	echo '	h3 { font-family: sans-serif; color: #0098aa; background-color: #FFF; font-size: 40%; font-weight: bold; border: 0; margin: 0; padding: 0; }' >> fastqc_summary.html
	echo '	.TFtable tr:nth-child(even){ background: #D2DADC; }'	>> fastqc_summary.html	
	echo '	</style>' >> fastqc_summary.html
	echo '	</head>' >> fastqc_summary.html
	echo '	<h1> Summary of Fastqc Reports' >> fastqc_summary.html
	echo '	<br/>' >> fastqc_summary.html
	echo '	<br/>' >> fastqc_summary.html
	echo '	<br/> </h1>' >> fastqc_summary.html
	echo '	<body> ' >> fastqc_summary.html
	echo '	<table border="1" cellpadding="10" bgcolor="white" class="TFtable">' >> fastqc_summary.html
	echo '	<tr>' >> fastqc_summary.html
	echo '	<td><b>Fastq File Name</b></td>' >> fastqc_summary.html
	echo '    <td><b>Basic Statistics</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per base sequence quality</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per sequence quality scores</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per base sequence content</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per base GC content</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per sequence GC content</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per base N content</b></td>' >> fastqc_summary.html
	echo '    <td><b>Sequence Length Distribution</b></td>' >> fastqc_summary.html
	echo '    <td><b>Sequence Duplication Levels</b></td>' >> fastqc_summary.html
	echo '    <td><b>Overrepresented sequences</b></td>' >> fastqc_summary.html
	echo '    <td><b>Kmer Content</b></td>' >> fastqc_summary.html
	echo ' </tr>' >> fastqc_summary.html

	WW=0
	for w in *_fastqc
	do
	WW=`expr $WW + 1`
	echo ' <tr>' >> fastqc_summary.html
	echo '	<div>' >> fastqc_summary.html
	echo '	<ul>' >> fastqc_summary.html
	cd $w
	ww=`echo $w | sed 's/_fastqc/_report.html'/`
	embed_images.pl fastqc_report.html ../"$ww"
	csplit -f sum ../"$ww" /\<li\>/ {10} /class=\"main\"\>/
	test=`grep 'Kmer graph' ../"$ww"`	
	if [[ -n $test ]];
	then 	
	csplit ../"$ww" /'<div class="module">'/ {10} /'alt="Kmer graph"'/
	else
	csplit ../"$ww" /'<div class="module">'/ {10} /'<p>No overrepresented Kmers</p>'/ 
	fi
	seq_file=`grep title\> xx00 | sed 's/<head><title>//'  | sed 's/FastQC\ Report<\/title>//'`
	echo "	<td><b>$seq_file</b></td>" >> ../fastqc_summary.html
	sumRe=( sum01 sum02 sum03 sum04 sum05 sum06 sum07 sum08 sum09 sum10 sum11 )
	for ((RR=0; RR < 11; RR += 1))
	do
	cellcontent=`cat "${sumRe[$RR]}" | sed s/#M/\#\$WW-M/`
	echo "	<td><b>$cellcontent</b></td>" >> ../fastqc_summary.html
	done		
	cd ..
	echo ' </tr>' >> fastqc_summary.html
	done
	echo '</table>' >> fastqc_summary.html
	echo '	<br/>' >> fastqc_summary.html
	echo '	<br/>' >> fastqc_summary.html
	echo '	<br/>' >> fastqc_summary.html
	echo '	<br/>' >> fastqc_summary.html
	echo '<table border="1" cellpadding="10" bgcolor="white" class="TFtable">' >> fastqc_summary.html
	echo '	<tr>' >> fastqc_summary.html
	echo '    <td><b>Basic Statistics</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per base sequence quality</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per sequence quality scores</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per base sequence content</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per base GC content</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per sequence GC content</b></td>' >> fastqc_summary.html
	echo '    <td><b>Per base N content</b></td>' >> fastqc_summary.html
	echo '    <td><b>Sequence Length Distribution</b></td>' >> fastqc_summary.html
	echo '    <td><b>Sequence Duplication Levels</b></td>' >> fastqc_summary.html
	echo '    <td><b>Overrepresented sequences</b></td>' >> fastqc_summary.html
	echo '    <td><b>Kmer Content</b></td>' >> fastqc_summary.html
	echo ' </tr>' >> fastqc_summary.html
	WW=0
	for w in *_fastqc
	do
	WW=`expr $WW + 1`
	echo ' <tr>' >> fastqc_summary.html	
	cd $w
	
	seq_file=`grep title\> xx00 | sed 's/<head><title>//'  | sed 's/FastQC\ Report<\/title>//'`
	
	xxRe=( xx01 xx02 xx03 xx04 xx05 xx06 xx07 xx08 xx09 xx10 xx11 )	
	for ((SS=0; SS < 10; SS += 1))
	do
	graphstuff=`cat "${xxRe[$SS]}" | sed 's/"indented" src=/"indented" height="320" width="400" src=/g'`
	graphcontent="<h3 id=$WW-M$SS > $seq_file $graphstuff </h3>"
	echo "	<td><b>$graphcontent</b></td>" >> ../fastqc_summary.html
	done
	test=`grep 'Kmer graph' ../"$ww"`	
	if [[ -n $test ]];
	then
	graphstuff11=`cat "xx11" | sed 's/"indented" src=/"indented" height="320" width="400" src=/g'`
	graphcontent11="<h3 id=$WW-M10 > $seq_file $graphstuff11 "'"</h2></p></div></h3>'
    else
	graphstuff11=`cat "xx11"`
	graphcontent11="<h3 id=$WW-M10 > $seq_file $graphstuff11 "'"alt="[OK]"> Kmer Content</h2><h3><p>No overrepresented Kmers</p></div></h3>'
	fi
	echo "	<td><b>$graphcontent11</b></td>" >> ../fastqc_summary.html
	cd ..
	echo ' </tr>' >> fastqc_summary.html
	done
	echo '</table>' >> fastqc_summary.html					
	echo '</body>' >> fastqc_summary.html
	echo '</html' >> fastqc_summary.html
	
	echo "The summary file for all the FASTQC reports has been created." >> HTProcess_Reads/HTProcess.log
	rm -r *_fastqc
	if [[ -n "$Dir1" ]];
    then 
	rm validate*
	fi
	mkdir intermediate_files
	mv *report.html intermediate_files/
	mv *fastqc.zip intermediate_files/
	rm manifest_file.txt
	mv trimlog* intermediate_files/


echo "HTProcess Trimmomatic is finished `date`" >> HTProcess_Reads/HTProcess.log
mv HTProcess_Reads HTProcess_Reads_T1
