#! /bin/bash

#export PATH=/usr/local2/rogerab/HTProcess:${PATH}
#export PATH=/usr/local2/FastQC:${PATH}

while getopts a:b:c:d:e:f:g:u:r: option
do
        case "${option}"
        in
                a) Library_name=${OPTARG};;
                b) Library_num=${OPTARG};;
                c) condition=${OPTARG};;
                d) Dir1=$OPTARG;;
                e) Dir2=${OPTARG};;
                f) pair_spacing=${OPTARG};;
                g) pair_sd=${OPTARG};;
                u) pair_type=$OPTARG;;
                r) DirS=$OPTARG;;
      esac
done
#export PATH=/home/rogerab/bin2:${PATH}
        PN1=`ls $Dir1 | wc -w`
        PN2=`ls $Dir2 | wc -w`
        SN=`ls $DirS | wc -w`
        if [ "$PN1">0 ];
        then 
		if [ "$PN2">0 ];
		then
        if [ "$PN1" == "$PN2" ];
        then
        pairing=paired
        else echo "There is a problem with pairing of read files in Directory 1 and Directory 2";
        fi
        fi
		fi
		
		if [[ "$pairing" == paired ]];
		then
		if [[ "$SN">0 ]];
		then 
		pairing="paired_and_unpaired"
		else
		pairing="paired-only"
		fi
		fi
		
		if [[ "$SN">0 ]];
		then 
        if [ "$PN1" == 0 ];
		then
		pairing="unpaired"
		fi
		fi

mkdir HTProcess_Reads
touch HTProcess_Reads/manifest_file.txt
echo "HTProcess1_Reads" > HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "Library_name=$Library_name" >> HTProcess_Reads/manifest_file.txt
echo "Library_num=$Library_num" >> HTProcess_Reads/manifest_file.txt
echo "condition=$condition" >> HTProcess_Reads/manifest_file.txt
        if [ "$PN1">0 ];
        then 
		echo "pair_spacing=$pair_spacing" >> HTProcess_Reads/manifest_file.txt
		echo "pair_sd=$pair_sd" >> HTProcess_Reads/manifest_file.txt
		echo "pair_type=$pair_type" >> HTProcess_Reads/manifest_file.txt
		fi
touch temp.txt
echo "pairing=$pairing" >> HTProcess_Reads/manifest_file.txt
echo "HTPROCESS1 `date`" >> HTProcess_Reads/HTProcess.log
echo "The pairing is $pairing" >> HTProcess_Reads/HTProcess.log
library_max=0


        if [ "$PN1">0 ];
        then 
		if [ "$PN2">0 ]; 
		then
m=

		echo Reads1 >> temp.txt
		for x in $Dir1/*
		do
		X=$(basename $x)
		XX=`echo $X | sed 's/.fastq//' | sed 's/.txt//'`
		XXN="$XX"'_fastqc'
		echo "!XXX $X" >> temp.txt
		m=`expr $m + 1`
		fastqc -t 4 $x -o ./
		wait
		E=`cat $XXN/fastqc_data.txt | grep -Po 'Illumina \d.\d' | sed 's/Illumina //'`
		M=`cat $XXN/fastqc_data.txt | grep -Po 'Sequence length\t[0-9].*' | sed 's/Sequence length\t//' | sed 's/[0-9]*-//'`
		echo "encoding_$X=$E" >> HTProcess_Reads/manifest_file.txt
		echo "max_len_$X=$M"  >> HTProcess_Reads/manifest_file.txt
		if [ $M > $library_max ] ;
		then
		library_max=$M
		fi
		done
echo "fastqc is finished testing $m files in the first paired read directory." >> HTProcess_Reads/HTProcess.log

r=
    echo Reads2 >> temp.txt
		for y in $Dir2/*
		do
		Y=$(basename $y)
		YY=`echo $Y | sed 's/.fastq//' | sed 's/.txt//'`
		YYN="$YY"'_fastqc'
                echo "!YYY $Y" >> temp.txt                     
		r=`expr $r + 1`
		fastqc -t 4 $y -o ./
		wait
		E=`cat $YYN/fastqc_data.txt | grep -Po 'Illumina \d.\d' | sed 's/Illumina //'`
		M=`cat $YYN/fastqc_data.txt | grep -Po 'Sequence length\t[0-9].*' | sed 's/Sequence length\t//' | sed 's/[0-9]*-//'`
		echo "encoding_$Y=$E" >> HTProcess_Reads/manifest_file.txt
		echo "max_len_$Y=$M"  >> HTProcess_Reads/manifest_file.txt
		if [ $M > $library_max ] ;
		then
		library_max=$M
		fi
		done
echo "fastqc is finished testing $r files in the second paired read directory."  >> HTProcess_Reads/HTProcess.log
fi
fi
	
		if [[ "$SN">0 ]];
		then 

p=

    	echo ReadsS >> temp.txt
		for z in $DirS/*
		do
		Z=$(basename $z)
		ZZ=`echo $Z | sed 's/.fastq//' | sed 's/.txt//'`
		ZZN="$ZZ"'_fastqc'
		
                echo "!ZZZ $Z" >> temp.txt                     
		p=`expr $p + 1`
		fastqc -t 4 $z -o ./
		wait
		E=`cat $ZZN/fastqc_data.txt | grep -Po 'Illumina \d.\d' | sed 's/Illumina //'`
		M=`cat $ZZN/fastqc_data.txt | grep -Po 'Sequence length\t[0-9].*' | sed 's/Sequence length\t//' | sed 's/[0-9]*-//'`
		echo "encoding_$Z=$E" >> HTProcess_Reads/manifest_file.txt
		echo "max_len_$Z=$M"  >> HTProcess_Reads/manifest_file.txt
		if [ $M > $library_max ] ;
		then
		library_max=$M
		fi
		done
fi

echo "library_max=$library_max"  >> HTProcess_Reads/manifest_file.txt
echo "fastqc is finished testing $p files in the directory for single reads."  >> HTProcess_Reads/HTProcess.log
 
    if [ "$PN1">0 ];
    then 
    if [ "$PN2">0 ];
	then
    if [ "$r" == "$m" ] ;
    then
	echo "Reads1 and Reads2 have the same number of files."  >> HTProcess_Reads/HTProcess.log
	echo "Testing for valid pairing."  >> HTProcess_Reads/HTProcess.log
	echo "Paired Reads" >> HTProcess_Reads/manifest_file.txt
	for ((d=0; d < $m; d += 1))
	do
		xarray=( $Dir1/* )
		yarray=( $Dir2/* )
		e=${xarray[$d]}
		f=${yarray[$d]}
		FastqPairedEndValidator.pl $e $f > "validate$d"
	V=`cat "validate$d" | grep -Po "properly ordered"`
	if [ "$V" == "properly ordered" ] ;
	then
	echo "!PPP $(basename $e),$(basename $f)"  >> HTProcess_Reads/manifest_file.txt
	echo "$(basename $e),$(basename $f) $V"  >> HTProcess_Reads/HTProcess.log
		cp $Dir1/* HTProcess_Reads/
		cp $Dir2/* HTProcess_Reads/
	else echo "There is a problem with the pairing of pair $e,$f." >> HTProcess_Reads/HTProcess.log
	fi
	done
	else echo "There is a problem with the pairing of the files in Reads1 and Reads2."
	fi
	fi
	fi
		
	if [[ "$SN">0 ]];
	then 
	cp $DirS/* HTProcess_Reads/
	fi
 	
 	grep -v 'encoding_\*=' HTProcess_Reads/manifest_file.txt | grep -v 'max_len_\*=' > tempA
 	cat tempA > HTProcess_Reads/manifest_file.txt
	rm tempA	
		
	echo "All Trim settings have been set to trim settings 1. Edit them on manifest_file.txt to customize trimming." >> HTProcess_Reads/HTProcess.log


cat temp.txt >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
echo "!!!TRIM SETTINGS!!!" >> HTProcess_Reads/manifest_file.txt
echo "......................................................" >> HTProcess_Reads/manifest_file.txt
grep '!PPP' HTProcess_Reads/manifest_file.txt | sed 's/!PPP/!PairTrim 1/' >> HTProcess_Reads/manifest_file.txt

if [[ "$SN">0 ]];
then 
grep '!ZZZ' HTProcess_Reads/manifest_file.txt | sed 's/!ZZZ/!SingleTrim 1/' >> HTProcess_Reads/manifest_file.txt
fi		
rm temp.txt

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
	mkdir individual_reports
	mv *report.html individual_reports/
	mv *fastqc.zip individual_reports/
	echo "HTPROCESS-FASTQC FINISHED `date`" >> HTProcess_Reads/HTProcess.log
