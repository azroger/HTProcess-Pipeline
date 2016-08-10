#! /bin/bash


# Created by Roger Barthelson to run in the iPlant DE
# A step of the HTProcess pipeline
# Runs cufflinks on all BAM files in input HTProcess_BAM directories 


#export PATH=$PATH:/usr/local2/cufflinks-2.2.0.Linux_x86_64:/usr/local2/samtools-0.1.18/:/usr/local2/bowtie2-2.1.0/
#export PATH=$PATH:/usr/local2/tophat-2.0.11.Linux_x86_64:/usr/local2/samtools-0.1.19:/usr/local2/bowtie2-2.1.0

while getopts a:b:c:d:e:f:g:h:k:l:m:n:o:p:q:r:s:t:u:v:x:y:z option
do
        case "${option}" in
        
                a) assemMode=${OPTARG};;
                b) multi_read=${OPTARG};;
                c) library_type=${OPTARG};;
                d) compatible_hits=${OPTARG};;
                e) min_isoform=${OPTARG};;
                f) pre_mrna=${OPTARG};;
                g) no_eff_length=${OPTARG};;
                h) no_length_corr=${OPTARG};;
                m) option1=${OPTARG};;
                n) max_intron=${OPTARG};;
                r) reference=${OPTARG};;
                s) single_assembly=${OPTARG};;
                v) max_overhang=${OPTARG};;
        esac
done

command1=''
#if [[ -n "${frag_bias}" ]] ;
#then
#command1="$command1 -b $frag_bias"
#fi
if [[ -n "${multi_read}" ]] ;
then
command1="$command1 -u"
fi
if [[ -n "${library_type}" ]] ;
then
command1="$command1 --library-type $library_type"
fi
if [[ -n "${compatible_hits}" ]] ;
then
command1="$command1 --compatible-hits-norm"
fi
if [[ -n "${no_eff_length}" ]] ;
then
command1="$command1 --no-effective-length-correction"
fi
if [[ -n "${no_length_corr}" ]] ;
then
command1="$command1 --no-length-correction"
fi
#if [[ -n "${label}" ]] ;
#then
#command1="$command1 -L $label"
#fi
if [[ -n "${pre_mrna}" ]] ;
then
command1="$command1 -j $pre_mrna"
fi
if [[ -n "${min_isoform}" ]] ;
then
command1="$command1 -F $min_isoform"
fi
if [[ -n "${option1}" ]] ;
then
command1="$command1 $option1"
fi
if [[ -n "${max_intron}" ]] ;
then
command1="$command1 -I $max_intron"
fi
if [[ -n "${max_overhang}" ]] ;
then
command1="$command1 --overhang-tolerance $max_overhang"
fi
gene_fpkm=''
Iso_fpkm=''
mkdir HTProcess_Cufflinks
cp HTProcess_BAM*/manifest_file.txt ./
mkdir HTProcess_Cufflinks/intermediate_files
cp manifest_file.txt HTProcess_Cufflinks/
cp HTProcess_BAM*/HTProcess.log HTProcess_Cufflinks/
echo "......................................................" >>HTProcess_Cufflinks/HTProcess.log
echo "......................................................" >> HTProcess_Cufflinks/manifest_file.txt
echo "" >> HTProcess_Cufflinks/manifest_file.txt
echo "" >> HTProcess_Cufflinks/manifest_file.txt
echo "......................................................" >> HTProcess_Cufflinks/manifest_file.txt
echo "HTPROCESS Cufflinks is starting `date`" >> HTProcess_Cufflinks/HTProcess.log

if [ -e "${reference}" ] ;
then
echo "......................................................" >> HTProcess_Cufflinks/manifest_file.txt
echo "!!!Reference-GTF-is-Defined!!!" >> HTProcess_Cufflinks/manifest_file.txt
echo "!GGG $reference" >> HTProcess_Cufflinks/manifest_file.txt
echo "......................................................" >> HTProcess_Cufflinks/manifest_file.txt
fi
if [ "${assemMode}" = 'refOnly'  ] ;
then
command1="$command1 --GTF ${reference}"
echo "Cufflinks will use Reference Annotation ONLY!" >> HTProcess_Cufflinks/HTProcess.log
fi
if [ "${assemMode}" = 'guided' ] ;
then
command1="$command1 -g ${reference}"
fi
if [ "${assemMode}" = 'noRef' ] ;
then
command1="$command1"
fi

#echo "command1 says this  $command1" > mycommandr
# echo "assembly_mode is  ${assemMode}" >> mycommandr
#echo "the GTF file is  ${reference}" >> mycommandr

Assembl='!!BAM_file'
if [ "${single_assembly}" == 'yes' ]; then Assembl='!!!MergedBAM' ; fi

count=`grep -c $Assembl manifest_file.txt`
# Library_name=`grep 'Library_name' manifest_file.txt | sed 's/Library_name//' | sed 's/\=//'`
# condition=`grep 'condition' manifest_file.txt | sed 's/condition//' | sed 's/\=//'`
array=(`grep $Assembl manifest_file.txt | sed s/"$Assembl"//`);
echo "......................................................" >> HTProcess_Cufflinks/manifest_file.txt
echo "!!!Cufflinks Output files!!!" >> HTProcess_Cufflinks/manifest_file.txt
for ((i = 0; i < "$count"; i += 1));
do
bamfile=`echo "${array[$i]}" | sed 's/,.*//'`
singletest=`echo $bamfile | grep 'TrmS_'`
if [[ -z "${singletest}" ]] ;
then
cufflinksout=`echo $bamfile`
newname="$cufflinksout"'.cfflnx'
command2="$command1 -o CuffOut"
echo "Bamfile equals $bamfile" >> HTProcess_Cufflinks/HTProcess.log
cufflinks -p 5 --no-update-check $command2 HTProcess_BAM*/$bamfile
echo "command: cufflinks -p 5 $command2 --no-update-check HTProcess_BAM*/$bamfile" >> HTProcess_Cufflinks/HTProcess.log
wait
mv CuffOut/genes.fpkm_tracking HTProcess_Cufflinks/"$newname"'_genes.fpkm_tracking'
echo '!!!GENES_FPKM' "$newname"'_genes.fpkm_tracking' >> temp
mv CuffOut/isoforms.fpkm_tracking HTProcess_Cufflinks/"$newname"'_isoforms.fpkm_tracking'
echo '!!!ISOFORMS_FPKM' "$newname"'_isoforms.fpkm_tracking' >> temp
mv CuffOut/transcripts.gtf HTProcess_Cufflinks/$newname'.transcripts.gtf'
mv CuffOut/skipped.gtf HTProcess_Cufflinks/intermediate_files/$newname.'skipped.gtf'
echo "!!!Assembly-GTF_file $newname.transcripts.gtf" >> temp
fi
done

#SINGLES/Orphans
Assembl='!!OSBAM-file'
if [ "${single_assembly}" == 'yes' ]; then Assembl='NOTATHING' ; fi
count=`grep -c $Assembl manifest_file.txt`
array=(`grep $Assembl manifest_file.txt | sed s/"$Assembl"//`);
for ((i = 0; i < "$count"; i += 1));
do
bamfile=`echo "${array[$i]}" | sed 's/,.*//'`

cufflinksout=`echo $bamfile`
newname="$cufflinksout"'.cfflnx'
command2="$command1 -o CuffOut"
echo "Bamfile equals $bamfile" >> HTProcess_Cufflinks/HTProcess.log
cufflinks -p 5 --no-update-check $command2 HTProcess_BAM*/$bamfile
echo "command: cufflinks -p 5 --no-update-check $command2 HTProcess_BAM*/$bamfile" >> HTProcess_Cufflinks/HTProcess.log
mv CuffOut/genes.fpkm_tracking HTProcess_Cufflinks/"$newname"'_genes.fpkm_tracking'
echo '!!!GENES_FPKM' "$newname"'_genes.fpkm_tracking' >> temp
mv CuffOut/isoforms.fpkm_tracking HTProcess_Cufflinks/"$newname"'_isoforms.fpkm_tracking'
echo '!!!ISOFORMS_FPKM' "$newname"'_isoforms.fpkm_tracking' >> temp
mv CuffOut/transcripts.gtf HTProcess_Cufflinks/$newname'.transcripts.gtf'
mv CuffOut/skipped.gtf HTProcess_Cufflinks/intermediate_files/$newname.'skipped.gtf'
echo "!!!Assembly-GTF_file $newname.transcripts.gtf" >> temp
#fi
done
echo "......................................................" >> HTProcess_Cufflinks/manifest_file.txt
echo "" >> HTProcess_Cufflinks/manifest_file.txt
echo "" >> HTProcess_Cufflinks/manifest_file.txt
echo "......................................................" >> HTProcess_Cufflinks/manifest_file.txt
cat temp >> HTProcess_Cufflinks/manifest_file.txt
condition=`grep 'condition' manifest_file.txt | sed 's/condition//' | sed 's/\=//'`
rm manifest_file.txt
rm temp
echo "Finished Running Cufflinks." >> HTProcess_Cufflinks/HTProcess.log
mv HTProcess_Cufflinks HTProcess_Cufflinks'_'"$condition"

