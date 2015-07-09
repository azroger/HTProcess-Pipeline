# HTProcess-Pipeline
This is a repository of my scripts for running HTProcess apps in the iPlant Collaborative's Discovery Environment. This
pipeline is intended to simplify analysis of multiple read files at once. All analysis through this pipeline can be done in
whole or in part in the pipeline, but all HTProcess apps require that all previous steps be run in the pipeline. The first step
required is HTProcess-prepare_directories_and_run_fastqc, which is a Discovery Environment (DE) workflow, i.e. 2 separate apps
connected as a workflow. It starts by creating the necessary directories with the app
HTProcess-prepare directories, and the directories are used as inputs for HTProcess_fastqc, which runs FastQC from
Brabraham Bioinformatics (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/). This app and the combined workflow create
the "HTProcess_Reads" directory, the output to be used as the obligate input in the next application, HTProcess_trimmomatic. All
HTProcess apps require usually a complete directory coming from a previous HTProcess step. The critical component in each
such directory is the manifest file (manifest_file.txt), which documents the contents of the directory. Entries in the first
manifest file are derived from the required inputs for the first critical step (HTProcess_fastqc or
HTProcess-prepare_directories_and_run_fastqc) and from the outputs from FastQC. Some of these entries must be accurate for
proper operation of downstream HTProcess pipeline apps.
