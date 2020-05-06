---
output:
  word_document: default
  pdf_document: default
  html_document: default
---
Log on to Seawulf and go to the following directory
---------------------------------------------------

    cd /home/shared/inbre-group1/SCiSM/phylogenome
    module load Anaconda3/5.1.0
    source activate 
    # Create a conda environment and then activate it 
    source activate phylogenomics

Prokaryote Taxonomic assignment using Kaiju
-------------------------------------------

Go to <a href="https://github.com/bioinformatics-centre/kaiju/blob/master/README.md" class="uri">https://github.com/bioinformatics-centre/kaiju/blob/master/README.md</a> for detailed instructions. All code regarding Kaiju was taken from their Github page and adapted to our needs accordingly
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Downloading and compiling Kaiju
-------------------------------

Kaiju can be downloaded directly from GitHub either as a compressed archive or using the git command line client:
=================================================================================================================

    cd /home/shared/inbre-group1/SCiSM/phylogenome

    git clone https://github.com/bioinformatics-centre/kaiju.git

This will create the directory kaiju in the current directory. Kaiju is written in C/C++11 for Linux and does not depend on additional libraries. For compiling Kaiju and its associated programs, type:
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    cd kaiju/src
    make

Creating the reference database and index
-----------------------------------------

Before classification of reads, Kaiju’s database index needs to be built from the reference protein database. You can either create a local index based on the currently available data from GenBank, or download one of the indexes used by the Kaiju web server.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

For creating a local index, the program kaiju-makedb in the bin/ directory will download a source database and the taxonomy files from the NCBI FTP server, convert them into a protein database and construct Kaiju’s index (the Burrows-Wheeler transform and the FM-index) in one go.
========================================================================================================================================================================================================================================================================================

The downloaded files are several GB in size. It is therefore recommended to run the program in a directory with at least 100 GB of free space.
==============================================================================================================================================

Example usage:
==============

Make a bash file to run Kaiju using the refseq (Completely assembled and annotated reference genomes of Archaea, Bacteria, and viruses from the NCBI RefSeq database).
----------------------------------------------------------------------------------------------------------------------------------------------------------------------

    kaiju-makedb -s refseq

The taxon identifiers must be contained in the NCBI taxonomy files nodes.dmp and names.dmp. Then, Kaiju’s index is created using the programs kaiju-mkbwt and kaiju-mkfmi. For example, if the database FASTA file is called proteins.faa, then run:
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    kaiju-mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o refseq kaiju_db_refseq.faa
    kaiju-mkfmi refseq 

Make a directory and make a soft link of the files you need to run the taxonomic assignment (all trim and paired end fastq files)
=================================================================================================================================

    mkdir unassembled_reads
    cd unassembled_reads
    ln -s /home/shared/inbre-group1/SCiSM/trim/*_trimpe.fastq.gz /home/shared/inbre-group1/SCiSM/phylogenome/unassembled_reads

Sample Korr
-----------

Make a bash file to run your job
================================

    nano korr_kaiju.sh


    #!/bin/bash
    #SBATCH -t 74:00:00
    #SBATCH --nodes=2 --ntasks-per-node=2
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu

    kaiju -z 6 -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i Korr_all_R1_trimpe.fastq.gz -j Korr_all_R2_trimpe.fastq.gz -o Korr_trimpe_Kaiju.out -E 0.05 

Sumbmit your job
----------------

    sbatch korr_kaiju.sh

Output format: Kaiju will print one line for each read or read pair. The default output format contains three columns separated by tabs.
----------------------------------------------------------------------------------------------------------------------------------------

1.either C or U, indicating whether the read is classified or
unclassified. 2.name of the read 3.NCBI taxon identifier of the assigned
taxon

Creating input file for Krona: The program kaiju2krona can be used to convert Kaiju’s tab-separated output file into a tab-separated text file, which can be imported into Krona. It requires the nodes.dmp and names.dmp files from the NCBI taxonomy for mapping the taxon identifiers from Kaiju’s output to the corresponding taxon names.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i Korr_trimpe_Kaiju.out -o Korr_trimpe_Kaiju.out.krona

The file kaiju.out.krona can then be imported into Krona and converted into an HTML file using Krona’s ktImportText program:
----------------------------------------------------------------------------------------------------------------------------

\#\#We need to install Krona first.

    conda install -c auto krona

    ktImportText -o Korr_trimpe_Kaiju.out.krona.html Korr_trimpe_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/unassembled_reads/Korr_trimpe_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

Repeat the above steps for each sample. Remeber to load the Anaconda module and the activate your conda environment, otherwise the job won\`t run unless you do it within your bash file. either way works equally well.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

\#\#para

    nano para_kaiju.sh


    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=2 --ntasks-per-node=2
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    kaiju -z 6 -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i Para_R1_trimpe.fastq -j Para_R2_trimpe.fastq -o Para_trimpe_Kaiju.out -E 0.05 

    sbatch para_kaiju.sh

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i Para_trimpe_Kaiju.out -o Para_trimpe_Kaiju.out.krona

    ktImportText -o Para_trimpe_Kaiju.out.krona.html Para_trimpe_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/unassembled_reads/Para_trimpe_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

para\_neg
---------

    nano para_neg.sh


    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=2 --ntasks-per-node=2
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    kaiju -z 6 -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i Para_neg_R1_trimpe.fastq -j Para_neg_R2_trimpe.fastq -o Para_neg_trimpe_Kaiju.out -E 0.05 

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i Para_neg_trimpe_Kaiju.out -o Para_neg_trimpe_Kaiju.out.krona

    ktImportText -o Para_neg_trimpe_Kaiju.out.krona.html Para_neg_trimpe_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/unassembled_reads/Para_neg_trimpe_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

Sip
---

    nano sip_kaiju.sh

    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=2 --ntasks-per-node=2
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    kaiju -z 6 -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i Sip_R1_trimpe.fastq -j Sip_R2_trimpe.fastq -o sip_trimpe_Kaiju.out -E 0.05 

    sbatch sip_kaiju.sh

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i sip_trimpe_Kaiju.out -o sip_trimpe_Kaiju.out.krona

    ktImportText -o sip_trimpe_Kaiju.out.krona.html sip_trimpe_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/unassembled_reads/sip_trimpe_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

sk
--

    ln -s /home/shared/inbre-group1/SCiSM/trim/Sk_R1_trimpe.fastq /home/shared/inbre-group1/SCiSM/trim/Sk_R2_trimpe.fastq .

    nano sk_kaiju.sh


    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    kaiju -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i Sk_R1_trimpe.fastq -j Sk_R2_trimpe.fastq -o sk_trimpe_Kaiju.out -E 0.05 

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i sk_trimpe_Kaiju.out -o sk_trimpe_Kaiju.out.krona

    ktImportText -o sk_trimpe_Kaiju.out.krona.html sk_trimpe_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/unassembled_reads/sk_trimpe_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

Porkaryotic taxonomic assignment on the Assembled genomes
---------------------------------------------------------

Korr
----

    mkdir ../assembled_reads
    cd ../assembled_reads

    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/korr/Korr_trim_assembly.fa .
    nano korr_assembled_kaiju.sh

    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu

    kaiju -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i Korr_trim_assembly.fa -o Korr_trim_assembly_Kaiju.out -E 0.05 

    sbatch korr_assembled_kaiju.sh

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i Korr_trim_assembly_Kaiju.out -o Korr_trim_assembly_Kaiju.out.krona

    ktImportText -o Korr_trim_assembly_Kaiju.out.krona.html Korr_trim_assembly_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads/Korr_trim_assembly_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

Para
----

    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/para/para_trim_assembly.fa .
    nano para_assembled_kaiju.sh


    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu

    kaiju -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i para_trim_assembly.fa -o Para_trim_assembly_Kaiju.out -E 0.05 

    sbatch para_assembled_kaiju.sh

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i Para_trim_assembly_Kaiju.out -o Para_trim_assembly_Kaiju.out.krona

    ktImportText -o Para_trim_assembly_Kaiju.out.krona.html Para_trim_assembly_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads/Para_trim_assembly_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

Sip
---

    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/sip/sip.final.contigs.fa
    nano sip_assembled_kaiju.sh


    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu

    kaiju -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i sip.final.contigs.fa -o sip_trim_assembly_Kaiju.out -E 0.05 

    sbatch para_assembled_kaiju.sh

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i sip_trim_assembly_Kaiju.out -o sip_trim_assembly_Kaiju.out.krona

    ktImportText -o sip_trim_assembly_Kaiju.out.krona.html sip_trim_assembly_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads/sip_trim_assembly_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

sk
--

    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/sk/sk.final.contigs.fa
    nano sk_assembled_IDBA_kaiju.sh


    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu

    kaiju -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -f /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/refseq/refseq.fmi -i sk.final.contigs.fa -o sk.final.contigs.fa_IDBA_Kaiju.out -E 0.05 

    sbatch sk_assembled_IDBA_kaiju.sh

    kaiju2krona -t /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/nodes.dmp -n /home/shared/inbre-group1/SCiSM/phylogenome/kaijudb/names.dmp -i sk.final.contigs.fa_IDBA_Kaiju.out -o sk.final.contigs.fa_IDBA_Kaiju.out.krona

    ktImportText -o sk.final.contigs.fa_IDBA_Kaiju.out.krona.html sk.final.contigs.fa_IDBA_Kaiju.out.krona

Download html file to local
---------------------------

    scp -r matias_gomez@seawulf.uri.edu:/home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads/sk.final.contigs.fa_IDBA_Kaiju.out.krona.html /Users/matiasgcvc/Desktop

PHYLOGENOMICS
-------------

Making phylogenetic inferences from assembled genomes 16s only for metagenome
-----------------------------------------------------------------------------

<a href="https://github.com/HRGV/phyloFlash" class="uri">https://github.com/HRGV/phyloFlash</a>
-----------------------------------------------------------------------------------------------

<a href="https://github.com/HRGV/Marmics_Metagenomics/wiki/Module-4---Metagenome-quality-control-and-compositional-analysis" class="uri">https://github.com/HRGV/Marmics_Metagenomics/wiki/Module-4---Metagenome-quality-control-and-compositional-analysis</a>
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    cd /home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads
    module load Anaconda3/5.1.0 
    source activate phylogenomics

    conda install phyloflash
    mkdir phyloflash
    cd phyloflash/

The job failed at a point that have not beeen able to make sense of, I\`ll try to create and environment for phylo flash only.
------------------------------------------------------------------------------------------------------------------------------

    module unload Anaconda3/5.1.0 
    module load Anaconda3/5.3.0 
    conda create -n phyloflash
    source activate phyloflash
    conda install phyloflash

It failed again
---------------

See <a href="https://github.com/HRGV/phyloFlash/issues/90" class="uri">https://github.com/HRGV/phyloFlash/issues/90</a> to solve the issue. I had to edit the script to eliminate the line of code that generated the error. We\`ll see if it works or not
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

\#\#Downloading database automatically

\#\#To create a suitable database, just run

    phyloFlash_makedb.pl --remote

Now trying to run a job on the unassembled reads as required by phyloflash
--------------------------------------------------------------------------

    cd unassembled_reads/
    mkdir 16s_phyloflash
    cd 16s_phyloflash/

Korr
----

     mkdir korr
     cd korr

    nano korr_phyloflash.sh

    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    phyloFlash.pl -dbhome /home/shared/inbre-group1/SCiSM/phylogenome/phyloflash/138 -lib korr -read1 /home/shared/inbre-group1/SCiSM/trim/Korr_all_R1_trimpe.fastq.gz -read2 /home/shared/inbre-group1/SCiSM/trim/Korr_all_R2_trimpe.fastq.gz -CPUs 10 -readlength 150 -almosteverything 

    sbatch korr_phyloflash.sh 

Para
----

     mkdir para 
     cd para

    nano para_phyloflash.sh

    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    phyloFlash.pl -dbhome /home/shared/inbre-group1/SCiSM/phylogenome/phyloflash/138 -lib para -read1 /home/shared/inbre-group1/SCiSM/trim/Para_R1_trimpe.fastq.gz -read2 /home/shared/inbre-group1/SCiSM/trim/Para_R2_trimpe.fastq.gz -CPUs 10 -readlength 150 -almosteverything 

    sbatch para_phyloflash.sh 

Sip
---

     mkdir sip
     cd sip

    nano sip_phyloflash.sh

    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    phyloFlash.pl -dbhome /home/shared/inbre-group1/SCiSM/phylogenome/phyloflash/138 -lib sip -read1 /home/shared/inbre-group1/SCiSM/trim/Sip_R1_trimpe.fastq.gz -read2 /home/shared/inbre-group1/SCiSM/trim/Sip_R2_trimpe.fastq.gz -CPUs 10 -readlength 150 -almosteverything 

    sbatch sip_phyloflash.sh 

SK
--

     mkdir sk
     cd sk

    nano sk_phyloflash.sh

    #!/bin/bash
    #SBATCH -t 72:00:00
    #SBATCH --nodes=1 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    phyloFlash.pl -dbhome /home/shared/inbre-group1/SCiSM/phylogenome/phyloflash/138 -lib sk -read1 /home/shared/inbre-group1/SCiSM/trim/Sk_R1_trimpe.fastq.gz -read2 /home/shared/inbre-group1/SCiSM/trim/Sk_R2_trimpe.fastq.gz -CPUs 24 -readlength 150 -almosteverything 

    sbatch sk_phyloflash.sh

Although the program failed because of the errors shown here (<a href="https://github.com/HRGV/phyloFlash/issues/90" class="uri">https://github.com/HRGV/phyloFlash/issues/90</a>), we got the files that are neeeded to do phylogenies based on assembled SSUs from all samples
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

<a href="https://github.com/HRGV/Marmics_Metagenomics/wiki/Module-4---Metagenome-quality-control-and-compositional-analysis" class="uri">https://github.com/HRGV/Marmics_Metagenomics/wiki/Module-4---Metagenome-quality-control-and-compositional-analysis</a>
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Phylogeny based comparison using phyloFlash assembled SSUs
----------------------------------------------------------

Another way to complement the barplots based on the hit counts is to generate a collective phylogenetic tree of the assembled SSU sequences and their reference sequences. This aproach does not need R and might be easier to implement. For this, you need to collect all SSU assemblies, collect the references, deduplicate the references, align the sequences and generate a tree. A quick workflow for this based on the easy to install or already available tools uses vsearch, mafft and fasttree to generate this tree in nwk format that can be explored e.g. with itol. This workflow assumes that you have already collected all zipped phyloFlash result files you want to compare in a single folder.
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Let’s extract all the spades assemblies of the SSU genes from the phyloFlash zips as well as the hits in the SILVA database
---------------------------------------------------------------------------------------------------------------------------

    cat korr/*spades_rRNAs.final.fasta para/*spades_rRNAs.final.fasta sip/*spades_rRNAs.final.fasta sk/*spades_rRNAs.final.fasta > pf_collection.fasta

    cat korr/*all.dbhits.NR97.fa para/*all.dbhits.NR97.fa sip/*all.dbhits.NR97.fa sk/*all.dbhits.NR97.fa > pf_collection_dbhits.fasta

The database hits will be redundant, so let’s cluster them at approximately species level and only keep one per species level clade
-----------------------------------------------------------------------------------------------------------------------------------

    vsearch --cluster_smallmem pf_collection_dbhits.fasta --usersort --id 0.98 --centroids pf_collection_dbits.NR98.fasta --notrunclabels

Now we can collect the non-redundant reference hits and the spades assembled SSUs
---------------------------------------------------------------------------------

      cat pf_collection_dbits.NR98.fasta pf_collection.fasta > pf_collection_dbhits_and_spades.fasta

We have to align the sequences, mafft should do the trick quite well. Assuming that you have &lt;100 sequences to compare, mafft-xinsi that considers the secondary structure of the SSU likely would be the best choice, just replace mafft with mafft-xinsi in the command and output label…
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

      mafft --adjustdirection --thread -1 pf_collection_dbhits_and_spades.fasta > pf_collection_dbhits_and_spades_mafft.fasta

fasttree -nt -gtr -quote pf\_collection\_dbhits\_and\_spades\_mafft.fasta &gt; pf\_collection\_dbhits\_and\_spades\_mafft\_fasttree.nwk
---------------------------------------------------------------------------------------------------------------------------------------

    conda install -c bioconda fasttree  
      
      fasttree -nt -gtr -quote pf_collection_dbhits_and_spades_mafft.fasta > pf_collection_dbhits_and_spades_mafft_fasttree.nwk

Copy from server to local (do it from a local terminal window)
--------------------------------------------------------------

    scp -r matias_gomez@seawulf.uri.edu://home/shared/inbre-group1/SCiSM/phylogenome/unassembled_reads/16s_phyloflash/pf_collection_dbhits_and_spades_mafft_fasttree.nwk /Users/matiasgcvc/Desktop

Load file in R to visualize it
------------------------------

    library(ape)
    tree<-read.tree(file = "pf_collection_dbhits_and_spades_mafft_fasttree.nwk")

    plot.phylo(tree, cex = 0.5, font = 2, adj = 0, edge.width = 1, show.tip.label = F )

    plot.phylo(tree, cex = 0.7, font = 2, adj = 0, edge.width = 1, show.tip.label = T )

    plot.phylo(tree, cex = 0.7, font = 2, adj = 0, edge.width = 1, type = "unrooted", show.tip.label = F  )

    tree$tip.label

Then you can use a tree visualization program to edit it more easily or keep working in R. I used Itol webserver to make better looking trees: <a href="https://itol.embl.de/" class="uri">https://itol.embl.de/</a>
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Ciliate whole metagenome phylogenomics
--------------------------------------

REALPHY
-------

The REALPHY pipeline requires the user to provide a set of DNA sequences for each taxon to be included in the phylo- genetic tree. This set will typically consist of short-sequence reads but may also include larger sequences, such as fully or partially assembled genomes.
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

\#The program requires several other programs in order to run. These
are: JAVA downloadable from:
<a href="http://www.java.com/" class="uri">http://www.java.com/</a>
(works with version 1.6.0\_20 64bit) Bowtie2 downloadable from:
<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" class="uri">http://bowtie-bio.sourceforge.net/bowtie2/index.shtml</a>
(works with version 2.0.0-beta7)

\#A few more programs can be run by REALPHY if desired (Note: Either
TREE-PUZZLE, RAxML, PHYML or dnapars is needed to build a tree):

    module load Java/11.0.2
    module load Anaconda3/5.1.0
    source activate phylogenomics
    conda install -c bioconda realphy
    conda install -c bioconda bowtie2
    conda install -c bioconda raxml
    conda install -c bioconda samtools

Download a closely related ciliate genome to use as a reference for phylogenomics analysis
------------------------------------------------------------------------------------------

\#<a href="https://www.ncbi.nlm.nih.gov/Traces/wgs/AFXY02?display=contigs" class="uri">https://www.ncbi.nlm.nih.gov/Traces/wgs/AFXY02?display=contigs</a>
\#<a href="https://www.ncbi.nlm.nih.gov/Traces/wgs/AFXF02?display=contigs" class="uri">https://www.ncbi.nlm.nih.gov/Traces/wgs/AFXF02?display=contigs</a>

    cd /home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads
    wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/AF/XF/AFXF02/AFXF02.1.fsa_nt.gz
    gunzip AFXF02.1.fsa_nt.gz
    mkdir realphy
    cd realphy
    mv ../AFXF02.1.fsa_nt .
    mv AFXF02.1.fsa_nt AFXF02.1.fsa_nt.fa

    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/korr/korr_megahit_result/korr.final.contigs.fa .
    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/sip/assembly_results/sip_contig.fa .
    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/sip/sip_megahit_result/sip.final.contigs.fa .
    ln -s  /home/shared/inbre-group1/SCiSM/assembly/IDBA/sk/sk_megahit_result/sk.final.contigs.fa

    mkdir results

Soft link to assembled genomes to working directory
===================================================


    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/korr/korr_megahit_result/korr.final.contigs.fa .
    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/sip/assembly_results/sip_contig.fa .
    ln -s /home/shared/inbre-group1/SCiSM/assembly/IDBA/sip/sip_megahit_result/sip.final.contigs.fa .
    ln -s  /home/shared/inbre-group1/SCiSM/assembly/IDBA/sk/sk_megahit_result/sk.final.contigs.fa

The above programs can be installed anywhere on the computer, however the location and name of the executable needs to be specified in a config file called config.txt.
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

The config file then needs to be placed into the designated output folder (which is specified as a command line parameter) before REALPHY can be run.
-----------------------------------------------------------------------------------------------------------------------------------------------------

    cd results
    nano config.txt

    BOWTIE2 ~/.conda/envs/phylogenomics/bin/bowtie2
    RAXML ~/.conda/envs/phylogenomics/bin/raxmlHPC

Usage: java -Xmx\[available RAM in MB\]m -jar RealPhy\_v1.12.jar \[Sequence folder\] \[Output folder\] \[Options\]
==================================================================================================================

\#Sequence folder needs to contain fasta files ending with .fas, .fna,
.fasta or .fa, genbank files ending in .gbk or .gb and short read files
in fastq format ending in .fastq or fastq.gz. The output folder needs to
contain a file called “config.txt”, which contains information about the
location of the required executables such as bowtie2.

Options:
========

\#-ref \[sequence file name (without extension or path!)\]
default=random; Possible values: The file name of a sequence data set
without the extension (.fas, .fasta, .fa, .fna, .fastq, .fastq.gz, .gb
or .gbk). Sets the reference sequence.

-treeBuilder \[integer\] default=4;
-----------------------------------

0=Do not build a tree;
----------------------

1=treepuzzle;
-------------

2=raxml
-------

3=max. parsimony (dnapars)
--------------------------

4=PhyML (default)
-----------------

-config \[string\] this specifies the location of the config.txt. If not set it is assumed that the config.txt is in the working directory
------------------------------------------------------------------------------------------------------------------------------------------


    nano realphy.sh

    #!/bin/bash
    #SBATCH -t 200:00:00
    #SBATCH --nodes=4 --ntasks-per-node=1
    #SBATCH --mail-type=BEGIN  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=END  --mail-user=matias_gomez@uri.edu
    #SBATCH --mail-type=FAIL  --mail-user=matias_gomez@uri.edu
    REALPHY_v112 /home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads/realphy /home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads/realphy/results -ref AFXF02.1.fsa_nt -treeBuilder 2 -config /home/shared/inbre-group1/SCiSM/phylogenome/assembled_reads/realphy/results

    sbatch --mem-per-cpu=40G realphy.sh

The job failed due to Java running out of memory. Whoever continues with this pipeline, you\`ll have to solve this issue. Sorry about this and good luck.
---------------------------------------------------------------------------------------------------------------------------------------------------------
