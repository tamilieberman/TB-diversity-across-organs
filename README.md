
TB-diversity-across-organs
================


Code and intermediate files for reproducing the figures from "Genomic diversity in autopsy samples reveals within-host dissemination of HIV-associated M. tuberculosis ", by Tami D Lieberman, Douglas Wilson, Reshma Misra, Lealia L Xiong, Prashini Moodley, Ted Cohen, and Roy Kishony
 <br> <br>

If you are applying these scripts to your own data, it is strongly recommended that you investigate your raw data carefully and adjust parameters to suit your genome, coverage, etc. If you are having trouble understanding these scripts or how to use them, would like additional functionality/flexibility, or have other questions, feel free to contact Tami Lieberman (Email address is easy to find at other locations.). <br>

If you find any of these scripts helpful, please cite: <br>
XXXXXX


Option 1 for reproducing the figures (recommended): Start with pre-processed data: <br>
------------------------------------------------------------

1) Download subject_folders, containing a directory for each subject, each containing candidate_mutation_table.mat <br>
2) Download the scripts and tools directories <br>
3) Run identify_de_novo_muts.m to call de novo mutations in each subject (and gather other information for each subject), producing de_novo_muts.mat for each subject <br>
4) Run analyze_de_novo_muts.m to analyze these mutations and generate the figures <br>


Option 2 for reproducing the figures: Start with called mutations: <br>
------------------------------------------------------------

Run step 4 only from above, ensuring de_novo_muts.mat are downloaded for each subject in subject_folders. <br>


Option 3: Generate the processed data <br>
------------------------------------------------------------
1) Download raw fastq files for each sample from the SRA (BioProject PRJNA323744) <br>
2) Process each fastq file using the commands in for_making_vcf_files.txt, processing each sample in its own directory <br>
3) Modify the sample_names.csv file for each subject in subject_folders to point to the correct location of the processed files. <br>
4) Run build_candidate_mutation_table_tb(SUBJECTFOLDER) for each subject <br>
5) Continue as in Option 1 <br>


Functions that you may find useful for your other data<br>
------------------------------------------------------------
*build_candidate_mutation_table_tb.m* <br>
Starting with a  .pileup file (summarizing the alignment file in a text format) and .vcf file (must contain FQ scores for EVERY position on the genome — see for_making_vcf_files.txt) grabs useful information from each potential variant position. Adjust parameters at the top of this file to suit your needs. <br> <br>

*identify_de_novo_muts.m* <br>
Starting with the output of build_candidate_mutation_table_tb, identify candidate mutations. <br> <br>

*clickable_snp_table.m* <br>
Use to investigate the raw data at each genomic position. Requires data structures generated in identify_de_novo_muts.m <br> <br>

*find_genotypes.m* and *assign_genotypes.m* <br>
Assigns genotypes given a matrix of derived mutation frequencies across samples <br> <br>
