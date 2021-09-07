*miREval 2.0 Command-line Version (Beta)*

**MIREVAL REFERENCE FILES**    
miREval needs several parameters and reference files to make the full use. Most of them can be downloaded from the corresponding online database. miREval also has small tools to help user handling the raw reference and make them miREval suitable.     
    
I) -d <NAME_OF_BLAST_DATABASE>       
BLAST Databese is made using BLAST+ toolkit. First you need to have a FASTA file for the whole genome. You can either    
  1) Download the reference genome from UCSC (See http://hgdownload.soe.ucsc.edu/downloads.html)     
  2) Or use you own specific reference genome    
     
Before building the BLAST Databese. PLEASE make sure the following:    
  1) check the FASTA title for your genome alignment. The title should ONLY contain the CHROMOSOME NAME.(e.g.: "chrM" or "M" are both acceptable)     
     ***NOTE***: title in other format such like NCBI format (e.g.:"gi|1110804074|ref|NC_004193.1") will probably raise error when you want to compare you inquiry sequence with miRBase or rFam data.    
  2) Sometimes the reference is by chromosome set (one fasta file per chromosome). In this case, you need to merge them into a single FASTA file.     
     
Now you can build the BLAST database for your interested genome using BLAST+ toolkit by running the following:    
```
  > formatdb -p F -i MERGED_FASTA_FILE -n NAME_OF_BLAST_DATABASE -o T
```
     
II) -v <Chromosome-size>    
miREval need to know the size of each chromosome.       
  1) If your genome is downloaded from UCSC, this information can be easily got using "fetchChromSizes" under "bin/" of miREval:     
```
  > /YOUR_PATH/mireval/bin/fetchChromSizes <db> > OUTPUT_FILE_NAME 
```  
  ***NOTE***: <db> is the name of UCSC database (e.g.:hg19, mm10, etc...)    
     
  2) If you're using your own genome, you can use "countChromSizes.py" under "bin/" of miREval:
```    
  > /YOUR_PATH/mireval/bin/countChromSizes.py -i MERGED_FASTA_FILE -o OUTPUT_FILE_NAME
```
    
III) -T <transcripts_table>     
miREval will look whether the inquiry overlaps any transcripts annotated in this table. The table should be in UCSC table format.       
  1) For a genome recorded in UCSC, one can use the table browser to download this table.      
     a) In the UCSC table browser, select your genome and assembly     
     b) group: Gene and Gene Predictions; track: RefSeq Genes     
     c) output format: all fields from selected table; you can give an output file name and click "get output"    
          
  2) If you're using your own genome, you can make your own transcripts table of the same format (example: /YOUR_PATH/mireval/genomes/caenorhabditis_elegans/trans/trans.txt).      
  ***NOTE***: if you don't have infomation of the columns related to exons, you can put 0 there.     
     
     
IV) -C <bigwig>, -S <bigwig>     
miREval can use the evidence of evolutionary conservation to facilitate miRNA prediction.      
There are two kinds of conservation scores: phylogenetic conservation and phylogenetic shadowing. The former one reflects the conservation between a bundle of organisms. The latter one is the convervation among very closely related species. These two kinds of data are fed to miREval by -C (for phylogenetic conservation) and -S (for phylogenetic shadowing) respecitvely.    
     
Both kinds of conservation data can be downloaded from UCSC (See http://hgdownload.soe.ucsc.edu/downloads.html) through the link "Basewise conservation scores (phyloP) of N genomes with the species".     

***NOTE***: miREval needs conservation socres to be in bigwig format. For some species on UCSC, the scores are in wigFix format and some of them are annotated per chromosome. Please do the following to convert them:     
     
  1) download score files into a certain folder    
  2) In that folder, merge all the per-chromosome score to a single wigFix file by running:    
```  
  >  zcat *.gz | gzip -c > score.wigFix.gz
```
  3) convert wigFix to bigwig using the tool "wigToBigWig" under "bin/" of miREval (the chromosome-size file mentioned in Part II is required):    
```
  > /YOUR_PATH/mireval/bin/wigToBigWig /PATH_TO_MERGED_WIGFIX/score.wigFix.gz CHROMOSOME_SIZE_FILE OUTPUT_BIGWIG
```
    
V) -M <miRBase_gff3_file>     
miREval can use a list of known miRNA of the species to facilitate miRNA prediction. This file should be in gff3 format of miRBase.     
You can download the gff3 file of your insterested genome from miRBase (ftp://mirbase.org/pub/mirbase/CURRENT/genomes/)     
      
***NOTE***: The gff3 annotation you downloaed must use the same version of assembly of your BLAST database file. The assembly information of gff3 annotation is at the beginning of the file. If the assembly is not the same as your BLAST database file, you can either:     
  1) check the previous miRBase version at ftp://mirbase.org/pub/mirbase/     
  2) create your own known miRNA annotation for your insteresed species according to a typical miRBase gff3 one.    
     
     
VI) -a <float_number>, -b <float_number>, -c <float_number>, -d <float_number>     
When both conservation scores (either phylogenetic conservation or phylogenetic shadowing) and known miRNA gff3 annoatation are avaible, miRBase can facilitate miRNA predict by determine whether the inquiry sequence is above a certain conservation thresold (See http://mimirna.centenary.org.au/mireval/tutorial/circle.html for details).    
     
There are two kinds of thresolds for phylogenetic conservation or phylogenetic shadowing, respectively.      
  1) -a: upper thresold of phylogenetic conservation score, determined by the average phylogenetic conservation score of MATURE miRNAs of the species;     
  2) -b: lower thresold of phylogenetic conservation score, determined by the average phylogenetic conservation score of miRNA PRECURSORs of the species;    
  3) -x: upper thresold of phylogenetic shadowing score, determined by the average phylogenetic shadowing score of MATURE miRNAs of the species;     
  4) -y: lower thresold of phylogenetic shadowing score, determined by the average phylogenetic shadowing score of miRNA PRECURSORs of the species;     
     
You can use the tool "cons_thresold.py" under "bin/" of miREval to calculate these values:    
```
  > /YOUR_PATH/mireval/bin/cons_thresold.py -v CHROMOSOME_SIZE_FILE -C PHYLOGENETIC_CONSERVATION_BIGWIG -S PHYLOGENETIC_SHADOWING_BIGWIG -M MIRBASE_GFF3
```
     
It will print out 4 lines on the screen like the following, whose value are for -a, -b, -x, -y respectively:     
***NOTE***: if any of the 4 values is "NA", DO NOT use that corresponding parameter.      
      
      
VII) -R <fRam_location_index>     
miREval can search the ncRNA database Rfam to facilitate miRNA prediction. The database is included in the miREval package.     
Rfam annotation is a bundle of gff3 files that are not classfied by species. In order to use the correct annotation, an index file indicating which files should be used for the interested species must be made beforehand.     
      
This can be done by using the tool "rfam_chr.py" under "bin/" of miREval:
```
  > /YOUR_PATH/mireval/bin/rfam_chr.py -s SPECIES_NAME -o OUTPUT_FILE
```  
  ***NOTE***: 1) SPECIES_NAME should be expressed like this format: mus_musculus, homo_sapiens, busseola_fusca, etc...     
              2) the output file might be empty which means the species is not annotated by Rfam.     
       
       
**OUTPUT FILES**
miREval can be run in two different manner:     
  1) a "light" version only predicts the potential miRNA within inquiry sequences without other genomic information:    
```
  > /YOUR_PATH/mireval/mireval2 -i INQUIRY_FASTA_FILE -o OUTPUT_FOLDER
```
     
  2) a full version predicting the potential miRNA within inquiry sequences and ploting other genomic information:    
```
  > /YOUR_PATH/mireval/mireval2 -i INQUIRY_FASTA_FILE -d BLAST_DATABASE -v CHROMOSOME_SIZE -o OUTPUT_FOLDER [Opt...]
```
    
miREval output is saved to the folder difined by "-o". You can directly open the html file "index.html" within your output folder. A heatmap indicates you the prediction result for each inquiry sequence. Generally speaking, the deeper the red, the more likely the inquiry contains potential miRNA (See images in "help" folders for details).      
Depending on how you run the "mireval2", the last step is a little bit different:     
  1) If you run the "light" version of miREval (without "-d" option), you can directly open the link of inquiry sequence which you're insterested in.     
        
  2) If miREval is applied with a reference genome, you might also want to visualize the genomic elements around a certain inquiry sequence. Before you click the link of that sequence you can run the following command:     
```
  > /YOUR_PATH/mireval/mireval2-diagram -t MIREVAL_RESULT_FOLDER -n THE_NUMBER_OF_YOUR_INSTERESTED_SEQ
```  
  e.g.: if the "-o" option when you run "mireval2" is set to "/YOUR_PATH/result1" and you're insterested in the 3rd inquiry within your input FASTA file, you should run:
```
  > /YOUR_PATH/mireval/mireval2-diagram -t /YOUR_PATH/result1 -n 3
```
  and then you can click the link to view your result.


