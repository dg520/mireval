miREval 2.0 Command-line Version (Beta)

**DEPENDENCIES**
Befor install miREval you should have the following:
1. Perl5
   Please make sure the perl interpreter is at /usr/bin/perl
2. Python 2.X (>=2.6)
   Please make sure the python interpreter is at /usr/bin/python
3. BLAST+
   Please make sure /usr/bin/blastn is executable
4. Numpy module for Python
5. BioPython
6. BioPerl


**INSTALLATION**
1. Extract mireval.tar.gz:
   > tar -zxvf mireval.tar.gz

2. Set the mireval environment by running:
   > chmod a+x /YOUR_PATH_OF_MIREVAL/*.py
   > chmod a+x /YOUR_PATH_OF_MIREVAL/bin/*.py

3. Download and compile SVM-light:
   > cd /YOUR_PATH_OF_MIREVAL/bin/svm_light
   > make (or make all)
   Check whether the above command compiles the system and creates the two executables: "svm_learn" and "svm_classify". The following command should be executable:
   > /YOUR_PATH_OF_MIREVAL/bin/svm_light/svm_classify

4. Download and install MEME:
   a) 
   > wget http://ebi.edu.au/ftp/software/MEME/4.9.1/meme_4.9.1_1.tar.gz
   > tar -zxvf meme_4.9.1_1.tar.gz
   > cd meme_4.9.1_1
   > ./configure --prefix=/PATH_TO_MIREVAL/mireval/bin/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
   > make
   > make test
   > make install
   b) "fimo" should be executable from:
   > /YOUR_PATH_OF_MIREVAL/bin/meme/bin/fimo

5. check and install perl module needed by Circos-0.64 (NOTE: do not need to install Circos-0.64, it's already with miREval. Just check if you have all the perl module it needs):
   a) check modules by running:
   > /YOUR_PATH_OF_MIREVAL/bin/circos-0.64/bin/test.modules
   b) install the modules you do not have through CPAN:
   > perl -MCPAN -e shell
   ...
   > cpan[1] >install Set::IntSpan
   >install Math::Bezier
   >install MODULES_NAMES_YOU_DONT_HAVE
   ...


**TEST**
Run the following command to test whether miREval is installed successfully.

> cd /YOUR_PATH_OF_MIREVAL
> python ./mireval2 -i example/example.txt -d genomes/caenorhabditis_elegans/blastdb -v genomes/caenorhabditis_elegans/chrom.sizes -o test -m -T genomes/caenorhabditis_elegans/trans/trans.txt -C genomes/caenorhabditis_elegans/cons/score.bw -M genomes/caenorhabditis_elegans/mir/mir.gff3 -l 89 -a 0.77 -b 0.9
> ./mireval2_diagram.py -t test -n 1
> ./mireval2_diagram.py -t test -n 2

This would create a folder called "test" under "/YOUR_PATH_OF_MIREVAL". You can open "index.html" inside the folder to view the prediction result. 

***NOTE***: there is no phylogenetic shadowing information for Caenorhabditis elegans, thus parameters "-S", "-x", "-y" is not applied. 
            -R can be applied if you first run the tool "rfam_chr.py" under "bin/"
            See TUTORIAL for details.
