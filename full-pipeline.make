SHELL=/bin/bash -o pipefail

##################################################
#
# Step 0. Preamble: set up paths and variable
#
##################################################

# do not leave failed files around
.DELETE_ON_ERROR:

# The programs we will install must be on the PATH
export PATH := ./DAZZ_DB:./DALIGNER:./nanocorrect:./poaV2:./wgs-8.2/Linux-amd64/bin/:./samtools/:$(PATH)

# Download links for programs that are not on github
CA_LINK=http://downloads.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-8.2/wgs-8.2-Linux_amd64.tar.bz2
POA_LINK=http://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz

# Parameters to control execution
NC_PROCESS=8
NP_PROCESS=8
THREADS=4

##################################################
#
# Step 1. Download the input data from the ENA
# and unpack it into fast5 files
#
##################################################

# Download rules
ERX708228.tar:
	wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3.tar -O ERX708228.tar

ERX708229.tar:
	wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_32_LomanLabz_K12_His_tag.tar -O ERX708229.tar

ERX708230.tar:
	wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag.tar -O ERX708230.tar

ERX708231.tar:
	wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_39.tar -O ERX708231.tar

# Untar into directories named x.fast5 and move the files into the top directory
%.fast5/: %.tar
	mkdir -p $@
	cd $@; tar -xf ../$<
	find $@ -name "*.fast5" -exec mv {} $@ \;

##################################################
#
# Step 2. Download and install the programs
# required to run the assembly
#
##################################################

# Install nanocorrect & dependencies
nanocorrect.version:
	git clone https://github.com/jts/nanocorrect.git
	ln -s nanocorrect/poa-blosum80.mat
	-cd nanocorrect; git log | head -1 > ../$@

# Install Python dependencies
pythonlibs.version:
	pip install pysam > $@
	pip install cython >> $@
	pip install numpy==1.8.1 >> $@
	pip install h5py==2.3.0 >> $@
	pip install cython >> $@
	pip install poretools >> $@
	pip install biopython >> $@

# Install poa
poa.version:
	wget $(POA_LINK)
	tar -xzf poaV2.tar.gz
	cd poaV2; make CFLAGS='-O3 -g -DUSE_WEIGHTED_LINKS -DUSE_PROJECT_HEADER -I.' poa
	ln -s poaV2/poa
	echo $(POA_LINK) > $@

# Install DALIGNER
daligner.version:
	git clone https://github.com/thegenemyers/DALIGNER.git
	cd DALIGNER; git checkout 549da77b91395dd; make
	echo "549da77b91395dd" > $@

# Install DAZZ_DB
dazz_db.version:
	git clone https://github.com/thegenemyers/DAZZ_DB
	cd DAZZ_DB; git checkout 8cb2f29c4011a2c2; make
	echo "8cb2f29c4011a2c2" > $@

# Install Celera Assembler
ca.version:
	wget $(CA_LINK)
	tar -xjf wgs-8.2-Linux_amd64.tar.bz2
	wget http://www.cbcb.umd.edu/software/PBcR/data/convertFastaAndQualToFastq.jar
	mv convertFastaAndQualToFastq.jar wgs-8.2/Linux-amd64/bin/
	echo $(CA_LINK) > $@

# Install samtools
samtools.version:
	git clone --recursive https://github.com/samtools/htslib.git
	cd htslib; make
	git clone --recursive https://github.com/samtools/samtools.git
	cd samtools; make
	-cd samtools; git log | head -1 > ../$@

# Install nanopolish, automatically downloading libhdf5
nanopolish.version:
	git clone --recursive https://github.com/jts/nanopolish.git
	cd nanopolish; make libhdf5.install nanopolish
	-cd nanopolish; git log | head -1 > ../$@

# Install bwa
bwa.version:
	git clone https://github.com/lh3/bwa.git
	cd bwa; make
	-cd bwa; git log | head -1 > ../$@

#
test_path: FORCE
	-LAshow 2> test
	-fasta2DB 2>> test
	-poa 2>> test
FORCE:

##################################################
#
# Step 3. Correct the raw nanopore data using
# nanocorrect. This step takes a long time.
#
##################################################

# Export 2D reads to fasta files using poretools
raw.reads.fasta: ERX708228.fast5 ERX708229.fast5 ERX708230.fast5 ERX708231.fast5 pythonlibs.version
	poretools fasta --type 2D ERX708228.fast5/ > $@
	poretools fasta --type 2D ERX708229.fast5/ >> $@
	poretools fasta --type 2D ERX708230.fast5/ >> $@
	poretools fasta --type 2D ERX708231.fast5/ >> $@

# Run nanocorrect in parallel
%.corrected.fasta: %.fasta samtools.version pythonlibs.version nanocorrect.version daligner.version dazz_db.version poa.version
	make -f nanocorrect/nanocorrect-overlap.make INPUT=$< NAME=$*
	samtools faidx $<
	python nanocorrect/makerange.py $< | parallel --progress -P $(NC_PROCESS) 'python nanocorrect/nanocorrect.py $* {} > $*.{}.corrected.fasta'
	cat $*.*.corrected.fasta > $@
	#rm $*.*.corrected.fasta

##################################################
#
# Step 4. Run the celera assembler on the 
# corrected reads.
#
##################################################

# prepare the input into celera assembler. 
# we want to use the twice-corrected data here, so two prereqs
assembly.input.fastq: raw.reads.corrected.fasta raw.reads.corrected.corrected.fasta ca.version
	java -Xmx1024M -jar ./wgs-8.2/Linux-amd64/bin/convertFastaAndQualToFastq.jar raw.reads.corrected.corrected.fasta > $@

assembly.frg: assembly.input.fastq
	fastqToCA -technology sanger -libraryname assembly -reads $< > $@

# Download the spec file we use
trimmedends_ovlErrorRate0.06.spec:
	wget --no-check-certificate https://raw.githubusercontent.com/jts/nanopore-paper-analysis/936e98732c92d7432d369adf079dda1fd436dcab/trimmedends_ovlErrorRate0.06.spec

# Run the assembly
celera-assembly/9-terminator/asm.scf.fasta: trimmedends_ovlErrorRate0.06.spec assembly.frg
	runCA -d celera-assembly -p asm -s trimmedends_ovlErrorRate0.06.spec assembly.frg

draft_genome.fasta: celera-assembly/9-terminator/asm.scf.fasta
	ln -s $< $@

##################################################
#
# Step 5. Polish the draft assembly with nanopolish
#
##################################################

# preprocess the fasta file for nanopolish
raw.reads.np.fasta: raw.reads.fasta nanopolish.version
	nanopolish/consensus-preprocess.pl $< > $@

# index the draft assembly for bwa
draft_genome.fasta.bwt: draft_genome.fasta
	bwa index $<

# index the draft assembly for faidx
draft_genome.fasta.fai: draft_genome.fasta
	samtools faidx $<

# align reads to draft assembly
reads_to_draft.sorted.bam: draft_genome.fasta draft_genome.fasta.bwt raw.reads.np.fasta bwa.version samtools.version
	bwa/bwa mem -t $(THREADS) -x ont2d draft_genome.fasta raw.reads.np.fasta | samtools view -Sb - | samtools sort - $@

# index the bam file
reads_to_draft.sorted.bam.bai: reads_to_draft.sorted.bam
	samtools/samtools index $<

# run nanopolish
polished_genome.fasta: draft_genome.fasta draft_genome.fasta.fai reads_to_draft.sorted.bam.bai raw.reads.np.fasta
	python nanopolish/nanopolish_makerange.py draft_genome.fasta | parallel --progress -P $(NP_PROCESS) \
        nanopolish/nanopolish consensus -o nanopolish.{1}.fa -r raw.reads.np.fasta -b reads_to_draft.sorted.bam -g draft_genome.fasta -w {1} -t $(THREADS)
	python nanopolish/nanopolish_merge.py draft_genome.fasta nanopolish.scf*.fa > $@
