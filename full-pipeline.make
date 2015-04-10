SHELL=/bin/bash -o pipefail

all: poretools.version \
     nanocorrect.version \
     daligner.version \
     poa.version \
     dazz_db.version \
     samtools.version

# These programs need to be a on the PATH for everything to work
export PATH := ./DAZZ_DB:./DALIGNER:./nanocorrect:./poaV2:$(PATH)

# Paths to programs that come as tar balls rather than git repos
CA_PATH=http://downloads.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-8.2/wgs-8.2-Linux_amd64.tar.bz2
POA_PATH=http://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz
PROCESS=8

#
# Download raw data, rename to ENA accession
#
ERX708228.tar:
	wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_20_LomanLabz_PC_Ecoli_K12_R7.3.tar -O ERX708228.tar

ERX708229.tar:
	wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_32_LomanLabz_K12_His_tag.tar -O ERX708229.tar

ERX708230.tar:
	wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_33_LomanLabz_PC_K12_0.4SPRI_Histag.tar -O ERX708230.tar

ERX708231.tar:
	wget ftp.sra.ebi.ac.uk/vol1/ERA411/ERA411499/oxfordnanopore_native/flowcell_39.tar -O ERX708231.tar

#
# Untar into directories named x.fast5
# The find command moves the deeply nested fast5 files into the top directory
#
%.fast5/: %.tar
	mkdir -p $@
	cd $@; tar -xf ../$<
	find $@ -name "*.fast5" -exec mv {} $@ \;

#
# Check that poretools is installed and can be run
#
poretools.version:
	poretools -v 2> poretools.version

#
# Install nanocorrect & dependencies
#
nanocorrect.version:
	git clone https://github.com/jts/nanocorrect.git
	ln -s nanocorrect/poa-blosum80.mat
	-cd nanocorrect; git log | head -1 > ../$@

#
# Install Python dependencies
#
pythonlibs.version:
	pip install pysam > $@
	pip install poretools >> $@

#
# Install poa
#
poa.version:
	wget $(POA_PATH)
	tar -xzf poaV2.tar.gz
	cd poaV2; make poa
	ln -s poaV2/poa
	echo $(POA_PATH) > $@

#
# Install DALIGNER
#
daligner.version:
	git clone https://github.com/thegenemyers/DALIGNER.git
	cd DALIGNER; git checkout 549da77b91395dd; make
	echo "549da77b91395dd" > $@

#
# Install DAZZ_DB
#
dazz_db.version:
	git clone https://github.com/thegenemyers/DAZZ_DB
	cd DAZZ_DB; git checkout 8cb2f29c4011a2c2; make
	echo "8cb2f29c4011a2c2" > $@

#
# Install Celera Assembler
#
ca.version:
	wget $(CA_PATH)
	tar -xjf wgs-8.2-Linux_amd64.tar.bz2
	echo $(CA_PATH) > $@

#
# Install samtools
#
samtools.version:
	git clone --recursive https://github.com/samtools/htslib.git
	cd htslib; make
	git clone --recursive https://github.com/samtools/samtools.git
	cd samtools; make
	-cd samtools; git log | head -1 > ../$@

#
test_path: FORCE
	-LAshow 2> test
	-fasta2DB 2>> test
	-poa 2>> test
FORCE:

#
# Export reads with poretools
#
raw.reads.fasta: ERX708228.fast5 ERX708229.fast5 ERX708230.fast5 ERX708231.fast5
	poretools fasta --type 2D ERX708228.fast5/ > $@
	poretools fasta --type 2D ERX708229.fast5/ >> $@
	poretools fasta --type 2D ERX708230.fast5/ >> $@
	poretools fasta --type 2D ERX708231.fast5/ >> $@

#
# Run error correction
#
%.corrected.fasta: %.fasta samtools.version
	make -f nanocorrect/nanocorrect-overlap.make INPUT=$< NAME=$*
	samtools faidx $<
	python nanocorrect/makerange.py $< | head -1 | parallel --progress -P $(PROCESS) 'python nanocorrect/nanocorrect.py $* {} > $*.{}.corrected.fasta'
	cat $*.*.corrected.fasta > $@
	#rm $*.*.corrected.fasta

#
# Run CA on twice corrected reads
# We need both files here to make the error correction rule evalulated twice
#
assembly.input.fasta: raw.reads.corrected.fasta raw.reads.corrected.corrected.fasta
	echo "Running CA"
