# Bioinformatics Pipeline – Bolnick Lab 2015

> As practiced by Yoel Stuart, under the auspices of Jesse Weber’s expertise

> Small stylistic changes by Lucas Nell (June 2016)


# Notes about the data
We built 18 48-individual libraries and sequenced those libraries across 4 sequencing
runs, using 12 lanes, such that every library was sequenced in every lane. Thus, we
received 12 R1 and 12 R2 fastq.gz files for each of the 18 libraries, total, from the
four different sequencing runs.

# 1) CONCATENATION
We first wished to concatenate the fastq.gz files for each library. This should be done
by read, such that one is left with one *R1.fastq.gz and one *R2.fastq.gz file per
library.
```bash
cat SA14095/Lib02/Lib2Final_GCCAAT_L001_R1_001.fastq.gz \
    SA14095/Lib02/Lib2Final_GCCAAT_L002_R1_001.fastq.gz \
    > pool_demultiplex/Lib02/Lib02_R1andR2gzfiles/Lib02cat_R1_001.fastq.gz
```

# 2) PROCESS RADTAGS
We used `process_radtags` from STACKS to demultiplex the concatenated files.
For this, we created a textfile that provides the flex barcodes and the library index,
in 2 columns, that looked something like this, except for all 48 barcodes.
```
GCATG    CGATGT
AACCA    CGATGT
CGATC    CGATGT
TCGAT    CGATGT
TGCAT    CGATGT
CAACC    CGATGT
```

The `process_radtags` call for a single library calls the `*R1.fastq.gz` and the
`R2.fastq.gz` files.
```bash
process_radtags -P -p pool_demultiplex/Lib01/Lib01_R1andR2gzfiles \
    -b pool_demultiplex/Lib01/flexbarcodes_withindex_tw_Lib01.txt \
    –o pool_demultiplex/Lib01/Lib01_demulti_individuals \
    --renz_1 nlaIII --renz_2 mluCI \
    -r -i gzfastq --inline_index --disable_rad_check
```
Then remove the remainder files, which should be empty.
```bash
rm *rem.*
```

# 3) HOUSEKEEPING I
`process_radtags` generates demultiplexed `.fq` files for each individual, in a file
designated by the `–o` flag. Because the headers from R2 contain an "_2", they won’t
match their mate pair and will trip an error during the downstream `BWA` steps. Some
housekeeping is necessary to change header names and rename the changed files.

```bash
for f in `ls *2.fq`
do
    ff=${f%2.fq}fixed_2.fq
    sed 's/_2$/_1/g' $f > $ff
done
```

# 4) DOWNLOAD AND INDEX THE GENOME FOR BWA
The latest genome assembly can be downloaded from Ensemble using `wget`. Choose the
"unmasked, DNA, top-level" genome. The genome then needs to be indexed in `BWA` to
facilitate genome mapping; this step adds hash tags and builds a Burrow-Wheeler’s
table. Make sure the genome is readable.
```bash
FTP_DIR=ftp://ftp.ensembl.org/pub/release-78/fasta/gasterosteus_aculeatus/dna
wget $FTP_DIR/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa.gz
gunzip Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa.gz
bwa index Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa
```


# 5) BWA ALN
This step creates `.sai` files, which area a SAM Alignment Index, that are used in
concert with the `.fastq` files (from `process_radtags`) to create a single `.sam`
file for each individual using `sampe`.
The following code creates an output on the terminal window that is copied and
pasted into a shell script for running on a cluster.
```bash
FOLDER=/work/02631/yestuart/sticklegenome/ensemble_stickle_v078
for R1 in ./*1.fq
do
    name=`echo $R1 | sed 's/.1.fq\+//'`
    R2=${R1%1.fq}fixed_2.fq
    echo bwa aln $FOLDER/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa \
        $R1 '>' $name.1.sai '2>' $name.1.log
    echo bwa aln $FOLDER/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa \
        $R2 '>' $name.2.sai '2>' $name.2.log
done
```


# 6) BWA SAMPE
This step creates a `.sam` file, with sequences aligned to the genome.
The following code creates an output on the terminal window that is copied and pasted
into a shell script for running on the cluster.
```bash
for R1 in `ls *1.fq`
do
    name=${R1%1.fq}1
    samname=${R1%.1.fq}_Lib05_45
    R2=${R1%1.fq}fixed_2
    R2sai=${R1%1.fq}2.sai

    echo "bwa sampe -P -r '@RG\tID:"$samname'\tPL:Illumina\tLB:'$samname'\tSM:'$samname"'" \
    yestuart/sticklegenome/ensemble_stickle_v078/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa \
    $name.sai $R2sai $R1 $R2.fq '>' $samname.sam
done
```

# 7) HOUSEKEEPING II
In preparation for making `.bam` files, download and unpack STAMPY into the cluster.
Make sure STAMPY is in your path.

# 8) DOWNLOAD AND INDEX THE GENOME FOR STAMPY
Download the genome again, this time indexing and creating a hash table with STAMPY
requirements.
```bash
#index for stampy
/work/02631/yestuart/Stampy/stampy-1.0.23/stampy.py \
    -G Gasterosteus_aculeatus.BROADS1.dna.toplevel \
    Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa
#create hash table for stampy
    /work/02631/yestuart/Stampy/stampy-1.0.23/stampy.py \
    -H Gasterosteus_aculeatus.BROADS1.dna.toplevel \
    -g Gasterosteus_aculeatus.BROADS1.dna.toplevel
```

# 9) MAKE SORTED.STAMPY.BAM FILES WITH STAMPY
This step relies on a shell script developed by J. Weber that takes `.sam` files,
turns them to .bam files using Samtools, and runs them through Stampy for another
alignment: `"Sam_to_bam_to_stampy_new.sh"`. This creates a `sorted.stampy.bam` file.

# 10) SORT THE SORTED.STAMPY.BAM FILES

This script creates a shell file which can then be run on the cluster.

```bash
for R1 in `ls -S *sorted.stampy.bam`
do
    echo samtools sort $R1 $R1.sorted >> sort_stampy_bams.sh
done
```

# 11) INDEX THE SORTED.STAMPY.BAM FILES
This script creates a shell file which can then be run on the cluster.
```bash
for R1 in `ls -S *sorted.stampy.bam.sorted.bam`
do
    echo samtools index $R1 >> z_index_stampy_bams.sh
done
```

# 12) CREATE A BAMLIST
The bamlist is used by `mpileup`. We wanted population-level SNP calling, so we
created bamlists by each watershed. The following text creates the bamlist.
```bash
for SAM in `ls -S *bam.sorted.bam`; do
    echo $SAM >> z_bam_list.txt
done
chmod a+rx z_bam_list.txt
```

# 13) INDEX THE GENOME FOR SAMTOOLS
```bash
samtools faidx Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa
```

# 14) RUN MPILEUP
Calling the bamlist, and the `samtools`-appropriate genome, run `mpileup` and output
accordingly. This generates the genotype probabilities. Here is where you want to
be smart about your filtering.
```bash
ref=sticklegenome/ensemble_stickle_v078_forsamtools/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa
samtools mpileup -C 50 -E -S -D -u -I \
    -f $ref -b z_bam_list.txt > Thor.mpileup.output.bcf
```

# 15) GENERATE VCF FILES FROM BCF FILES OUTPUTTED FROM MPILEUP
Use `bcftools` to generate a `.vcf` file, which will have the SNP calls.
```bash
bcftools view Boot.mpileup.output.bcf -v -c -g > BootLS.vcf
```

# 16) FILTER FOR MEAN MIN AND MAX DEPTH OF COVERAGE
Filtered max-mean of 75, and min-mean of 1.
```bash
vcftools --vcf h_bcftools.viewSNPS.ThorPool3b.vcf \
    --out i_vcftools.meandepthfilter.1.75.ThorPool3b \
    --min-meanDP 1 --max-meanDP 75 --recode
```

17) FILTER FOR INDIVIDUAL X SITE MIN AND MAX DEPTH OF COVERAGE
Filtered each site at each individual for minimum depths of 4, 6, and 8. Max depth of
100.
```bash
vcftools --vcf i_vcftools.meandepthfilter.1.75.ThorPool3b.recode.vcf \
    --out j_filtered.mean1.75.ThorPool3b.maxDP100.minDP8 \
    --minDP 8 --maxDP 100 --recode
vcftools --vcf i_vcftools.meandepthfilter.1.75.ThorPool3b.recode.vcf \
    --out j_filtered.mean1.75.ThorPool3b.maxDP100.minDP6 \
    --minDP 6 --maxDP 100 --recode
vcftools --vcf i_vcftools.meandepthfilter.1.75.ThorPool3b.recode.vcf \
    --out j_filtered.mean1.75.ThorPool3b.maxDP100.minDP4 \
    --minDP 4 --maxDP 100 --recode
```

# 18) FILTER FOR COMPLETENESS
Filtered for completeness of 0.5, 0.8, and 0.9.
```bash
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP4.recode.vcf \
    --out k_completeness50percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP4 \
    --max-missing 0.5 --recode
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP4.recode.vcf \
    --out k_completeness80percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP4 \
    --max-missing 0.8 --recode
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP4.recode.vcf \
    --out k_completeness90percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP4 \
    --max-missing 0.9 --recode
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP6.recode.vcf \
    --out k_completeness50percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP6 \
    --max-missing 0.5 --recode
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP6.recode.vcf \
    --out k_completeness80percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP6 \
    --max-missing 0.8 --recode
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP6.recode.vcf \
    --out k_completeness90percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP6 \
    --max-missing 0.9 --recode
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP8.recode.vcf \
    --out k_completeness50percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP8 \
    --max-missing 0.5 --recode
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP8.recode.vcf \
    --out k_completeness80percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP8 \
    --max-missing 0.8 --recode
vcftools --vcf j_filtered.mean1.75.ThorPool3b.maxDP100.minDP8.recode.vcf \
    --out k_completeness90percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP8 \
    --max-missing 0.9 --recode
```

# 19) GENERATE 012 FILES
.012 has 0, 1, and 2 for homozygote consensus, heterozygote, and homozygote alternate
allele, relative to the genome.
```bash
vcftools --vcf  k_completeness50percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP4.recode.vcf \
    --out l_ThorPool3b.50pcomp.DP4 --012
vcftools --vcf  k_completeness80percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP4.recode.vcf \
    --out l_ThorPool3b.80pcomp.DP4 --012
vcftools --vcf  k_completeness90percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP4.recode.vcf \
    --out l_ThorPool3b.90pcomp.DP4 --012
vcftools --vcf  k_completeness50percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP6.recode.vcf \
    --out l_ThorPool3b.50pcomp.DP6 --012
vcftools --vcf  k_completeness80percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP6.recode.vcf \
    --out l_ThorPool3b.80pcomp.DP6 --012
vcftools --vcf  k_completeness90percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP6.recode.vcf \
    --out l_ThorPool3b.90pcomp.DP6 --012
vcftools --vcf  k_completeness50percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP8.recode.vcf \
    --out l_ThorPool3b.50pcomp.DP8 --012
vcftools --vcf  k_completeness80percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP8.recode.vcf \
    --out l_ThorPool3b.80pcomp.DP8 --012
vcftools --vcf  k_completeness90percent.filtered.mean1.75.ThorPool3b.maxDP100.minDP8.recode.vcf \
    --out l_ThorPool3b.90pcomp.DP8 --012
```
