Sample Matching for large-scale Proteomics *(SMAP)*
-------------------------------------------------------

A pipeline to validate and correct sample identity based on a
combination of concordance and specificity scores. SMAP first detects
variant peptides from multiplexed isobaric labeling-based quantitative
proteomics data using the proteogenomics approach, and then infers
allelic information for each sample based on its expression level of the
variant peptides.

**Highlights**

use mutant peptide to get the peptide's chromosome and start position

use the peptide's chromosome and start position to extract the
corresponding reference peptides in the file of

extract quantification data from both reference and mutation peptides
from the file of quantity results of peptides.

Sex information is not necessary but having both genetics-based and
reported sexes will help identify true IDs.

Software requirement
--------------------

Linux operating system

Get started
-----------

Users can use the following script to detect and correct samples.

**perl pipeline.pl &gt; outcome.txt**


Input file examples
-------------------

SNP and protein database 

The variant data were re-annotated using the genome annotation tool
ANNOVAR (Wang et al., 2010) based on the human reference genome (hg19)
genome assembly. A customized protein database was constructed by
appending human UniprotKB database (UniProt, 2021) with <span
id="_Hlk65489796" class="anchor"></span>SNPs with the amino acid
sequences of nonsynonymous variants. The identified variant peptides are
furthered filtered using the target-decoy strategy for false discovery
rate (FDR) evaluation (Junmin Peng, 2003). Using minor allele frequency
(MAF) &lt; 1% as a filter cutoff value.

**mutation\_events.txt**

  mutation                        gene      aa\_change   peptide\_counts   peptide\_sequences     scan\_counts   best\_scans
  ------------------------------- --------- ------------ ----------------- ---------------------- -------------- ------------------------------------
  chr3:133524717-133524717.G\_C   SRPRB     p.V9L        1                 LADGGGAGGTFQPYLDTLR    1              psy\_dwp\_b01\_f12.22198.1.2.spout
  chr2:234183368-234183368.A\_G   ATG16L1   p.T317A      1                 SVSSFPVPQDNVDAHPGSGK   1              psy\_dwp\_b01\_f01.14515.1.3.spout
  chr1:223954080-223954080.A\_C   CAPN2     p.K490Q      1                 RQDIQSDGFSIETCK        2              psy\_dwp\_b01\_f34.13389.1.3.spout

**reference\_peptides.bed (empty in original folder)**

**mutation\_peptides.bed**

**raw\_quan\_HH\_SNP\_b05\_psm\_nonzero.txt**

Database=/home/groupdirs/wanglab/ling/project/genotype\_db\_search/gnm\_gnm\_SNP\_b05/intermediate/sum\_accepted\_PSMs/misc/idsum.db

Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;group;subgroup;unique;tryptic;pos;precursor\_peak\_intensity\_percentage;RT;K\_y1;R\_y1;charge;endingAA;nLabels;mz126
(S1);mz127N (S2);mz127C (S3);mz128N (S4);mz128C (S5);mz129N (S6);mz129C
(S7);mz130N (S8);mz130C (S9);mz131N (S10);mz131C (S11);sig126
(S1);sig127N (S2);sig127C (S3);sig128N (S4);sig128C (S5);sig129N
(S6);sig129C (S7);sig130N (S8);sig130C (S9);sig131N (S10);sig131C (S11)

R.SPTSSVATPSSTISTPTK.R;sp|Q8WXI2|CNKR2\_HUMAN;/home/groupdirs/wanglab/ling/project/genotype\_db\_search/gnm\_gnm\_SNP\_b05/intermediate/psy\_dwp\_b05\_f05/psy\_dwp\_b05\_f05.1/psy\_dwp\_b05\_f05.20837.1.2.spout;2207.20996;2207.2170532033;-3.15571477306367;47.52;0.743897306397306;0/0;0;103430;1;1;2;AA325toAA342;100.00%;3529.38;41753.5546875;0;2;K;2;126.127690687793;127.124745028234;127.130993545698;128.128086033376;128.134342180288;129.131404150176;129.137660297088;130.134722266975;130.140978413887;131.13805564267;131.144327048477;97487.88;82986.48;80333.15;81907.46;79867.44;84530.66;117866.24;83086.55;84191.40;84290.17;94521.85;

R.SPTSSVATPSSTISTPTK.R;sp|Q8WXI2|CNKR2\_HUMAN;/home/groupdirs/wanglab/ling/project/genotype\_db\_search/gnm\_gnm\_SNP\_b05/intermediate/psy\_dwp\_b05\_f05/psy\_dwp\_b05\_f05.1/psy\_dwp\_b05\_f05.20962.1.2.spout;2207.21314;2207.2170532033;-1.58307875043492;62.26;0.819627369097334;0/0;0;103430;1;1;2;AA325toAA342;100.00%;3547.2;43743.43359375;0;2;K;2;126.127705946688;127.124737398787;127.131016434041;128.128055515586;128.134357439183;129.131419409071;129.137675555982;130.13473752587;130.141024190572;131.13808616046;131.144357566267;96007.54;96167.68;64648.51;88606.83;76993.10;93138.50;106253.97;96484.06;88902.38;86105.93;100864.16;

**raw\_quan\_HH\_SNP\_b05\_psm\_zero.txt**

Database=/home/groupdirs/wanglab/ling/project/genotype\_db\_search/gnm\_gnm\_SNP\_b05/intermediate/sum\_accepted\_PSMs/misc/idsum.db

Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;group;subgroup;unique;tryptic;pos;precursor\_peak\_intensity\_percentage;mz126
(S1);mz127N (S2);mz127C (S3);mz128N (S4);mz128C (S5);mz129N (S6);mz129C
(S7);mz130N (S8);mz130C (S9);mz131N (S10);mz131C (S11);sig126
(S1);sig127N (S2);sig127C (S3);sig128N (S4);sig128C (S5);sig129N
(S6);sig129C (S7);sig130N (S8);sig130C (S9);sig131N (S10);sig131C (S11)

K.ADTSQEICSPR.L;cu|84681\_76308|cu;/home/groupdirs/wanglab/ling/project/genotype\_db\_search/gnm\_gnm\_SNP\_b05/intermediate/psy\_dwp\_b05\_f12/psy\_dwp\_b05\_f12.1/psy\_dwp\_b05\_f12.9058.1.2.spout;1492.72371;1492.7262584131;-1.54132392774583;19.36;0.496384297520661;0/0;2;179267;1;0;2;AA12toAA22;100.00%;126.127799060325;127.124838214515;127.131109620779;128.128293734481;128.134473587369;129.131459335423;129.137761259479;130.134792784221;130.14107944938;131.138172009707;131.14444341597;8806.566;8118.488;8576.511;7348.662;8772.410;7472.893;8777.004;3674.616;3787.724;9678.674;10963.492;

K.ADTSQEICSPR.L;cu|84681\_76309|cu;/home/groupdirs/wanglab/ling/project/genotype\_db\_search/gnm\_gnm\_SNP\_b05/intermediate/psy\_dwp\_b05\_f12/psy\_dwp\_b05\_f12.1/psy\_dwp\_b05\_f12.9058.1.2.spout;1492.72371;1492.7262584131;-1.54132392774583;19.36;0.496384297520661;0/0;2;179267;2;0;2;AA12toAA22;100.00%;126.127799060325;127.124838214515;127.131109620779;128.128293734481;128.134473587369;129.131459335423;129.137761259479;130.134792784221;130.14107944938;131.138172009707;131.14444341597;8806.566;8118.488;8576.511;7348.662;8772.410;7472.893;8777.004;3674.616;3787.724;9678.674;10963.492;

Output file
-----------

**MS\_SNPmatrix\_head.txt**


This is a temporary file to save the genotypes in each sample, the columns are including: \#CHROM, POS, ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, and sample ids.

**Mutant\_Reference\_quan\_results.txt**

*This is a temporary file to save the quantity results.*

chr1:22839498-22839498.A\_G ZBTB40 p.N848S 1 FAASSTLK 1 psy\_dwp\_b01\_f19.8838.1.2.spout cu|5378\_41890|cu /home/xwang4/2018/Psycho\_genotype\_2/Database\_search/Batch1/gnm\_SNP\_MS/intermediate/psy\_dwp\_b01\_f19/psy\_dwp\_b01\_f19.1/psy\_dwp\_b01\_f19.8838.1.2.spout 167018.90 103700.16 163867.60 105834.09 154465.75 118336.61 114133.33 80552.52 93904.62 296129.15 111908.22 0.401093476598089 0.1073754608744 0.386475472781999 0.117274168354891 0.342862906800241 0.175269879671094 0.155772033360017 0 0.0619366765312176 1 0.145450367231365 nonzero

chr19:622336-622336.T\_G POLRMT p.E555A 1 QYWEALGAPEALR 1 psy\_dwp\_b01\_f12.20429.1.2.spout cu|5378\_56811|cu /home/xwang4/2018/Psycho\_genotype\_2/Database\_search/Batch1/gnm\_SNP\_MS/intermediate/psy\_dwp\_b01\_f12/psy\_dwp\_b01\_f12.1/psy\_dwp\_b01\_f12.20429.1.2.spout 29395.04 74782.12 90367.44 45756.24 46073.51 49550.12 51208.19 42921.40 58650.97 60005.98 39541.38 0 0.74438729654729 1 0.268337805302071 0.273541307214412 0.330560712715917 0.357754492196469 0.221843981867205 0.479822509856919 0.502045843693212 0.16640873575585 nonzero

**Inferred\_results.txt**

*SMAP infers sample allelic information based on the relative expression*
*level of each sample, then infers allelic information for each sample*
*based on the scaled intensities of the identified variant peptides*

*There are two parts in this file, first is the counting numbers, the*
*second part is the scores.*

  Part 1                                 
  -------- ------------ ---------------- -------------------
  id       sample\_id   matched Number   unmatched Number
  24       2015-1426    71               81
  24       2016-867     68               84
  24       2015-737     74               78
  Part2                                  
  Id       Sample id    score            $$\text{Cscore}$$
  21       2015-767     2.818182         0.495037
  22       2015-921     3.666667         0.509091
  23       2016-926     3.666667         0.43459

**output.txt**

This file saved all the results of counting steps, you could print to
the screen if no saving in the file.

Data preparation
----------------

**Identify variant peptides in the proteomics data** 


All the steps in the analysis of a MS-based proteomic data are used the
JUMPg program.JUMPg program (Li et al., 2016) was used to identify
variant peptides, which takes MS-based proteomic data as input. JUMPg
constructs a customized database using genomic variant files (e.g., VCF
files), in which JUMP performs preprocessing, tag generation, MS/MS
pattern matching, and scoring as previously reported (Wang et al.,
2014). <span id="_Hlk65489866" class="anchor"></span>The variant data
were re-annotated using the genome annotation tool ANNOVAR (Wang et al.,
2010) based on the human reference genome (hg19) genome assembly. A
customized protein database was constructed by appending human UniprotKB
database (UniProt, 2021) with SNPs with the amino acid sequences of
nonsynonymous variants. The identified variant peptides are furthered
filtered using the target-decoy strategy for false discovery rate (FDR)
evaluation (Junmin Peng, 2003). Using minor allele frequency (MAF) &lt;
1% as a filter cutoff value.

Contact information
-------------------

For any questions please contact [*xusheng.wang@und.edu*](mailto:xusheng.wang@und.edu)

[Packages](https://github.com/users/jart/packages?repo_name=cosmopolitan)

No packages published 
----------------------

Junmin Peng, J.E.E., Carson C Thoreen, Larry J Licklider, Steven P Gygi.
(2003). Evaluation of multidimensional chromatography coupled with
tandem mass spectrometry (LC/LC-MS/MS) for large-scale protein analysis
the yeast proteome.pdf&gt;. J Proteome Res *2*, 43-50.

Li, Y., Wang, X., Cho, J.H., Shaw, T.I., Wu, Z., Bai, B., Wang, H.,
Zhou, S., Beach, T.G., Wu, G.*, et al.* (2016). JUMPg: An Integrative
Proteogenomics Pipeline Identifying Unannotated Proteins in Human Brain
and Cancer Cells. J Proteome Res *15*, 2309-2320.

UniProt, C. (2021). UniProt: the universal protein knowledgebase in
2021. Nucleic acids research *49*, D480-D489.

Wang, K., Li, M., and Hakonarson, H. (2010). ANNOVAR: functional
annotation of genetic variants from high-throughput sequencing data.
Nucleic acids research *38*, e164.

Wang, X., Li, Y., Wu, Z., Wang, H., Tan, H., and Peng, J. (2014). JUMP:
a tag-based database search tool for peptide identification with high
sensitivity and accuracy. Mol Cell Proteomics *13*, 3663-3673.
