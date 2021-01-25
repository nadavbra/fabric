What is FABRIC?
===============

FABRIC (Functional Alteration Bias Recovery In Coding-regions) is a framework for detecting genes showing functional alteration bias in some evolutionary context. For example, in the context of cancer genomics it can be used to detect alteration promoting genes (namely genes affected by mutations that are significantly more harmful than expected at random). Cancer alteration promoting genes are strong candidates for cancer driver genes. Likewise, in the context of population genetic variation it can be used to detect alteration rejecting genes (namely genes affected by variants that are significantly less harmful than expected at random). Such genes are likely the product of negative selection, and are expected to harbor important functions.

The framework relies on a machine-learning prediction model for assessing the functional impact of genetic variants. All the examples shown here are based specifically on `FIRM <https://github.com/nadavbra/firm>`_, a specific prediction model designed exactly for that task, but FABRIC can use any other model. Note that while in FIRM a score of 0 indicates a harmless mutation and a score of 1 indicates a harmful mutation, in FABRIC it is the other way around.

Importantly, FABRIC is not sensitive (in terms of false discoveries) to the accuracy of the underlying prediction model, as it relies on precise statistical calculations. Specifically, it compares each gene's observed effect scores (calculated by applying the prediction model over the observed mutations in the gene) to its gene-specific background effect score distribution expected at random (calculated by applying the same prediction model over all possible mutations in the gene).

To see how the results of FABRIC look like, you can explore `The FABRIC Cancer Portal <http://fabric-cancer.huji.ac.il/>`_, a catalouge of FABRIC's summary statistics across the entire human coding genome (~18,000 protein-coding genes) and across all 33 `TCGA <https://portal.gdc.cancer.gov/>`_ cancer types and pan-cancer. 

For more details on FABRIC you can refer to our paper `Quantifying gene selection in cancer through protein functional alteration bias <https://doi.org/10.1093/nar/gkz546>`_, published in Nucleic Acids Research (2019).

Or, if you are more a video person, you can watch this talk on YouTube (originally given at ISMB 2020):

.. image:: https://img.youtube.com/vi/GUPmRZiLMUw/0.jpg
   :target: https://www.youtube.com/watch?v=GUPmRZiLMUw
   
   
Installation
============

FABRIC is a pure Python3 package.

To install FABRIC, simply run:

.. code-block:: sh

   pip install fabric-genetics
   
FABRIC requires (and will atomatically isntall) the following common packages:

* numpy
* scipy
* pandas
* statsmodels


Usage
=====

If you want to analyze your own list of variants/mutations with FABRIC (to quantify the bias of affected genes towards more or less functional damage), simple run:

.. code-block:: sh

   fabric --input-variants-file=your_variant_list.csv --possible-variant-effects-file=all_cds_snp_effects_GRCh38(|hg19).csv --genes-file=genes_GRCh38(|hg19).csv --output-file=fabric_output.csv
   
The files :code:`all_cds_snp_effects_GRCh38.csv` and :code:`genes_GRCh38.csv` (or the equivalent files for version *hg19* of the human reference genome) can be obtained from `ftp://ftp.cs.huji.ac.il/users/nadavb/firm_data/ <ftp://ftp.cs.huji.ac.il/users/nadavb/firm_data/>`_ (you can also create them from scratch using `FIRM <https://github.com/nadavbra/firm>`_'s :code:`firm_list_all_possible_cds_snps` command, but that is a long process).
Make sure that you are using the files matching the version of the reference genome that is relevant for your list of variants (provided by the :code:`--input-variants-file` argument). If you wish to use the effect scores of a tool other than FIRM, the arguments :code:`--possible-variant-effects-file` and :code:`--genes-file` can accept any CSV files, provided that they are at the right format (more on that in the "Using other effect scores" section).

The CSV file provided by :code:`--input-variants-file` can list any set of independent variants (one in each row). Each variant is recognized by its locus (chromosome & position) and the reference & alternative DNA sequences. By default, FABRIC expects the following column names in the input CSV: *chrom*, *pos*, *ref* and *alt*. However, you can specify any othr column names by using the arguments :code:`--chrom-column`, :code:`--pos-column`, :code:`--ref-column` and :code:`--alt-column`. Other columns in the input CSV file are simply ignored. Use the flag :code:`--input-variants-tab-delimiter:code:` if the input CSV file is separated by tabs rather than commas. 

Important: Input variants have to be independent! FABRIC ignores complex (non single-nucleotide) variants or variants that are not affecting the CDS of a recognized protein-coding genes.

By default, FABRIC analyzes all the input variants jointly (per gene), which provides maximal statistical power. If you want FABRIC to analyze groups of variants indepndently, you can include another column in the input CSV that specifies a unique label for each group of variants that should be analyzed as an indpeendent group. You can then ask FABRIC to split variants according to this column by setting the --analyze-by-project-column argument to be the name of that column. Even after providing this column, FABRIC will still also include a combined analysis of all variants together (unless you also add the --skip-combined-analysis flag). When analyzing different sub-projects in this way, you can also add the flag --analyze-diff to also analyze and output the differences between projects (per gene).

The output CSV file(s) produced by FABRIC will contain the summary statistics of each of the analyzed protein-coding genes, including the p-value, q-value and z-value of each gene. Negative z-values indicate potential positive selection (i.e. more functional damage than would be expected at random), while positive z-values indicate potential negative selection (i.e. less functional damage than would be expected at random).

For more options, run:

.. code-block:: sh

   fabric --help


Example 1: Analyzing cancer somatic mutations from TCGA 
-----------

In this example, we will analyze ~3M somatic mutations from 33 cancer types obtained from the TCGA (which is exactly the same dataset analyzed in  `The FABRIC Cancer Portal <http://fabric-cancer.huji.ac.il/>`_). You can download the relevant dataset (gdc_combined.csv) from ftp://ftp.cs.huji.ac.il/users/nadavb/fabric_examples/gdc_combined.maf. You can also generate this file yourself through the Jupyter Notebook provided in this GitHub repository (go to the "Combine GDC's downloaded tar file into a single MAF file" section in that notebook).

To analyze these mutations through a combined (pan-cancer) analysis, simply run:

.. code-block:: sh

   fabric --input-variants-file=gdc_combined.maf --possible-variant-effects-file=all_cds_snp_effects_GRCh38.csv --genes-file=genes_GRCh38.csv --output-file=gdc_pan_cancer_fabric_results.csv --input-variants-tab-delimiter --chrom-column=Chromosome --pos-column=Start_Position --ref-column=Tumor_Seq_Allele1 --alt-column=Tumor_Seq_Allele2
   
Recall that the files all_cds_snp_effects_GRCh38.csv and genes_GRCh38.csv can be taken from ftp://ftp.cs.huji.ac.il/users/nadavb/firm_data/.
   
If you want to also include a separate analysis for each of the 33 cancer types, run instead:

.. code-block:: sh

   fabric --input-variants-file=gdc_combined.maf --possible-variant-effects-file=all_cds_snp_effects_GRCh38.csv --genes-file=genes_GRCh38.csv --output-dir=gdc_fabric_results --analyze-by-project-column=tcga_project --analyze-diff --input-variants-tab-delimiter --chrom-column=Chromosome --pos-column=Start_Position --ref-column=Tumor_Seq_Allele1 --alt-column=Tumor_Seq_Allele2
   
Since this is going to analyze 33 TCGA projects independently, it's going to take a long time to run, so it's recommended to run it with nohup or a similar tool.


Example 2: Analyzing genetic variants in the healthy human population from ExAC
-----------

In this example, we will analyze ~9M variants sequenced from the exomes of ~60K individuals obtained from ExAC (http://exac.broadinstitute.org/). The file is available at:
ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz.

Since this is a VCF file, we will first need to convert it into CSV using the vcf_to_csv tool installed by FABRIC. Simply run:

.. code-block:: sh

   vcf_to_csv --vcf-file=ExAC.r1.sites.vep.vcf.gz --output-csv-file=exac_variants.csv --only-pass
   
The --only-pass flag is used to only retrieve variants passing the quality-control filter in the VCF file (i.e. with "PASS" in the FILTER field).

After you have convereted the data into CSV format, you can run FABRIC over this dataset:

.. code-block:: sh

   fabric --input-variants-file=exac_variants.csv --possible-variant-effects-file=all_cds_snp_effects_hg19.csv --genes-file=genes_hg19.csv --output-file=exac_fabric_results.csv
   
Recall that the files all_cds_snp_effects_hg19.csv and genes_hg19.csv can be taken from ftp://ftp.cs.huji.ac.il/users/nadavb/firm_data/.


Using other effect scores
=====
    

Cite us
=======

If you use FABRIC as part of work contributing to a scientific publication, we ask that you cite our paper: Nadav Brandes, Nathan Linial, Michal Linial, Quantifying gene selection in cancer through protein functional alteration bias, Nucleic Acids Research, gkz546, https://doi.org/10.1093/nar/gkz546
