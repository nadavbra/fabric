What is FABRIC?
===============

FABRIC (Functional Alteration Bias Recovery In Coding-regions) is a framework for detecting genes showing functional alteration bias in some evolutionary context. For example, in the context of cancer genomics it can be used to detect alteration promoting genes (namely genes affected by mutations that are significantly more harmful than expected at random). Cancer alteration promoting genes are strong candidates for cancer driver genes. Likewise, in the context of population genetic variation it can be used to detect alteration rejecting genes (namely genes affected by variants that are significantly less harmful than expected at random). Such genes are likely the product of negative selection, and are expected to harbor important functions.

The framework relies on a machine-learning prediction model for assessing the functional impact of genetic variants. Currently, it is coupled with FIRM (https://github.com/nadavbra/firm), a specific prediction model designed exactly for that task. However, it should be relatively easy to expand the module to other models as well. Note that while in FIRM a score of 0 indicates a harmless mutation and a score of 1 indicates a harmful mutation, in FABRIC it is the other way around.

Importantly, FABRIC is not sensitive (in terms of false discoveries) to the accuracy of the underlying prediction model, as it relies on precise statistical calculations. Specifically, it compares each gene's observed effect scores (calculated by applying the prediction model over the observed mutations in the gene) to its gene-specific background effect score distribution expected at random (calculated by applying the same prediction model over all possible mutations in the gene).

To see how the results of FABRIC look like, you can explore `The FABRIC Cancer Portal <http://fabric-cancer.huji.ac.il/>`_, a catalouge of FABRIC's summary statistics across the entire human coding genome (~18,000 protein-coding genes) and across all 33 `TCGA <https://portal.gdc.cancer.gov/>`_ cancer types and pan-cancer. 

For more details on FABRIC you can refer to our paper `Quantifying gene selection in cancer through protein functional alteration bias <https://doi.org/10.1093/nar/gkz546>`_, published in Nucleic Acids Research (2019).

Or, if you are more a video person, you can watch this talk on YouTube (originally given at ISMB 2020):

.. image:: https://img.youtube.com/vi/GUPmRZiLMUw/0.jpg
   :target: https://www.youtube.com/watch?v=GUPmRZiLMUw


Usage
=====

Preparation
-----------

After installtion (see below), you will need to prepare the background effect score distributions for all the protein-coding genes in the version of the human reference genome you plan to work with (GRCh38 or hg19). This is done by running the following command:

        set_gene_bg_scores --reference-genome=<GRCh38 or hg19> --gene-dataset-dir=<DIR PATH> --gene-bg-scores-dir=<DIR PATH>

Where:

* --gene-dataset-dir should be given a directory where you wish to store a CSV file with the details of all the genes involved in the analysis.
* --gene-bg-scores-dir should be given the directory in which you want to create the background effect scores.

For example, if you work with GRCh38 and your working directory is ~/fabric_data, you will run the command:

        set_gene_bg_scores --reference-genome=GRCh38 --gene-dataset-dir=~/fabric_data/ --gene-bg-scores-dir=~/fabric_data/gene_bg_scores/


This will create the directory ~/fabric_data/gene_bg_scores/GRCh38/ and fill it with JSON files (one for each gene, identified by its UniProt ID) containing the background effect scores of all analyzed genes. 
Note that this should be a very lengthy process, so it is recommended to run it in the background, and to use a machine with many CPUs. 

**Important note**: You can actually spare this very long process by using pre-calculated background effect score distributions that are available via FTP at: ftp://ftp.cs.huji.ac.il/users/nadavb/fabric_data/. For a particular version of the human reference genome, you will need to get the file genes_<version>.csv and put it in your gene dataset directory (e.g. ~/fabric_data/), and also the file gene_bg_scores/<version>.tar.gz that you will need to extract into your gene background scores directory (e.g. for ~/fabric_data/gene_bg_scores/ you will need to end up with the JSON files inside the directory ~/fabric_data/gene_bg_scores/<version>/). If you use these pre-calculated scores, you will also need to make sure that geneffect works with the same gene annotations that created these files, so you will need to work with somewhat older versions of genenames and UniProt that are found in the data/ subdirectory at this FTP site (refer to your geneffect's configuration file, typically found at ~/.geneffect_config.py, for instructions where to update these data files). **It is absolutely crucial to make sure that geneffect and firm are configured in the exact same way they were configured when these effect scores were originally calculated, or else you may end up with false positives. Therefore, if you choose to use the pre-calculated scores, it is highly recommended that you install fabric using the automatic installation script (see below).**


Analyzing cancer data (MAF format)
----------------------------------

FABRIC can be easily used to analyze mutations provided in the MAF format. 

For example, you can download all cancer somatic mutations provided by TCGA from the GDC data portal (https://portal.gdc.cancer.gov/repository) by applying the following filters:

1) Select "MuTect2 ..." as the "Workflow Type"
2) Select "MAF" as the "Data Format"
3) Select "open" Access

Download the resulted query in TSV format to obtain a .tar.gz file with the requested MAF files. You can then process this tar file with the code provided in the "TCGA & ExAC Analysis" Jupyter Notebook provided in this GitHub repository (go to the "Combine GDC's downloaded tar file into a single MAF file" section in that notebook) in order to obtain a final MAF file (gdc_combined.maf) you can work with.

Of course, this is just an example to allow you to test FABRIC. If you want, you can start right off with your own MAF file.
FABRIC accepts as a valid MAF file any tab-delimited file with the following headers: 

1) **Variant_Type** (to filter only "SNP" mutations)
2) **Tumor_Seq_Allele1** (assumed to be the reference nucleotide of each SNP)
3) **Tumor_Seq_Allele2** (assumed to be the alternative nucleotide of each SNP)
4) **tcga_project** (assumed to be the TCGA project each mutation belongs to; required only when running a per-project analysis, and can be omitted in a pan-cancer analysis).



**Step 1: Calculate the effect scores of the observed mutations**

First, FABRIC needs to calculate the effect scores of all the mutations in your MAF file. This is done by running the following command:

        set_maf_gene_effect_scores --reference-genome=<GRCh38 or hg19> --gene-dataset-dir=<DIR PATH> --maf-file=<FILE PATH> --effect-scores-output-csv-file=<FILE PATH>

Where:

* --gene-dataset-dir should be the same directory you provided when running set_gene_bg_scores.
* --maf-file should be given your input MAF file.
* --effect-scores-output-csv-file should be given the file where you want to save the output effect scores calculated by the framework.

For example, if you want to analyze the TCGA dataset obtained above (gdc_combined.maf), assuming your working directory is ~/tcga_and_exac_analysis, you will run the command: 

        set_maf_gene_effect_scores --reference-genome=GRCh38 --gene-dataset-dir=~/fabric_data/ --maf-file=~/tcga_and_exac_analysis/gdc_combined.maf --effect-scores-output-csv-file=~/tcga_and_exac_analysis/gdc_effect_scores.csv

The reference genome GRCh38 is provided because that's the version used by TCGA.
When finished, the file ~/tcga_and_exac_analysis/gdc_effect_scores.csv will be created with the effect scores of all the relevant mutations in gdc_combined.maf.

Note that when analyzing millions of mutations, this can be a pretty lengthy process as well. It is recommended to run it in the background, and to use a machine with many CPUs. 


**Step 2: Run the analysis**

After all the preparations are finished (the background and observed effect scores have been calculated), you are ready to run the actual analysis. This is done by running the following command:

        analyze_maf_genes --reference-genome=<GRCh38 or hg19> --gene-dataset-dir=<DIR PATH> --maf-file=<FILE PATH> --effect-scores-csv-file=<FILE PATH> --gene-bg-scores-dir=<DIR PATH> --output-dir=<DIR PATH> [--only-combined] [--analyze-diff]

Where:

* --gene-dataset-dir, --maf-file, --effect-scores-csv-file and --gene-bg-scores-dir should be the same directories and files as in the previous commands.
* --output-dir should be given the directory in which you want to save the results.
* --only-combined and --analyze-diff are optional flags. Use the first to run only a combined (pan-cancer) analysis; use the second to also run an analysis of differences cross cancer types.

For example, to keep working on the TCGA dataset, you will run the command: 

        analyze_maf_genes --reference-genome=GRCh38 --gene-dataset-dir=~/fabric_data/ --maf-file=~/tcga_and_exac_analysis/gdc_combined.maf --effect-scores-csv-file=~/tcga_and_exac_analysis/gdc_effect_scores.csv --gene-bg-scores-dir=~/fabric_data/gene_bg_scores/ --output-dir=~/tcga_and_exac_analysis/gdc_results/ --analyze-diff

When finished, the directory ~/tcga_and_exac_analysis/gdc_results/ will be created and filled with the analysis results. A CSV file with the results of all analyzed genes will be created for each TCGA project (cancer-type), and another CSV (combined.csv) for the combined (pan-cancer) analysis. Since you also provided the --analyze-diff flag, a diff.csv file will be created as well.


Analyzing genetic variation (VCF format)
----------------------------------------

In addition to MAF format, FABRIC can also process a list of variants given in VCF format. 

For example, you can download all the variants observed in the healthy human population from ExAC (http://exac.broadinstitute.org/) in VCF format. The file is available at:
ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz 


**Step 1: Parse the VCF file and calculate the effect scores of the observed variants**

First, FABRIC needs to parse the VCF file, and to calculate the effect scores of all the relevant variants in it. This is done by running the following command:

        create_vcf_dataset --reference-genome=<GRCh38 or hg19> --gene-dataset-dir=<DIR PATH> --vcf-file=<FILE PATH> --output-csv-file=<FILE PATH> [--only-pass]

Where:

* --gene-dataset-dir is the same as with the MAF format.
* --vcf-file should be given the input VCF file.
* --output-csv-file should be given the file where you want to save the output processed dataset with the effect scores calculated by the framework.
* --only-pass is an optional flag to take only variants with a "PASS" filter.

For example, if you want to analyze the ExAC dataset obtained above (ExAC.r1.sites.vep.vcf.gz), assuming your working directory is ~/tcga_and_exac_analysis, you will run the command: 

        create_vcf_dataset --reference-genome=hg19 --gene-dataset-dir=~/fabric_data/ --vcf-file=~/tcga_and_exac_analysis/ExAC.r1.sites.vep.vcf.gz --output-csv-file=~/tcga_and_exac_analysis/exac_variants.csv --only-pass

The reference genome hg19 is provided because that's the version used by ExAC (ExAC's VCF headers contains: assembly=GRCh37.p13).
When finished, the file ~/tcga_and_exac_analysis/exac_variants.csv will be created with the processed dataset.

Note that when analyzing millions of variants, this can be a pretty lengthy process. It is recommended to run it in the background, and to use a machine with many CPUs. 


**Step 2: Run the analysis**

After all the background and observed effect scores have been calculated, you are ready to run the actual analysis. This is done by running the following command:

        analyze_vcf_genes --reference-genome=<GRCh38 or hg19> --gene-dataset-dir=<DIR PATH> --input-csv-file=<FILE PATH> --gene-bg-scores-dir=<DIR PATH> --output-csv-file=<FILE PATH>

Where:

* --gene-dataset-dir and --gene-bg-scores-dir should be the same as in the previous commands.
* --input-csv-file should be the output of the create_vcf_dataset command.
* --output-csv-file should be given the file where you want to save the results.

For example, to keep working on the ExAC dataset, you will run the command: 

        analyze_vcf_genes --reference-genome=hg19 --gene-dataset-dir=~/fabric_data/ --input-csv-file=~/tcga_and_exac_analysis/exac_variants.csv --gene-bg-scores-dir=~/fabric_data/gene_bg_scores/ --output-csv-file=~/tcga_and_exac_analysis/exac_results.csv

When finished, the file ~/tcga_and_exac_analysis/exac_results.csv will be created with the results of all analyzed genes.


Analyzing other types of data
-----------------------------

FABRIC is currently equipped with commandline scripts for processing data provided only in the MAF or VCF format. However, the API provided by the fabric Python module is quite generic, and can be used to write your own custom code to handle any kind of data. 


Installation
============

Dependencies:

* numpy
* scipy
* pandas
* biopython
* scikit-learn
* statsmodels
* geneffect (https://github.com/nadavbra/geneffect)
* firm (https://github.com/nadavbra/firm)


Automatic installation (using the installation script)
----------

    >>> wget https://raw.githubusercontent.com/nadavbra/fabric/master/install_fabric.sh
    >>> chmod a+x install_fabric.sh
    >>> ./install_fabric.sh
    
And follow the script's instructions.
    
    
Manual installation
----------

Make sure that geneffect and firm are properly installed and configured.

Clone the project and run:

    python setup.py install
    

Cite us
=======

If you use FABRIC as part of work contributing to a scientific publication, we ask that you cite our paper: Nadav Brandes, Nathan Linial, Michal Linial, Quantifying gene selection in cancer through protein functional alteration bias, Nucleic Acids Research, gkz546, https://doi.org/10.1093/nar/gkz546
