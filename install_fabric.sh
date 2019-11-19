#! /bin/tcsh -e

# Determine the Python version to use.

# If a the environment variable PYTHON_CMD exists, then just use it.
if ($?PYTHON_CMD) then
    set python_cmd = $PYTHON_CMD
    goto determine_final_python_version
endif

ask_default_python_version:
printf "Do you want to install geneffect in your default Python version (`python --version |& cat`)? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set python_cmd = "python"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative version of Python: "
    set python_cmd = $<
else
    echo "Please specify either Yes or No."
    goto ask_default_python_version
endif

determine_final_python_version:
echo "Will install geneffect for the Python version `$python_cmd --version |& cat` ($python_cmd)"


# Determine a temp working directory.

# If a the environment variable PYTHON_CMD exists, then just use it.
if ($?TEMP_DIR) then
    set temp_dir = $TEMP_DIR
    goto determine_final_temp_dir
endif

ask_default_temp_dir:
printf "Do you want to use /tmp as a temporary working directory? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set temp_dir = "/tmp"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative directory path: "
    set temp_dir = $<
else
    echo "Please specify either Yes or No."
    goto ask_default_temp_dir
endif

determine_final_temp_dir:
set temp_dir = `echo $temp_dir | sed 's:/*$::'`
echo "Will use $temp_dir as a temporary working directory."


# Determine the data directory for geneffect.

# If a the environment variable DATA_DIR exists, then just use it.
if ($?DATA_DIR) then
    set data_dir = $DATA_DIR
    goto determine_final_data_dir
endif

ask_data_dir:
printf "Do you want to use ~/data as the directory for all the data files required by geneffect? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set data_dir = "~/data"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative directory path: "
    set data_dir = $<
else
    echo "Please specify either Yes or No."
    goto ask_data_dir
endif

determine_final_data_dir:
set data_dir = `echo $data_dir | sed 's:/*$::'`
echo "Will use $data_dir as the directory for the data files of geneffect."


# Determine the data directory for fabric.

# If a the environment variable FABRIC_DATA_DIR exists, then just use it.
if ($?FABRIC_DATA_DIR) then
    set fabric_data_dir = $FABRIC_DATA_DIR
    goto determine_final_fabric_data_dir
endif

ask_fabric_data_dir:
printf "Do you want to use ~/fabric_data as the directory for all the data files required by geneffect? Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set fabric_data_dir = "~/fabric_data"
else if ($response == "n" || $response == "no") then
    printf "Please specify an alternative directory path: "
    set fabric_data_dir = $<
else
    echo "Please specify either Yes or No."
    goto ask_fabric_data_dir
endif

determine_final_fabric_data_dir:
set fabric_data_dir = `echo $fabric_data_dir | sed 's:/*$::'`
echo "Will use $fabric_data_dir as the directory for the data files of fabric."
mkdir -p $fabric_data_dir


# Determine whether to use pre-calculated gene effect scores.

# If a the environment variable USE_PRECALCULATED_EFFECT_SCORES exists, then just use it.
if ($?USE_PRECALCULATED_EFFECT_SCORES) then
    set use_precalculated_effect_scores = $USE_PRECALCULATED_EFFECT_SCORES
    goto determine_whether_to_use_precalculated_effect_scores
endif

ask_whether_to_use_precalculated_effect_scores:
printf "Do you want to use use pre-calculated gene effect scores? Choosing that option will spare you a lot of CPU time, but it will override your geneffect data files to a somewhat older version. Yes/No (Default: Yes): "
set response = $<
set response = `echo $response | tr "[:upper:]" "[:lower:]"`

if ($response == "" || $response == "y" || $response == "yes") then
    set use_precalculated_effect_scores = "yes"
else if ($response == "n" || $response == "no") then
    set use_precalculated_effect_scores = "no"
else
    echo "Please specify either Yes or No."
    goto ask_whether_to_use_precalculated_effect_scores
endif

determine_whether_to_use_precalculated_effect_scores:

if ($use_precalculated_effect_scores == "yes") then
    echo "Will use pre-calculated gene effect scores."
else
    echo "Will not use pre-calculated gene effect scores."
endif


# Install general dependencies.

echo "Installing dependencies..."
$python_cmd -m pip install numpy scipy pandas biopython sklearn statsmodels cython interval_tree


# Install firm (which will install geneffect).

echo "Installing firm & geneffect..."
wget https://raw.githubusercontent.com/nadavbra/firm/master/install_firm.sh -O ~/install_firm.sh ${temp_dir}/install_firm.sh
chmod a+x ${temp_dir}/install_firm.sh
setenv PYTHON_CMD $python_cmd
setenv TEMP_DIR $temp_dir
setenv DATA_DIR $data_dir
${temp_dir}/install_firm.sh
rm -f ${temp_dir}/install_firm.sh


# Install fabric.

mkdir -p ${temp_dir}/fabric
git clone https://github.com/nadavbra/fabric.git ${temp_dir}/fabric
cd ${temp_dir}/fabric
$python_cmd ./setup.py install
cd -
rm -fr ${temp_dir}/fabric


# Downloaded the pre-calculated effect scores and set the other geneffect data files to compatible versions (if asked to do so).

if ($use_precalculated_effect_scores == "yes") then
    
    echo "Downloading the gene datasets..."
    wget ftp://ftp.cs.huji.ac.il/users/nadavb/fabric_data/genes_hg19.csv -O ${fabric_data_dir}/genes_hg19.csv
    wget ftp://ftp.cs.huji.ac.il/users/nadavb/fabric_data/genes_GRCh38.csv -O ${fabric_data_dir}/genes_GRCh38.csv
    
    mkdir -p ${fabric_data_dir}/gene_bg_scores
    
    echo "Downloading the background scores for GRCh38..."
    wget ftp://ftp.cs.huji.ac.il/users/nadavb/fabric_data/gene_bg_scores/GRCh38.tar.gz -O ${fabric_data_dir}/gene_bg_scores/GRCh38.tar.gz
    gunzip ${fabric_data_dir}/gene_bg_scores/GRCh38.tar.gz
    tar -xf ${fabric_data_dir}/gene_bg_scores/GRCh38.tar -C ${fabric_data_dir}/gene_bg_scores
    rm -f ${fabric_data_dir}/gene_bg_scores/GRCh38.tar
    
    echo "Downloading the background scores for hg19..."
    mkdir -p ${fabric_data_dir}/gene_bg_scores/hg19
    wget ftp://ftp.cs.huji.ac.il/users/nadavb/fabric_data/gene_bg_scores/hg19.tar.gz -O ${fabric_data_dir}/gene_bg_scores/hg19.tar.gz
    gunzip ${fabric_data_dir}/gene_bg_scores/hg19.tar.gz
    tar -xf ${fabric_data_dir}/gene_bg_scores/hg19.tar -C ${fabric_data_dir}/gene_bg_scores
    rm -f ${fabric_data_dir}/gene_bg_scores/hg19.tar
    
    echo "Replacing geneffect data files to older compatible versions..."
    wget ftp://ftp.cs.huji.ac.il/users/nadavb/fabric_data/data/uniprot/uniprot_human_reviewed.xml.gz -O ${data_dir}/uniprot/uniprot_human_reviewed.xml.gz
    wget ftp://ftp.cs.huji.ac.il/users/nadavb/fabric_data/data/genenames/non_alt_loci_set.json -O ${data_dir}/genenames/non_alt_loci_set.json
endif


echo "Succefully installed fabric."