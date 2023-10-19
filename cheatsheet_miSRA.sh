## Install
# Install anaconda
conda
conda install -c anaconda python 
conda create --name miSRA python=3.11 #creates an environment
pip3 install miSRA    
miSRA

## Use

conda activate miSRA  
mkdir miSRA_example
cd miSRA_example  
miSRA -c config.json
#miSRA --config config.json


# mv hairpin_hsa.fa.txt hairpin_hsa.fa
# mv mature_hsa.fa.txt mature_hsa.fa   
