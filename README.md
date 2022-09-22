# multi-omics-R-loops
This repository is for the paper "Multi-omics to characterize the functional relationships of R-loops with epigenetic modifications, RNAPII transcription and gene expression"

The python packages required to be used in this script include:

```
matplotlib                3.2.2                        
matplotlib-base           3.2.2            
numpy                     1.21.5           
pandas                    1.3.5            
python                    3.7.12          
python-dateutil           2.8.2              
python_abi                3.7                     
scikit-learn              1.0.2            
scipy                     1.7.3            
setuptools                59.8.0           
sqlite                    3.37.0               
statsmodels               0.13.2
catboost                  0.26.1
lightgbm                  3.2.1
xgboost                   1.5.0
```




Running environment installation:

1. Download and install Anaconda 3 on your system. Once complete, restart the shell to initialize the conda installation.
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh
```

2. Most of the packages can be setup through conda directly. To avoid version incompatibilities, it’s best to install most of the packages in their own conda environment.
```
conda create -n multi-omic python=3
```
3. Activate the created environment and install the packages available from conda.
```
conda activate multi-omic
conda install -c bioconda matplotlib scikit-learn statsmodels catboost lightgbm xgboost
```
How to use it?
```
multi-omics R-loop prediction model training and evaluation for one specific bin:
	python multi-bin-histone-predict-R-loop.py example_bin.bed  example_bin_result_1.txt > example_bin_result_2.txt
```
To note, example_bin.bed is a matrix file where a row represents one specific bin region of all protein coding TSSs or TTSs, non-coding TSSs or TTSs, and enhancers and column represent omics type, and cell represents the quantification of corresponding omics type in corresponding bins. 

Additionally, it generates five png files:
```
example_bin_1.png    correlation matrix plot for the bin
example_bin_2.png    scaled algorithm comparison based on mean squared error
example_bin_3.png    scaled algorithm comparison based on root mean squared error
example_bin_4.png    gridSearchCV implement for extra trees regressor models
example_bin_5.png    output the relative importance of epigenetics marks and predictive R squared
```
To note, to carry out on R-loop prediction based on single-omics or combinatorial-omics, the only thing you need to do is adjust the input in the multi-bin-histone-predict-R-loop.py.

Citation
```
If you find that tool useful in your research, please consider citing our paper.
```
