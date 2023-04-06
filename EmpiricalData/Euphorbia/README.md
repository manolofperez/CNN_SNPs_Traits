# ***Euphorbia* scripts and datasets**
This folder contains the empirical datasets of the plant *Euphorbia balsamifera* species complex  
used to test our deep learning approach. 

The folder contains the following files:

simulate_ms_Euphorbia.py - python script to simulate coalescent trees and segregating sites (saved as NumPy arrays)
for the *Euphorbia* dataset.

Train_Test_Predict.ipynb -  python notebook containing code and outputs for CNN training, cross-validation 
and prediction of the most likely model using empirical data.

Pred_Emp_Comb_BM_Predictions.txt

Pred_Emp_Comb_OU_Predictions.txt

Pred_Emp_SNP_Predictions.txt

Pred_Emp_traits_BM_Predictions.txt

Pred_Emp_traits_OU_Predictions.txt

input_SNPs.txt - Segregating sites from the empirical dataset.

input_cont.txt

input_traits_ade.txt

input_traits_bal.txt

input_traits_sep.txt

TrainedModels - Folder containing all the trained models resulting from running the file "Train_Test_Predict.ipynb".
