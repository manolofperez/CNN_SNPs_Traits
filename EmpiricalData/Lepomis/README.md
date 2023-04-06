# ***Lepomis* scripts and datasets**
This folder contains the empirical datasets of the fish *Lepomis megalotis* species complex  
used to test our deep learning approach. 

The folder contains the following files:

simulate_ms_Lepomis.py - python script to simulate coalescent trees and segregating sites (saved as NumPy arrays)
for the *Lepomis* dataset.

Train_Test_Predict_Lepomis.ipynb -  python notebook containing code and outputs for CNN training, cross-validation 
and prediction of the most likely model using empirical data.

Pred_Emp_Comb_BM_Predictions.txt - Predictions for the model trained with the Combined SNPs + BM traits.

Pred_Emp_Comb_OU_Predictions.txt - Predictions for the model trained with the Combined SNPs + OU traits.

Pred_Emp_SNP_Predictions.txt - Predictions for the model trained with SNPs.

Pred_Emp_traits_BM_Predictions.txt - Predictions for the model trained with the BM traits.

Pred_Emp_traits_OU_Predictions.txt - Predictions for the model trained with the OU traits.

input_SNPs.txt.gz - Segregating sites from the empirical dataset.

input_traits_Laqu.txt - Traits for *Lepomis aquilensis*.

input_traits_Laqu_Lpel.txt - Traits for putative hybrids between *Lepomis aquilensis* and *Lepomis peltastes*.

input_traits_Lmeg.txt - Traits for *Lepomis magalotis*.

input_traits_Loua.txt - Traits for *Lepomis* sp. Ouachita.

input_traits_Lozk.txt - Traits for *Lepomis* sp. Ozark.

input_traits_Lozk_Lmeg.txt - Traits for putative hybrids between *Lepomis* sp. Ozark and *Lepomis megalotis*.

input_traits_Lpel.txt - Traits for *Lepomis peltastes*.

input_traits_Lpel_Lmeg.txt - Traits for putative hybrids between *Lepomis peltastes* and *Lepomis megalotis*.

input_traits_Lsol.txt - Traits for *Lepomis solis*.

input_traits_Lsol_Lmeg.txt - Traits for putative hybrids between *Lepomis solis* and *Lepomis megalotis*.

TrainedModels - Folder containing all the trained models resulting from running the file "Train_Test_Predict_Lepomis.ipynb".
