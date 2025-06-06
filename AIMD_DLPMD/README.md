(A) Description
An example of Generating DLP model.

(B) Requirements
#python package
numpy 
mdtraj
json
os
shutil
dpdata

#xTB-version
xtb > 6.5.1

#cp2k-version
cp2k > 7.1

#deepmd-version
deepmd > 2.2.0

#Gaussian
Gaussian 16

(C) usage
examples:

(a) Generating AIMD files(coordinate files(.xyz),input files(.inp)) in folders (cut_*) and Performing AIMD simulations

Run:

   cd system_cut
   
   python system_cut.py
   
   sh sub_xtb.sh
   
   sh nptxyz_input.sh
   
   sh sub_npt.sh
   
   mkdir nvt
   
   python nvtmd_input.py
   
   cd nvt 
   
   cp ../sub_nvt.sh .
   
   sh sub_nvt.sh
   
(b) preparing training data

Run:

   cp ../sub_trans.sh .
   
   sh sub_trans.sh
   
(c) Searching transition states by Gaussian16 and labeling these structures corresponding to reaction pathways using cp2k.

 
(d) training model

Run:

   cd train/01.stage
   
   dp train input.json # only structures corresponding to reaction pathways in training dataset.
   
   dp freeze -o 01.model
   
   cd train.02.stage
   
   dp train input.json -t 01.model # structures corresponding to reaction pathways and AIMD trajectories in training dataset.
   
   dp freeze -o 02.model
   

