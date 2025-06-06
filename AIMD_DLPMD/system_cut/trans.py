import dpdata
x=dpdata.LabeledSystem('./',cp2k_output_name='./input_nvt.log', fmt='cp2k/aimd_output')
x.to('deepmd/npy', 'dpmd_npy')

