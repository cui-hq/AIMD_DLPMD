import mdtraj
import numpy as np
import json
import os
import shutil
traj=mdtraj.load("confout.gro",top="./N_3_RBI_ALK_Br_CCL.gro")
topology=traj.topology
BOX_CELL = traj.unitcell_lengths[0][0]*10
residue_list=[residue.index for residue in topology.residues]
datastore=[]
doc_list=[]
top_list=[]
top_MOL_index_list=[]

res_list=[residue for residue in topology.residues]
res_type_list=list(set([residue.name for residue in res_list]))

with open("param.json") as fp:
    jdata=json.load(fp)

datastore=[]
doc_list=[]
top_list=[]
top_MOL_index_list=[]

with open(jdata["sys_bond_top"][0], 'r') as top_file:
    for line in top_file:
        if (line.strip() != ""):
            datastore.append(line)

for i in datastore:
    if (i.split()[0] == "#include") :
        doc_list.append(eval(i.split()[1]))

for doc in doc_list:
    with open(doc, 'r') as top_file:
        for line in top_file:
            top_list.append(line)

top_MOL_index_list = [[i] for i in (np.where(np.array(top_list) == "[ moleculetype ]\n")[0] + 2)]

for index in top_MOL_index_list:
    if top_list[index[0]+5].strip() == "" or top_list[index[0]+5].split()[0] == ";":
        index.append([])
    else:
        for num in range(index[0],len(top_list)):
            if top_list[num] == "[ bonds ]\n":
                begin_num = num + 2
                break
        for num in range(begin_num,len(top_list)):
            if top_list[num].strip() == "":
                final_num = num - 1
                break
        bonds_list = []

        for num in range(begin_num,final_num+1):
            bonds_list.append([int(top_list[num].split()[0]),int(top_list[num].split()[1])])
        index.append(bonds_list)

break_atom_id=[]
break_residue_id=[]

for res_id in residue_list:
    atom_index = [atom.index for atom in topology.residue(res_id).atoms]
    for index_bonds in top_MOL_index_list:
        if topology.residue(res_id).name == top_list[index_bonds[0]].split()[0]:
            bonds_index_list = [[bond[0]+atom_index[0]-1,bond[1]+atom_index[0]-1] for bond in index_bonds[1]]

    
    for break_mol in jdata["break_bond"]:
        if topology.residue(res_id).name == break_mol:
            for bond_name in jdata["break_bond"][break_mol]:
                for bond_index in bonds_index_list:
                    if [topology.atom(bond_index[0]).name,topology.atom(bond_index[1]).name] == bond_name or [topology.atom(bond_index[1]).name,topology.atom(bond_index[0]).name] == bond_name:
                        distance_bond = mdtraj.compute_distances(traj,atom_pairs=[bond_index])[0]*10
                        if distance_bond[0] > 3:# BOX_CELL/2:
                            if topology.atom(bond_index[0]).name == bond_name[0]:
                                break_atom_id.append(bond_index[0])
                                break_residue_id.append(res_id)
                            elif topology.atom(bond_index[1]).name == bond_name[0]:
                                break_atom_id.append(bond_index[1])
                                break_residue_id.append(res_id)



#X_BOX=[(0,2),(1,3),(2,4),(3,5)]
#Y_BOX=[(1,3),(0,2),(2,4),(3,5)]
#Z_BOX=[(1,3),(0,2),(2,4),(3,5)]
X_BOX=[(0,2),(2,4)]
Y_BOX=[(0,2),(2,4)]
Z_BOX=[(0,2),(2,4)]

BOX_cut=[]
for i in X_BOX:
    for j in Y_BOX:
        for z in Z_BOX:
            BOX_cut.append((i,j,z))
#box_cut_off = 2
#print(BOX_cut)
#print(len(BOX_cut))
box_cut={}
for i in BOX_cut:
        box_cut[i] = []

for key in box_cut.keys():
    for res_id in residue_list:
        if res_id in break_residue_id:
            atom_index = [atom.index for atom in topology.residue(res_id).atoms]
            atom_res_delete=[]
            for atom_id in atom_index:
                if atom_id in break_atom_id:
                    atom_res_delete.append(atom_id)
            for atom_id in atom_res_delete:
                atom_index.remove(atom_id)
            res_coord_x = np.array([traj.xyz[0, id, 0] for id in atom_index])
            res_coord_y = np.array([traj.xyz[0, id, 1] for id in atom_index])
            res_coord_z = np.array([traj.xyz[0, id, 2] for id in atom_index])
            
            x_max = np.max(res_coord_x)
            x_min = np.min(res_coord_x)
            y_max = np.max(res_coord_y)
            y_min = np.min(res_coord_y)
            z_max = np.max(res_coord_z)
            z_min = np.min(res_coord_z)
            
            x_center = (x_max + x_min)/2
            y_center = (y_max + y_min)/2
            z_center = (z_max + z_min)/2


            if ((x_center > key[0][0] and x_center < key[0][1]) or ((x_center + 4) > key[0][0] and (x_center +4 < key[0][1]))) and \
            ((y_center > key[1][0] and y_center < key[1][1]) or ((y_center + 4) > key[1][0] and (y_center +4 < key[1][1]))) and \
            ((z_center > key[2][0] and z_center < key[2][1]) or ((z_center + 4) > key[2][0] and (z_center +4 < key[2][1]))):
                box_cut[key].append(res_id)


        else:
            atom_index = [atom.index for atom in topology.residue(res_id).atoms]
            res_coord_x = np.array([traj.xyz[0, id, 0] for id in atom_index])
            res_coord_y = np.array([traj.xyz[0, id, 1] for id in atom_index])
            res_coord_z = np.array([traj.xyz[0, id, 2] for id in atom_index])

            x_max = np.max(res_coord_x)
            x_min = np.min(res_coord_x)
            y_max = np.max(res_coord_y)
            y_min = np.min(res_coord_y)
            z_max = np.max(res_coord_z)
            z_min = np.min(res_coord_z)

            x_center = (x_max + x_min)/2
            y_center = (y_max + y_min)/2
            z_center = (z_max + z_min)/2

            if ((x_center > key[0][0] and x_center < key[0][1]) or ((x_center + 4) > key[0][0] and (x_center +4 < key[0][1]))) and \
                ((y_center > key[1][0] and y_center < key[1][1]) or ((y_center + 4) > key[1][0] and (y_center +4 < key[1][1]))) and \
                ((z_center > key[2][0] and z_center < key[2][1]) or ((z_center + 4) > key[2][0] and (z_center +4 < key[2][1]))):
                    box_cut[key].append(res_id)
 #       print(x_center)
 #       print(cut_coord)

for key in box_cut.keys():
    res_index = box_cut[key]
    res_name_list = [topology.residue(res_id).name for res_id in res_index]
    tot_charge = 0
    atom_tot_index=[]
#    print(res_index)
    for res_id in res_index:
        tot_charge += jdata["residue_charge"][topology.residue(res_id).name]
        atom_index = [atom.index for atom in topology.residue(res_id).atoms]
        atom_remove=[]
        for atom_id in atom_index:
            if atom_id in break_atom_id:
                atom_remove.append(atom_id)
                for charge_atom in jdata["break_atom"]:
                    if topology.atom(atom_id).name==charge_atom:
                        tot_charge = tot_charge - jdata["break_atom"][charge_atom]
        for atom_id in atom_remove:
            atom_index.remove(atom_id)

        atom_tot_index.extend(atom_index)
    
    for index in break_atom_id:
        if ((traj.xyz[0, index, 0] > key[0][0] and traj.xyz[0, index, 0] < key[0][1]) or ((traj.xyz[0, index, 0] +4) > key[0][0]  and (traj.xyz[0, index, 0] +4) < key[0][1])) and \
        ((traj.xyz[0, index, 1] > key[1][0] and traj.xyz[0, index, 1] < key[1][1]) or ((traj.xyz[0, index, 1] +4) > key[1][0]  and (traj.xyz[0, index, 1] +4) < key[1][1])) and \
        ((traj.xyz[0, index, 2] > key[2][0] and traj.xyz[0, index, 2] < key[2][1]) or ((traj.xyz[0, index, 2] +4) > key[2][0]  and (traj.xyz[0, index, 2] +4) < key[2][1])):
#            print(traj.xyz[0, index, :])
#            print(key)
            atom_tot_index.append(index)
            for charge_atom in jdata["break_atom"]:
                if topology.atom(index).name==charge_atom:
                    tot_charge = tot_charge + jdata["break_atom"][charge_atom]
#            print("-----")
#            print(index)
#            print(key)
#            print("----")
#    print(tot_charge)


    
    num_atoms = len(atom_tot_index)
#    print(key)
#    print(key[0][0],key[1][0],key[2][0])
    os.mkdir("cut_%s_%s_%s" %(str(key[0][0]),str(key[1][0]),str(key[2][0])))
    with open('coord_num.xyz','w+') as f:
        f.write(str(num_atoms))
        f.write("\n")
        f.write("\n")
        for index in atom_tot_index:
            coord_x = traj.xyz[0, index, 0]*10 - key[0][0]*10
            dX = coord_x - 10
            while (abs(dX) > 20):
                moveX=(-dX/abs(dX))*40
                coord_x += moveX
                dX+=moveX
            coord_y = traj.xyz[0, index, 1]*10 - key[1][0]*10
            dY = coord_y - 10
            while (abs(dY) > 20):
                moveY=(-dY/abs(dY))*40
                coord_y += moveY
                dY+=moveY
            coord_z = traj.xyz[0, index, 2]*10 - key[2][0]*10
            dZ = coord_z - 10
            while (abs(dZ) > 20):
                moveZ=(-dZ/abs(dZ))*40
                coord_z += moveZ
                dZ+=moveZ
            if topology.atom(index).name == "CL1" or topology.atom(index).name == "CL2" or topology.atom(index).name == "CL":
                f.write(" %s   %f  %f   %f\n" %( "Cl",coord_x,coord_y,coord_z))
            elif topology.atom(index).name == "BR":
                f.write(" %s   %f  %f   %f\n" %( "Br",coord_x,coord_y,coord_z))
            else:
                f.write(" %s   %f  %f   %f\n" %( topology.atom(index).name[0:1],coord_x,coord_y,coord_z))
        f.close()

    with open('input_npt.inp','w+') as f:
        r=open('head_npt.txt','r')
        a=r.readlines()
        for k in a:
            f.write(k)
        f.write("    CHARGE    %d\n" %( tot_charge))
        r=open('tail_npt.txt','r')
        a=r.readlines()
        for k in a:
            f.write(k)
    f.close()
    with open('xtb.sh','w+') as f:
        f.write("#!/bin/bash\n")
        f.write("xtb coord_num.xyz --input md.inp --omd --chrg %s --gfn 1" %(tot_charge))
    f.close()


    shutil.move("coord_num.xyz","cut_%s_%s_%s" %(str(key[0][0]),str(key[1][0]),str(key[2][0])))
    shutil.move("input_npt.inp","cut_%s_%s_%s/input.inp" %(str(key[0][0]),str(key[1][0]),str(key[2][0])))
    shutil.move("xtb.sh","cut_%s_%s_%s/xtb.sh" %(str(key[0][0]),str(key[1][0]),str(key[2][0])))
    shutil.copy("md.inp","cut_%s_%s_%s/md.inp" %(str(key[0][0]),str(key[1][0]),str(key[2][0])))

