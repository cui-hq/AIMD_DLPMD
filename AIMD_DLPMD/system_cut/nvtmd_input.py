import os
import shutil
import glob
files=glob.glob("nvt/*")
for file in files:
    b=os.system("tail -n $(head -n 1 %s/coord_num.xyz) %s/R_3-pos-1.xyz > %s/coord.xyz" %(file,file,file))
    b=os.system("tail -n 1 %s/R_3-1.cell | awk '{print $3}' > cell " %(file))
    b=os.system("grep CHARGE %s/input.inp | head -n 1 | awk '{print $2}' > charge " %(file))
    r=open('cell','r')
    a=r.readlines()
    cell=float(a[0].strip())
    r=open('charge','r')
    a=r.readlines()
    charge=float(a[0].strip())

    with open("input_nvt.inp", 'w+') as f:
        r=open('head_nvt.txt','r')
        a=r.readlines()
        for k in a:
            f.write(k)
        f.write("%7s%15.8f%15.8f%15.8f\n" %("A",cell,0,0))
        f.write("%7s%15.8f%15.8f%15.8f\n" %("B",0,cell,0))
        f.write("%7s%15.8f%15.8f%15.8f\n" %("C",0,0,cell))
        r=open('middle_nvt.txt','r')
        a=r.readlines()
        for k in a:
            f.write(k)
        f.write("    CHARGE    %1d\n" %(charge))
        r=open('tail_nvt.txt','r')
        a=r.readlines()
        for k in a:
            f.write(k)

    f.close()

    shutil.move(r"input_nvt.inp",r"%s/input_nvt.inp" %(file))
    os.system("rm %s/R_3*" %(file))
#    os.system("rm %s/input_md.inp" %(file))
