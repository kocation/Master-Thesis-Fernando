
import numpy as np
#TEST OF INTERPOLATION ROUTINE. COMPARE TO INTERP1 IN MATLAB


files=['FFA-W3-241_iced.txt','FFA-W3-301_iced.txt','FFA-W3-360_iced.txt','FFA-W3-480_iced.txt','FFA-W3-600_iced.txt','cylinder.txt']
#Initializing tables    
cl_tab=np.zeros([9,6])
cd_tab=np.zeros([9,6])
cm_tab=np.zeros([9,6])
aoa_tab=np.zeros([9])

f_stat=np.zeros([9,6])
cl_inv_tab=np.zeros([9,6])
cl_fs=np.zeros([9,6])
#Readin of tables. Only do this once at startup of simulation
for i in range(np.size(files)):
    aoa_tab[:],cl_tab[:,i],cd_tab[:,i],cm_tab[:,i],f_stat[:,i],cl_inv_tab[:,i],cl_fs[:,i] = np.loadtxt(files[i], skiprows=0).T

# Thickness of the airfoils considered
# NOTE THAT IN PYTHON THE INTERPOLATION REQUIRES THAT THE VALUES INCREASE IN THE VECTOR!

thick_prof=np.zeros(6)
thick_prof[0]=24.1;
thick_prof[1]=30.1;
thick_prof[2]=36;
thick_prof[3]=48;
thick_prof[4]=60;
thick_prof[5]=100;





def force_coeffs_10MW_iced(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab,f_stat,cl_inv,cl_fs):
    cl_aoa=np.zeros([1,6])
    cd_aoa=np.zeros([1,6])
    cm_aoa=np.zeros([1,6])

    fs_aoa=np.zeros([1,6])
    cl_inv_aoa=np.zeros([1,6])
    cl_fs_aoa=np.zeros([1,6])
    

    #Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cl_tab[:,i])
        cd_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cd_tab[:,i])
        cm_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cm_tab[:,i])

        fs_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,f_stat[:,i])
        cl_inv_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cl_inv_tab[:,i])
        cl_fs_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cl_fs[:,i])
    
    #Interpolate to current thickness:
    cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    cm=np.interp (thick,thick_prof,cm_aoa[0,:])

    fs_stat = np.interp (thick,thick_prof,fs_aoa[0,:])
    cl_inv = np.interp (thick,thick_prof,cl_inv_aoa[0,:])
    cl_fs = np.interp (thick,thick_prof,cl_fs_aoa[0,:])


    return cl, cd, cm , fs_stat, cl_inv, cl_fs


"""
# Lets test it:
angle_of_attack=10 # in degrees
thick = 42 # in percent !
[clift,cdrag,cmom]=force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab)

print('cl:',clift)
print('cd:',cdrag)
print('cm:',cmom)
"""
