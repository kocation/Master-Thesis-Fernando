#imports
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../files/')  # Adjust the relative path accordingly

from Interpolation_of_coeffs_iced import *  # for cl and cd
from scipy.interpolate import * 
from Python_load_Mann import * # function for importing turbulence files
from simulate_wind_speed import * #function for simulating wind speed
import time # for time stamps


def unsteady_BEM_function_iced(time_interval = 100, dt = 0.15, vs = [11.4], Tower_influ=0, Wind_shear=0, Dynamic_wake=0, Dynamic_stall=0, Turbulence=0, Controller = 0, Control_Vect = []):
    #areolastic dynamics
    #------------------------
    
    #vs = [15]#np.linspace(4,25,4*(25-4)+1) #wind speed vector, fill in ove value as [8] for running at 8 m/s 
    #, or change to a range for seeing how pitch changes as a function of wind speed

    Pitches = np.zeros(len(vs))
    Cps = np.zeros(len(vs))
    #start of code 
    t0 = time.time() #time stamp at start

    #Setp up option for code (tower influence, dynamic loads etc)
    Echo = 1 #Indicates total run time of the script as well change of pitch

    N = time_interval/dt #to achieve certain time in s
    N = int(N)

    # random oscilating income wind speed
    #np.random.seed(42)  # Seed for reproducibility
    #mean_wind_speed = vs[0] #mean wind speed
    #wind_speed = mean_wind_speed + np.random.uniform(-5, 5, N)  # Oscillations 


    # Parameters for wind
    I = 0.10  # Turbulence intensity
    l = 340.2  # Length scale
    V_10 = vs[0]  # Mean wind speed at 10 meters height
    f_cutout = 0.5  # Cutoff frequency (Hz)
    T_max = time_interval  # Time period for 10 minutes simulation (600 seconds)
    df = 1 / T_max  # Frequency step
    f_p = np.arange(df, f_cutout, df)  # Frequency range

    wind_speed = simulate_wind_speed(I, l, V_10, f_cutout, T_max, dt)


    #)set up constants and empty vectors
    if True: #in if statemt so that the entire setp up can be collapsed
        ang = 0 #yaw angle degree
        H=119 #hub height
        LS=7.1 #shaft length
        R=89.15 #balde radius
        B=3 #blade number
        rho = 1.225 #air density
        blade_file = np.genfromtxt("bladedat.txt")
        b_elm = blade_file[:,0]
        c = blade_file[:,2]
        twist = blade_file[:,1]
        t_c = blade_file[:,3]

        rt = np.array([H,0,0]) #tower heigh vector
        t = np.array([0]) #time counter
        

        omega = np.ones(int(N)+1) * 7.229*np.pi/30 #rotation speed

        
        #4096 same as turbulence box
        #N = np.pi*2/(dt*omega) #total steps to achieve a full rotation
        #N = N*B #Makes sure each blade has done a full rotation

        
        v = 0.2 #wind shear


        theta_base = 0
        THETA = np.zeros((B,int(N)+1)) #3,by time steps position matrix + init
        for j in(range(B)): #initial balde position set up
            THETA[j][0] = (j*2)*np.pi/B


        # Dynamic wake constants
        k_D = 0.6

        #set up angles
        theta_pitch = np.deg2rad(0) #deg to rad
        theta_cone = np.deg2rad(0) #deg to rad
        theta_yaw = np.deg2rad(ang) # deg to rad
        theta_tilt = np.deg2rad(0) # deg to rad

        #Matrices
        Wint_z = np.zeros((int(N)+1,B,len(b_elm))) #induced wind for dynamic wake
        Wint_y = np.zeros((int(N)+1,B,len(b_elm)))

        W_z = np.zeros((int(N)+1,B,len(b_elm))) #same as above
        W_y = np.zeros((int(N)+1,B,len(b_elm)))

        Wy = np.zeros((int(N)+1,B,len(b_elm))) #matrices for induced wind 
        Wz = np.zeros((int(N)+1,B,len(b_elm)))

        Wz_quasi = np.zeros((int(N)+1,B,len(b_elm))) #matrices for quasi induced wind 
        Wy_quasi = np.zeros((int(N)+1,B,len(b_elm)))

        pn = np.zeros((int(N)+1,B,len(b_elm))) #matrices for element forces 
        pt = np.zeros((int(N)+1,B,len(b_elm)))

        Cl = np.zeros((int(N)+1,B,len(b_elm))) #matrix for Cl
        fs = np.zeros((int(N)+1,B,len(b_elm)))

        A = np.zeros((int(N)+1,B,len(b_elm))) #filler matrix for debugging
        #ff = np.zeros((int(N)+1,B,len(b_elm))) #for prandtl tip loss

        x = np.zeros((int(N)+1,B,len(b_elm))) #matrices for coordinates
        y = np.zeros((int(N)+1,B,len(b_elm)))
        z = np.zeros((int(N)+1,B,len(b_elm)))

        """
        Vturb_x = np.zeros((int(N)+1,B,len(b_elm))) #matrices for turbulence intensities
        """
        #turbulence matricies
        Vturb_y = np.zeros((int(N)+1,B,len(b_elm))) #across
        Vturb_z = np.zeros((int(N)+1,B,len(b_elm))) #flow direction





        #Change of base matrices
        a1 = np.matrix([[1,0,0],
                    [0, np.cos(theta_yaw),np.sin(theta_yaw)],
                    [0,-1*np.sin(theta_yaw),np.cos(theta_yaw)]])

        a2 = np.matrix([[np.cos(theta_tilt),0,-1*np.sin(theta_tilt)],
                        [0,1,0],
                        [np.sin(theta_tilt),0,np.cos(theta_tilt)]])
        a3 = np.matrix([[1,0,0],
                        [0,1,0],
                        [0,0,1]])
        a12 = np.dot(np.dot(a3,a2),a1)
        a21 = np.transpose(a12)



        a34 = np.matrix([[np.cos(theta_cone),0,-np.sin(theta_cone)],
                        [0,1,0],
                        [np.sin(theta_cone),0,np.cos(theta_cone)]])
        a43 = np.transpose(a34)

        #turbuelence postition files
        dx = np.linspace(0,Ly,n2) #create a linspace with same definition as the simulation
        dy = np.linspace(0,Lz,n3) 

        #load turbulence files

        ut=load('sim1.bin',  N=(n1, n2, n3)) #z direction horizontal flow direction
        ushp = np.reshape(ut, (n1, n2, n3)) 

        vt=load('sim2.bin',  N=(n1, n2, n3)) #-y direction horizontal across
        vshp = np.reshape(vt, (n1, n2, n3)) 



        # wt=load('sim2.bin',  N=(n1, n2, n3)) #x directin vertical
        # wshp = np.reshape(wt, (n1, n2, n3)) 


        vvs = np.zeros((int(N)+1,B,len(b_elm))) #debug matrix
        theta_pitch = np.ones(N+1) * theta_pitch


    #if Controller:
    # controller set up
    #Constants 
    P_rated = 10.64e6 #Rated power
    T_rated = 1e6
    Kp = 1.5 #s
    KI = 0.64 #rad/rad
    K1 = np.deg2rad(14) #deg to rad
    omega_rated = 1.01 
    omega_ref = omega_rated + 0.01
    theta_dot = np.deg2rad(9) #maximum pitch angle speed possible
    theta_set_point = 0 
    
    #pitch limits
    theta_max = np.deg2rad(45) #max pitch angle
    theta_min = np.deg2rad(0) #min pitch angle

    K = P_rated/(omega_rated**3) #K val found from report
    
    I_rotor = 1.6*10**8 #moment of inertia of the rotor
    #Pi controller
    theta_P = np.zeros(N+1) 
    theta_I = np.zeros(N+1)
    M_aero = np.zeros(N+1)
    M_g = np.ones(N+1) * K*omega[0]**2 #generator moment



    #-----------------Start of loops
    for i in range(int(N)): #time loop
        #V0 = vs[0] #wind_speed[i] #wind speed
        V0 = wind_speed[i] 

        if Turbulence:
            fz = interp2d(dx, dy, ushp[i,:,:], kind='cubic') #interpolate stream wise turbulence
            fy = interp2d(dx, dy, vshp[i,:,:], kind='cubic') #interpolate cross wind turbulence


        t = np.append(t,t[i]+dt) #time vector

        try: #only works if the vector is set up correctly 
            if t[i]>=Control_Vect[0][0]: #if time is over time limit set, input angle control
                theta_pitch = Control_Vect[0][1]*np.pi/180 #set pitch to desired pitch
                Control_Vect.pop(0) #removes first time command to make the second command first in line
                if Echo:
                    print('Changed pitch to: ',theta_pitch*180/np.pi,' at ',np.round(t[i],2),' seconds.')
        except:
            pass
        
        

        for j in range(B): #Loop each blade
            
            THETA[j][i+1] = THETA[0][i] + omega[i]*dt +(j*2)*np.pi/B #position vector
            
            a23 = np.matrix([[np.cos(THETA[j][i]),np.sin(THETA[j][i]),0], #change of base matrix
                            [-np.sin(THETA[j][i]),np.cos(THETA[j][i]),0],
                            [0,0,1]])
            
            a14 = np.dot(np.dot(a34,a23),a12) #change from base 1 to 4
            a41 = np.transpose(a14) #other way arround
            
            for k in range(len(b_elm)-1): #Loop each element on blade
                
                r = (rt + np.matmul(a21,[0,0,-LS]) + np.matmul(a41,[b_elm[k],0,0])) #Position vector on point on blade
                x[i,j,k] = r[0,0] #save coordinates
                y[i,j,k] = r[0,1]
                z[i,j,k] = r[0,2]
                
                if Tower_influ and x[i,j,k]<= H:#cut off tower over hub heigh if option is activated
                    
                    a_t = 3.32 #tower radius
                else:
                    a_t = 0 #setting tower radius to 0 removes it's effects
                
                r_tower = np.sqrt(y[i,j,k]**2 + z[i,j,k]**2) #tower radius ref
                cos_theta = z[i,j,k]/r_tower 
                sin_theta = -y[i,j,k]/r_tower

                if Wind_shear: #apply wind shear effect if applicable
                    
                    U = V0*(x[i,j,k]/H)**v
                else:
                    U = V0
                
                Vr = cos_theta * U*(1-(a_t/r_tower)**2) #tower wind
                V_the = -U*(1+(a_t/r_tower)**2)*sin_theta
                
                #system 1
                V_x = 0
                V_z = Vr*cos_theta-V_the*sin_theta
                V_y = -Vr*sin_theta-V_the*cos_theta

                
                
                if Turbulence: #intrpolates the right turbulence to the right coordinate
                    Vturb_z[i,j,k] = fz(y[i,j,k]+180/2, x[i,j,k]+H)[0]
                    #Vturb_y[i,j,k] = fy(z[i,j,k], x[i,j,k])[0]

                Vo = np.matrix([[V_x],
                                [V_y - Vturb_y[i,j,k]],
                                [V_z + Vturb_z[i,j,k]]]) #wind matrix
                vvs[i,j,k] = V_z 
                
                #system 4

                V = np.dot(a14,Vo) #wind in other base
                Vy = V[1,0] 
                Vz = V[2,0] 
                
                

                Vrel_y =  Vy + Wy[i,j,k] - omega[i]*b_elm[k]*np.cos(theta_cone) # wind speeds 
                Vrel_z =  Vz + Wz[i,j,k] 
                Vrels = Vrel_y**2 + Vrel_z**2 #relativ wind

                flow_angle = np.arctan(Vrel_z/(-Vrel_y)) #find flow angle
                aot = flow_angle - (twist[k]*np.pi/180 + theta_pitch[i]) #angle of attack
                
                [Cl_inst,Cd,cmom,fs0,cl_inv,clfs] = force_coeffs_10MW_iced(aot*180/np.pi,t_c[k],aoa_tab,cl_tab,cd_tab,cm_tab,f_stat,cl_inv_tab,cl_fs) #interpolate cl cd
            
                if Dynamic_stall: #applies dynamic stall
                    tau = 4*c[k]/np.sqrt(Vrels)
                    fs[i+1,j,k] =  fs0 + (fs[i,j,k] - fs0)*np.exp(-dt/tau)
                    Cl[i+1,j,k] = fs[i+1,j,k]*cl_inv + (1-fs[i+1,j,k])*clfs
                else:
                    Cl[i+1,j,k] = Cl_inst
                    
                
                l = 0.5* rho * Vrels * c[k] * Cl[i+1,j,k] #lift and drag
                d = 0.5* rho * Vrels * c[k] * Cd
                
                
                pn[i+1,j,k] = l * np.cos(flow_angle) + d*np.sin(flow_angle) #elemet forces
                pt[i+1,j,k] = l * np.sin(flow_angle) - d*np.cos(flow_angle)
            

                a = -Wz[i,j,k]/V0 #axial induction
                
                A[i+1,j,k]=a
                
                if a <=1/3: #induction factor corrections
                    fg = 1
                #elif a >= 1/2:
                    #print('a >= 0.5')
                else:
                    fg = 1/4 * (5-3*a)
                
                
                Vo_Indu_wind = np.sqrt(Vy**2 + (Vz + fg*Wz[i,j,k])**2)
                inside_F = (-B/2)*(R-b_elm[k])/(b_elm[k]*np.sin(np.abs(flow_angle)))

                F = 2/np.pi * np.arccos(np.exp(inside_F)) #Prand't tip loss correction
                #ff[i+1,j,k] = F #debbugging 
                
                Wz_quasi[i+1,j,k] = -B*l*np.cos(flow_angle)/(4*np.pi*rho*b_elm[k]*F*Vo_Indu_wind) #quasy steady induced wind
                Wy_quasi[i+1,j,k] = -B*l*np.sin(flow_angle)/(4*np.pi*rho*b_elm[k]*F*Vo_Indu_wind)
                
                if not Dynamic_wake: #skip dynamic wake if time steps too low or not activated
                    Wz[i+1,j,k] = Wz_quasi[i+1,j,k]
                    Wy[i+1,j,k] = Wy_quasi[i+1,j,k]
                
                else: #dynamic wake
                    
                    tau_1 = 1.1/(1-1.3*min([a,0.5])) * R/V0
                    tau_2 = (0.39-0.26*(b_elm[k]/R)**2)*tau_1
                    #eq 7 forward differencing
                    
                    H_z = Wz_quasi[i+1,j,k] + k_D*tau_1 * (Wz_quasi[i+1,j,k] - Wz_quasi[i,j,k])/dt
                    H_y = Wy_quasi[i+1,j,k] + k_D*tau_1 * (Wy_quasi[i+1,j,k] - Wy_quasi[i,j,k])/dt
                    #eq 7 analytical solve
                    if i < 50: #if time step is too small skisp computation since convergence hasn't been achieved.
                        Wint_z[i+1,j,k] =  Wz_quasi[i+1,j,k]
                        Wint_y[i+1,j,k] =  Wy_quasi[i+1,j,k]
                        W_z[i+1,j,k] = Wint_z[i+1,j,k]
                        W_y[i+1,j,k] = Wint_y[i+1,j,k]
                    else:
                        Wint_z[i+1,j,k] = H_z + (Wint_z[i,j,k] - H_z)*np.exp(-dt/tau_1)
                        Wint_y[i+1,j,k] = H_y + (Wint_y[i,j,k] - H_y)*np.exp(-dt/tau_1)
                        #eq 8
                        W_z[i+1,j,k] = Wint_z[i+1,j,k] + (W_z[i,j,k] - Wint_z[i+1,j,k])*np.exp(-dt/tau_2)
                        W_y[i+1,j,k] = Wint_y[i+1,j,k] + (W_y[i,j,k] - Wint_y[i+1,j,k])*np.exp(-dt/tau_2)
                    Wz[i+1,j,k] = W_z[i+1,j,k]
                    Wy[i+1,j,k] = W_y[i+1,j,k]
                
                #add yaw effects here ?
        
        if Controller: #apply controller
            #Controller code
            GK = 1/(1+theta_pitch[i]/K1)
            theta_P[i+1] = GK*Kp*(omega[i] - omega_ref) #next theta P (proportional)
            theta_I[i+1] = theta_I[i] + GK*KI*(omega[i] - omega_ref)*dt #next theta I (integral)

            if theta_I[i+1] < theta_min:
                theta_I[i+1] = theta_min
            elif theta_I[i+1] > theta_max:
                theta_I[i+1] = theta_max

            theta_set_point = theta_P[i+1] + theta_I[i+1] #theta setp point
            theta_pitch[i+1] = theta_set_point #set next pitch to set point

            if theta_pitch[i+1] > theta_pitch[i] + theta_dot*dt:
                theta_pitch[i+1] = theta_pitch[i] + theta_dot*dt
            elif theta_pitch[i+1] < theta_pitch[i] - theta_dot*dt:
                theta_pitch[i+1] = theta_pitch[i] - theta_dot*dt
            
            if theta_pitch[i+1] > theta_max:
                theta_pitch[i+1] = theta_max
            elif theta_pitch[i+1] < theta_min:
                theta_pitch[i+1] = theta_min
            
            M_aero[i] = np.trapz(pt[i,0,:]*b_elm,b_elm)+np.trapz(pt[i,1,:]*b_elm,b_elm)+np.trapz(pt[i,2,:]*b_elm,b_elm) #Contributed effect of moment on all blade 

            # Pel = K * omega[i]**3
            # omega_rated = (Pel/K)**(1/3)
                
            if omega[i] <= omega_rated:
                M_g[i] = K*omega[i]**2
            else:
                M_g[i] = P_rated/omega_rated

            omega[i+1] = omega[i] + (M_aero[i] - M_g[i])/I_rotor *dt
            
            


    P = []
    T = []
    for i in range(len(t)):
        P.append(np.trapz(pt[i,0,:]*b_elm*omega[i],b_elm)+ np.trapz(pt[i,1,:]*b_elm*omega[i],b_elm)+ np.trapz(pt[i,2,:]*b_elm*omega[i],b_elm))
        T.append(np.trapz(pn[i,0,:],b_elm)+np.trapz(pn[i,1,:],b_elm)+np.trapz(pn[i,2,:],b_elm))

    # Initial proposed values (the rated ones)
    #P[0] = P_rated
    #T[0] = T_rated

    if Echo:
        print("\nRun time: ",np.round(time.time()-t0,2)," s\n")
        print("Models applied: ")
        if Wind_shear:
            print('\tWind shear')
        if Tower_influ:
            print('\tTower influence')
        if Dynamic_stall:
            print('\tDynamic stall')
        if Dynamic_wake:
            print('\tDynamic wake')
        if Turbulence:
            print('\tAtomospheric turbulence')
        if Controller:
            print('\tController')
        print("\n")
        

    output = locals()

    return output