"""Utility functions to run the Control Volume Heat c1b model

 Copyright (c) 2014, Gary Clow 

 Jan 2015.  Converted to Python: Brielle Kissack
 			Modified by Elchin Jafarov

 Functions:
 	feval(funcName, *args)
 	H2T(H, Cpave)
 	T2H(temp, Cpave)
 	K_eff(m,K,varep)
 	extract_str(s)
 	extract_num(s)
 	DE_coefs_1a(N,dtf,dt1mf,Dzrho,Kedz,cpave,QS,Jb,Jbn)
 	DE_coefs_1b(N,dtf,dt1mf,DzC,Kedz,QS,qb,qbn)
 	TDMA(N,aU,aP,aD,aUp,aPp,aDp,b,Hn,Hs)
 	init_T(N,t,spaceP,TpropP,QS,sourceP,bcP,testP,init_opt)
 	Z_grid(n,zfuL,zfdL,typL,nL,pos)
 	source_1p(z,sourceP)
 	source_2p(z,sourceP)
 	source_1(N,t,spaceP)
 	source_2(N,t,spaceP) not in there
 	set_experim()
 	Ksub_permafrost(T,phi,phi_i,phi_u,Kg25) Not sure complete or not
 	rho_permafrost(phi,phi_i,phi_u)
"""

import numpy as np
#import matplotlib.pyplot as plt

def feval(funcName, *args):  
    return eval(funcName)(*args)

def H2T(H, Cpave):
    
    """Converts enthalpy to temperature.
        
    The specific heat is allowed to be temperature dependent.
        
    Parameters
    ----------
    H : scalar  | Enthalpy.
    Cpave : ndarray  | Specific heat.
    
    Returns
    -------
    result : ndarray  | Temperature.
        
    Example
    --------
    >>> from c1b_util1 import H2T
    >>> H = 55.0; Cpave = np.array([0.5,1.0 ,0.75,0.8])
    >>> H2T(H,Cpave)
    Out[]: [ 110.        ,   55.        ,   73.33333333,   68.75      ]
        
    """
    
    Tref = 0  		# The reference temperature is assumed to be 0 C.
    temp = Tref + H/Cpave
    
    return temp

def T2H(temp, Cpave): 
      
    """Converts temperature to enthalpy. 
    
    The specific heat is allowed to be temperature dependent
   
    Parameters
    ----------
    T : scalar  | Temperature.
    Cpave : ndarray  | Specific heat.
    
    Returns
    -------
    result : ndarray  | Enthalpy.
    
    Example
    ------- 
    >>> from c1b_util1 import T2H   
    >>> T = 25.0; Cpave = np.array([0.5, 0.8,0.9,1.0])
    >>> T2H(T,Cpave) 
    Out[]: [ 12.5,  20. ,  22.5,  25. ]
   
    """ 
     
    Tref = 0            # The reference temperature is assumed to be 0 C.
    H = Cpave * (temp - Tref)  
 
    return H 
 
def K_eff(m,K,varep):
    
     """Find the effective conductivity at the interfaces.
     
     Parameters
     ----------
     m : scalar
     K : ndarray  | conductivity at the grid pts
     varep : ndarray  | fractional distance
     Returns
     -------
     result : ndarray  | Effective conductivity.
     
     Example
     -------
     >>> from c1b_util1 import K_eff
     >>> n = 4; a = np.ones(n); b = np.zeros(n) + 2 
     >>> K_eff(n,a,b)
     Out[]: array([ 1.,  1.,  1.,  0.])
     
     """
      
     Ke = np.zeros(m)
    
     for k in range(1,m):        
         Ke[k-1] = 1/((1-varep[k])/K[k-1] + varep[k]/K[k])
    
     return Ke
     
def extract_str(str):
    
     """Extract a string from a line that contains a comment.
     
     Parameters
     ----------
     s : string
     
     Returns
     -------
     result : string
     
     Example
     -------
     >>> from c1b_util1 import extract_str
     >>> str1 = "e3_layer % Lname" 
     >>> extract_str(str1)
     Out[]: 'e3_layer'
     
     """
     
     k = str.find('%')
     s2 = (str[0:k])
     
     return s2.strip()
     
def extract_num(str):
     
     """Extract a number from a line that contains a comment
     
     Parameters
     ----------
     s : string

     Returns
     -------
     result : float
     
     Example
     -------
     >>> from c1b_util1 import extract_num
     >>> str1 = "-10.5   % surface temperature"
     >>> extract_num(str1)
     Out[]: -10.5
     
     """

     k = str.find("%")
     x = float(str[0:k])
      
     return x
     
def DE_coefs_1a(N,dtf,dt1mf,Dzrho,Kedz,cpave,QS,Jb,Jbn):
    
     """Need to be careful with this function!!!
     
     Find DE coefs for CV code cv1a.m
     
     diffusion:               yes
     advection:               no
     composite materials:     yes
     melting:                 no
     source term:             yes
     
     Note:  for unconditional stability, aPp should be greater than 0.
     This can fail if step dt is too large.
     
     Parameters
     ----------
     N : integer    # Note, N is global parameter in the matlab code
     dtf : integer
     dt1mf : integer
     Dzrho : ndarray
     Kedz : ndarray
     cpave : ndarray
     QS : ndarray
     Jb : integer
     Jbn : integer
     
     Returns
     -------
     Result : ndarray
     
     Example
     -------
     >>> from c1b_util1 import DE_coefs_1a
     >>> N = 4; dtf = 1; dt1mf = 3; Jb = 7; Jbn = 8
     >>> Dzrho = np.zeros(N) + 4
     >>> Kedz = np.zeros(N) + 3
     >>> cpave = np.zeros(N) + 5
     >>> QS = np.zeros(N) + 6
     >>> DE_coefs_1a(dtf,dt1mf,Dzrho,Kedz,cpave,QS,Jb,Jbn)
     
     (array([ 0. ,  0.6,  0.6,  0.6]),
 	 array([ 0. ,  5.2,  5.2,  4.6]),
	 array([ 0. ,  0.6,  0.6,  0. ]),
	 array([ 0. ,  1.8,  1.8,  1.8]),
	 array([ 0. ,  0.4,  0.4,  2.2]),
	 array([ 0. ,  1.8,  1.8,  0. ]),
	 array([  0.,   6.,   6., -25.]))
     
     """

     aU = np.zeros(N)
     aP = np.zeros(N)
     aD = np.zeros(N)
     aUp = np.zeros(N)
     aPp = np.zeros(N)
     aDp = np.zeros(N)
     b = np.zeros(N)
     
     for k in range(1,N-1):
         
         aU[k] = dtf * Kedz[k] / cpave[k-1]
         aUp[k] = dt1mf * Kedz[k] / cpave[k-1]
         aD[k] = dtf * Kedz[k+1] / cpave[k+1]
         aDp[k] = dt1mf * Kedz[k+1] / cpave[k+1]
         aP[k] = Dzrho[k] + dtf * (Kedz[k]/cpave[k] + Kedz[k+1]/cpave[k])
         aPp[k] = Dzrho[k] - dt1mf * (Kedz[k]/cpave[k] + Kedz[k+1]/cpave[k])
         b[k] = QS[k]
         
     k = N - 1
     aU[k] = dtf * Kedz[k] / cpave[k-1]
     aUp[k] = dt1mf * Kedz[k] / cpave[k-1]
     aD[k] = 0
     aDp[k] = 0
     aP[k] = Dzrho[k] + dtf * Kedz[k]/cpave[k]
     aPp[k] = Dzrho[k] - dt1mf * Kedz[k]/cpave[k]
     b[k] = QS[k] - dtf*Jb - dt1mf*Jbn
     
     return aU,aP,aD,aUp,aPp,aDp,b
     
def DE_coefs_1b(N,dtf,dt1mf,DzC,Kedz,QS,qb,qbn):
    
     """Find DE coefs for CV code cv1b.m
     
     diffusion:            yes
     advection:            no
     composite materials:  yes
     melting:              yes
     source term:          yes
     
     Note:  for unconditional stability, aPp should be greater than 0.
     This will fail if time step dt is too large.
     
     Parameters
     ----------
     dtf : integer 
     dt1mf : integer
     Dzc : ndarray 
     Kedz : ndarray
     QS : ndarray
     qb : integer
     qbn : integer
     Returns
     -------
     result : ndarrays
     
     Example
     -------
     >>> from c1b_util1 import DE_coefs_1b
     >>> N = 4; dtf = 3; dt1mf = 2; qb = 5; qbn = 8
     >>> DzC = np.zeros(N) + 4
     >>> Kedz = np.zeros(N) 
     >>> QS = np.zeros(N) + 1
     >>> DE_coefs_1b(N,dtf,dt1mf,DzC,Kedz,QS,qb,qbn)
     Out[]: 
	(array([ 0.,  0.,  0.,  0.]),
 	array([ 0.,  4.,  4.,  4.]),
 	array([ 0.,  0.,  0.,  0.]),
 	array([ 0.,  0.,  0.,  0.]),
 	array([ 0.,  4.,  4.,  4.]),
 	array([ 0.,  0.,  0.,  0.]),
 	array([  0.,   1.,   1., -30.]))
     """
     
     aU = np.zeros(N)
     aP = np.zeros(N)
     aD = np.zeros(N)
     aUp = np.zeros(N)
     aPp = np.zeros(N)
     aDp = np.zeros(N)
     b = np.zeros(N)
     
     for k in range(1,N-1):
         aU[k] = dtf * Kedz[k]
         aUp[k] = dt1mf * Kedz[k]
         aD[k] = dtf * Kedz[k]
         aDp[k] = dt1mf * Kedz[k+1]
         aP[k] = DzC[k] + dtf * (Kedz[k] + Kedz[k+1])
         aPp[k] = DzC[k] - dt1mf * (Kedz[k] + Kedz[k+1])
         b[k] = QS[k] 
         k = N-1
     
     aU[k] = dtf * Kedz[k]
     aUp[k] = dt1mf * Kedz[k]
     aD[k] = 0
     aDp[k] = 0
     aP[k] = DzC[k] + dtf * Kedz[k]
     aPp[k] = DzC[k] - dt1mf * Kedz[k]
     b[k] = QS[k] - dtf*qb - dt1mf*qbn
     
     return aU,aP,aD,aUp,aPp,aDp,b
    
def TDMA(N,aU,aP,aD,aUp,aPp,aDp,b,Hn,Hs):
     
     """TDMA algorithm in the z-direction for either enthalpy-based or 
     temperature-based CV models.
     
     Sweep from top to bottom.
     
     Parameters:
     -----------
     aU : ndarray
     aP : ndarray
     aD : ndarray
     aUp : ndarray
     aPp : ndarray
     aDp : ndarray
     b : ndarray
     Hn : ndarray
     Hs : integer
     
     Returns
     -------
     result : ndarray
     
     Example
     -------
     >>> from c1b_util1 import TDMA
     >>> n = 4
     >>> aU = np.ones(n) + 2
     >>> aP = np.ones(n) 
     >>> aD = np.ones(n) + 1
     >>> aUp = np.ones(n) + 3 
     >>> aPp = np.ones(n) + 4 
     >>> aDp = np.ones(n) + 5 
     >>> b = np.ones(n) + 6 
     >>> Hn = np.ones(n) + 7 
     >>> Hs = 1
     >>> TDMA(N,aU,aP,aD,aUp,aPp,aDp,b,Hn,Hs)
     Out[]: array([   1.        ,    7.27272727,  -61.36363636, -105.09090909])
     """
     
     E = np.zeros(N)
     P = np.zeros(N)
     Q = np.zeros(N)
     H = np.zeros(N)
     
     for k in range(1,N-1):
         E[k] = aUp[k] * Hn[k-1] + aPp[k] * Hn[k] + aDp[k] * Hn[k+1] + b[k] 
     E[N-1] = aUp[N-1] * Hn[N-2] + aPp[N-1] * Hn[N-1] + b[N-1]
     
     H[0] = Hs
     
     P[N-1] = aU[N-1]/aP[N-1]
     Q[N-1] = E[N-1]/aP[N-1]
     
     for k in range(N-2,0,-1):
         fac = aP[k] - aD[k] * P[k+1]
         P[k] = aU[k]/fac
         Q[k] = (E[k] + aD[k] * Q[k+1])/fac
         
     for k in range(1,N):
         H[k] = P[k] * H[k-1] + Q[k]
         
     return H
     
def init_T(N,t,spaceP,TpropP,QS,sourceP,bcP,testP,init_opt):
	     
	"""Setup initial temperature field using the following methods,

	init_opt = 1 input initial field from a file
	init_opt = 2 calcuate field assuming steady-state conditions
	init_opt = 3 use analytic expressions

	Parameters
	----------
	N : integer   
	t : integer
	spaceP : list
	TpropP : list
	QS : ndarray
	sourceP : ndarray
	bcP : ndarray
	testP : list
	init_opt : integer
	Returns
	-------
	Result : ndarray 

	Example
	-------
	>>> from c1b_util1 import init_T
	>>> N = 4; t = 5
	>>> spaceP = [1,range(0,N),np.zeros(N) + 0.2,np.zeros(N) + 0.3]
	>>> sourceP = np.zeros(N)+4
	>>> TpropP = [np.zeros(N)+0.1,2.4,1e6,np.zeros(N)+0.5]
	>>> QS = np.zeros(N) + 3
	>>> bcP = np.zeros(N) + 5
	>>> testP = ['e2',0,0.1,2]
	>>> init_opt = 3
	>>> z = spaceP[1]
	>>> init_T(N,t,spaceP,TpropP,QS,sourceP,bcP,testP,init_opt)
	Out[]: array([ 141.56749883,  251.82037778,  337.68540625,  404.55715765])

	"""
	T = np.zeros(N); secyr=86400*365

	z = spaceP[1]; Dz = spaceP[2]; dz = spaceP[3]; 
	#NOTE: we cannot do zn = z[0:N], bc we end with N-1 element
	zn = z # grid points without the one embedded in lowest interface

	K = TpropP[0]; rho = TpropP[1]; cp = TpropP[2]; Ke = TpropP[3]
	diffu = K/(rho*cp)

	Ts = bcP[0]; Jb = bcP[1]; qb = -Jb

	if init_opt == 1:
		print ''
		Tzname = raw_input('Enter a filename: ') or 'default_T_init.txt'
		try:
			with open(Tzname) as file:
				pass
		except IOError as e:
			print 'Unable to open', Tzname, 'file'
		# First 3 rows of Tzname file has to be headers
		wk = np.loadtxt(Tzname, skiprows=3,unpack=False)
		zA = wk[:,0]   # depth (m)
		TzA = wk[:,1]  # T(z)

	elif init_opt == 2:
		J = np.zeros(N+1)
		J[N] = -qb
		print N+1, J
		for i in range(N-1,-1,-1):
			J[i]=J[i+1] - QS[i] # NOTE: there is a bug a Gary's code it says QS is an empty matrix
			print J[i]
 
		Kedz = Ke/dz # if this is zero then we are going to encounter div by 0 later
		T[0] = Ts
		print Ke, Kedz
		for i in range(1,N):
			T[i] = T[i-1] - J[i]/Kedz[i]

	elif init_opt == 3:

		experim = testP[0]
		T0 = testP[1]
		dT = testP[2]
		period = testP[3]
		print experim

		if experim == 'e1':
			print T0, zn, qb,K
			T = T0 + zn * qb / K 

		elif experim == 'e2':
			#S0 = sourceP[0]; h = sourceP[1]
			#T = T0 + (S0/K) * (h**2) * (1- np.exp(-zn/h))
			print 'This part does not work, bc of sourceP = [0], has to have more than 1 element'

		elif experim == 'e3':
			Tf = np.zeros(N+1)
			Tf[0] = T0
			for k in range(0,N):
				Tf[k+1] = Tf[k] + Dz[k] * qb/K[k]
			T[0] - Tf[0]
			for k in range(1,N):
				T[k] = 0.5*(Tf[k] + Tf[k+1])

		elif experim == 'e4':
			K0 = 2
			dKdT = -0.0140
			T[0] = T0
			for k in range(1,N):
				C = -(2*K0*T[k-1] + dKdT*T[k-1]**2 + 2*qb*dz[k])
				T[k] = (-K0 + np.sqrt(K0**2 - dKdT*C))/dKdT

		elif experim == 'e5':
			T = T0 + zn * (qb/K)

		elif experim == 'e6':
			lambda_ = secyr*period # period =0, getting div 0 on the next line
			w = 2*np.pi/lambda_
			k = np.sqrt(np.pi/(lambda_*diffu))
			T = T0 + dT * np.cos(k*zn - w*t + np.pi/2) * np.exp(-k*zn)

		elif experim == 'e7':
			T = T0 + zn * qb / K

		else:
			import sys
			sys.exit("init_T: sorry no analytic soln available when itest = 0")    

	return T  

def Z_grid1(zfuL, zfdL, typL, nL, pos):   

# Routine setups up the Control Volume Space Grid (Z-direction).

# The location and type of physical layers are input.  The space grid uses
# the boundaries between the physical layers for CV interfaces.  In addition,
# each physical layer can be broken into multiple control volumes (eg. 25 CVs)
#may be used to represent a 50-m thick sandstone layer).

# Physical Layers:

# Control Volume Space Grid:

# N = number of control volumes
# zf = CV interfaces                       (N+1)
# z = CV grid points                       (N+1) with one embedded in lowest interface
# Dz = width of CVs                        (N)
# dz = distance between grid points        (N)
# varep = fractional distance              (N)
# typ = material type @ each grid point    (N)

# Note: The upper physical layer contains a 1/2 CV, thus nL(1) must be a number
# like 1.5, 2.5, 4.5, etc.

# ----------------------------------------------------------------------------     
    
    n = 4 
    nlayers = len(typL)
# Test nL(1), see note above.

    nL = np.ones(n) + 1.5  
    
    if nL[0]%1 != 0.5: 
        print(' ') 
        print('Error: Z_grid.m') 
        print('Upper physical layer needs to include a 1/2-CV.  Thus, nL(1) ') 
        print('should be a number like 1.5, 2.5, 4.5, etc.') 
        #pause

# Define CV interface locations.
# Also define material type at each grid point.
    
    j = 1 
    Dz = (zfdL[0] - zfuL[0]) / nL[0]                # CV width for 1st physical layer 
    zf = [zfuL[0], zfuL[0] + 0.5*Dz]                # first CV 
    zfx = np.arange(zf[1] + Dz, zfdL[0]+Dz, Dz)     # remaining CVs in this phy. layer 
    zf = zf = np.concatenate((zf, zfx))           
    n = len(zf) - 1                                 # number of CVs so far                                                     
    for i in range(1,n): 
        typ[i] = typL[0] 
    nn = n 
    for j in range(1, nlayers): 
        Dz = (zfdL[j] - zfuL[j]) / nL[j]            # CV width for remaining layers 
        zfx = np.arange(zfuL[j] + Dz, zfdL[j], Dz) 
        zf = np.concatenate((zf,zfx)) 
        n = len(zf) - 1                            
        for i in (nn, n-1): 
            typ[i] = typL[j] 
        nn = n 
    N = len(zf) - 1                                 # total number of CVs                                                         
    typ[0] = typ[1]                                 # material type @ upper bnd 
    del Dz

# Define u & d interfaces for each of the N control volumes.  
    
    zfu = np.zeros(N) 
    zfd = np.zeros(N) 
    for i in range(0,N): 
        zfu[i] = zf[i] 
        zfd[i] = zf[i+1]
        
# Define grid point location at the center of each CV.
# Also embed a grid point within the upper & lower boundaries.
        
        
    z = np.zeros(N+1) 
    for i in range(1,N): 
        z[i] = 0.5*(zfu[i] + zfd[i]) 
    z[0] = zf[0] 
    z[N] = zf[N]
    
# Define width of control volumes.    
    
    Dz = np.zeros(N) 
    for i in range(0,N): 
        Dz[i] = zfd[i] - zfu[i]
        
# Define distance between grid points.        
        
    dz = np.zeros(N) 
    for i in range(1,N): 
        dz[i] = z[i] - z[i-1] 
    dz[0] = 10E-10
    
# Fractional distance of interface i to grid pt (i-1), from grid pt i.    
    
    varep = np.zeros(N) 
    for i in range(1,N): 
        varep[i] = (z[i] - zf[i])/ dz[i]

# > Display setup                
    
    return zf,z,Dz,dz,varep,typ,N

             
def Z_grid(n,zfuL,zfdL,typL,nL,pos):
     
     """Routine setups up the Contour Volume Space Grid (z-direction).
     Note: this subroutine is incomplete
     
     Parameters
     ----------
     zfuL : ndarray
     zfdL : ndarray
     typL : ndarray
     nL : ndarray
     pos : integer
     Returns
     -------
     Result : ndarray
     
     Example
     -------
     >>> from c1b_util1 import Z_grid 
     >>> n = 4; pos = 1
     >>> nL = np.ones(n) + 1.5
     >>> typL = np.ones(n) + 2
     >>> nlayers = len(typL)
     >>> zfdL = np.ones(n) + 1
     >>> zfuL = np.ones(n)
     >>> Z_grid(zfuL, zfdL, typL, nL, pos)
     Out[23]: 
	 (array([ 1. ,  1.2,  1.6,  2. ,  1.4,  1.8,  1.4,  1.8,  1.4,  1.8]),
 	 array([ 1. ,  1.4,  1.8,  1.7,  1.6,  1.6,  1.6,  1.6,  1.6,  1.8]),
	 array([ 0.2,  0.4,  0.4, -0.6,  0.4, -0.4,  0.4, -0.4,  0.4]),
	 array([  1.00000000e-09,   4.00000000e-01,   4.00000000e-01,
         -1.00000000e-01,  -1.00000000e-01,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00]),
	 array([ 0. ,  0.5,  0.5,  3. , -2. , -inf,  inf, -inf,  inf]),
	 array([ 3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  0.]),
	 9)
     """

     typ = np.zeros(10)
     nlayers = len(typL)

     nL = np.ones(n) + 1.5  
    
     if nL[0]%1 != 0.5:
         print(' ') 
         print('Error: Z_grid.m') 
         print('Upper physical layer needs to include a 1/2-CV.  Thus, nL(1) ') 
         print('should be a number like 1.5, 2.5, 4.5, etc.') 
    
     j = 1 
     Dz = (zfdL[0] - zfuL[0]) / nL[0]                 
     zf = [zfuL[0], zfuL[0] + 0.5*Dz]                 
     zfx = np.arange(zf[1] + Dz, zfdL[0]+Dz, Dz)      
     zf = zf = np.concatenate((zf, zfx))           
     n = len(zf) - 1                                                                                      
     for i in range(1,n): 
         typ[i] = typL[0] 
     nn = n 
     for j in range(1, nlayers): 
         Dz = (zfdL[j] - zfuL[j]) / nL[j]             
         zfx = np.arange(zfuL[j] + Dz, zfdL[j], Dz) 
         zf = np.concatenate((zf,zfx)) 
         n = len(zf) - 1                            
         for i in (nn, n-1): 
             typ[i] = typL[j] 
         nn = n 
     N = len(zf) - 1                                                                                          
     typ[0] = typ[1]                                  
     del Dz  
    
     zfu = np.zeros(N) 
     zfd = np.zeros(N) 
     for i in range(0,N): 
         zfu[i] = zf[i] 
         zfd[i] = zf[i+1]
        
     z = np.zeros(N+1) 
     for i in range(1,N): 
         z[i] = 0.5*(zfu[i] + zfd[i]) 
     z[0] = zf[0] 
     z[N] = zf[N]
    
     Dz = np.zeros(N) 
     for i in range(0,N): 
         Dz[i] = zfd[i] - zfu[i]
                
     dz = np.zeros(N) 
     for i in range(1,N): 
         dz[i] = z[i] - z[i-1] 
     dz[0] = 10E-10    
    
     varep = np.zeros(N) 
     for i in range(1,N): 
         varep[i] = (z[i] - zf[i])/ dz[i]                
    
     return zf,z,Dz,dz,varep,typ,N # Stops here ???
     
     Dz = (zfdL[0] - zfuL[0]) / nL[0]
     zf = [zfuL[0], zfuL[0] + 0.5*Dz]
     zfx = np.arange(zf[1] + Dz, zfdL[0]+Dz, Dz)
     zf = np.concatenate((zf, zfx))
     n = len(zf) - 1
     for i in range(1,n):
         typ[i] = typL[0]
     nn = n
     for j in range(1,nlayers):
         Dz = (zfdL[j] - zfuL[j]) / nL[j]
         zfx = np.arange(zfuL[j] + Dz, zfdL[j], Dz)
         zf = np.concatenate((zf,zfx))
         n = len(zf) - 1
         for i in (nn, n-1):
             typ[i] = typL[j]
         nn = n
     N = len(zf) - 1
     typ[0] = typ[1]
     del Dz
     zfu = np.zeros(N)
     zfd = np.zeros(N)
     for i in range(0,N):
         zfu[i] = zf[i]
         zfd[i] = zf[i+1]
     z = np.zeros(N+1)
     for i in range(1,N):
         z[i] = 0.5*(zfu[i] + zfd[i])
     z[0] = zf[0]
     z[N] = zf[N]
     Dz = np.zeros(N)
     for i in range(0,N):
         Dz[i] = zfd[i] - zfu[i]
     dz = np.zeros(N)
     for i in range(1,N):
         dz[i] = z[i] - z[i-1]
     dz[0] = 10E-10
     varep = np.zeros(N)
     for i in range(1,N):
         varep[i] = (z[i] - zf[i])/ dz[i]

     junkx = float('NaN') *np.ones(zfuL.shape)
     plt.plot(junkx,zfuL)
     plt.plot(junkx,zfdL)
     v = plt.axis([0,1,0,2])
     plt.gca().invert_yaxis()
     junk2 = 0.5*(v[1] - v[0] * np.ones(z.shape))
     plt.ylabel('Depth (m)')
     plt.title('Physical Layers (blue), CV interfaces (red), CV grid pts (black)')
# show physical layers
     for j in range(0,nlayers):
         plt.axhline(zfuL[j],0,1)
     plt.axhline(zfdL[j],0,1)
# show CV interfaces
     for i in range(0,N+1):
         plt.axhline(zf[i],0,1,color='r')
# show CV grid pts
     for i in range(0,N+1):
         plt.plot(junk2,z,'o',markerfacecolor = 'k')
     plt.axis([0,1,0.9,2.1])
     plt.gca().invert_yaxis()
     plt.show()

	
def source_1p(z,sourceP):
    
    """Source term at each grid point.
     
    Parameters
    ----------
    sourceP : ndarray
    Z : scalar
    Returns
    -------
	Result : ndarray
     
    Example
    -------
    >>> n = 1.0; sourceP = np.ones(n)*0.5
	>>> source_1p(1,sourceP)
    Out[]: 0.5
    """
    
    S0 = sourceP[0]  		
    S = S0*z
    
    return S

def source_2p(z,sourceP):
    
    """Source term at each grid point.
     
    Parameters
    ----------
    sourceP : ndarray
    Z : scalar
    Returns
    -------
	Result : ndarray
     
    Example
    -------
    >>> n = 1.0; sourceP = np.ones(n)*0.5
	>>> source_2p(1.5,sourceP)
    Out[]: 
    """
    
    S0 = sourceP[0]  
    h = sourceP[1]		
    S = S0*np.exp(-z/h)
    
    return S

def source_1(N,t,spaceP):
    
    """This function evaluates the source term during the setup stage.
    It:
	> Finds the integral of the source term across each CV (QS)
	> Finds the ratio Sfac = QS/[S(z(i))*Dz(i)] where z(i) = grid pt for CV(i)
	> Returns values for source parameters in cell array sourceP.

	For each CV,

	QS = \int S(z) dz = <S>*Dz	[W/m^2]

	where <S> is the mean value of S across the CV.  The purpose of Sfac is
	to allow us to estimate QS during the main time loop without having to do
	an integration, eg.

	QS(i) = <S(i)> * Dz(i) = [Sfac * S(z(i))] * Dz(i) .

	This assumes the shape of S(z) doesn't change with time.

	source version '1' (test).
     
    Parameters
    ----------
    spaceP : ndarray
    t : scalar
    Returns
    -------
	Result : ndarray
     
    Example
    -------
    >>> n = 4
	>>> spaceP = np.ones(n)*0.5
	>>> source_1(n,1,spaceP)
    Out[]: (array([ 0.,  0.,  0.,  0.]), array([ 1.,  1.,  1.,  1.]), 0)
    Note: Ask Gary why he passes t, since it is never used in by the function
    """
    
    S0 = 0  				#Set Source parameters
    sourceP = [S0]			#Note: Gary has {} sign for {S0} (cell array)
    						#	 	in python []
    QS = np.zeros(N)  		#Set integral
    Sfac = np.ones(N)	 	#Set Sfrac
    
    return QS,Sfac,sourceP
    

def set_experim():
    
    """This function sets an experiment    
    """
    
    regions = ['test','alaska']
    print 'Regions:'
    print '[0] test'
    print '[1] alaska'

    ireg = input('Type region#: ')
    region = regions[ireg]
    print region
    print '-' * 20

    print('Projects: ')

    if region == 'test':
        projects = ['easy','rock','ice']
        print '[0] easy'
        print '[1] rock'
        print '[2] ice'

    elif region == 'alaska':
        projects = ['active-layer','deep-permafrost']
        print '[0] active-layer'
        print '[1] deep-permafrost'

    iproj = input('Type project#: ')
    project = projects[iproj]
    print project
    print '-' * 10

    print('Experiments =')

    if project == 'easy':
        expers = ['e1','e2','e3','e4','e5','e6','e7']
        print '[0] e1 SS, UB = fixed Ts, LB = fixed Jb, source=0'
        print '[1] e2 SS with exponential source term'
        print '[2] e3 SS with composite material'
        print '[3] e4 SS with temp.-dependent K and Cp'
        print '[4] e5 instantaneous step change'
        print '[5] e6 periodic surface temp'
        print '[6] e7 triangle heating pulse'
        iexper = input('Type experiment#: ')
        experim = expers[iexper]

    elif project == 'active-layer':
        expers = ['dp1','dp2','dp3','dp4']
        print '[0] dp1 SS, UB = fixed Ts, LB = fixed Jb, source = 0'
        print '[1] dp2 triangle heating pulse'
        print '[2] dp3 instantaneous step change'
        print '[3] dp4 nominal Chukchi coast warming'
        iexper = input('Type experiment#: ')
        experim = expers[iexper]

    elif project == 'deep-permafrost':
        expers = ['dp1','dp2','dp3','dp4']
        print '[0] dp1 SS, UB = fixed Ts, LB = fixed Jb, source = 0'
        print '[1] dp2 triangle heating pulse'
        print '[2] dp3 instantaneous step change'
        print '[3] dp4 nominal Chukchi coast warming'
        iexper = input('Type experiment#: ')
        experim = expers[iexper]

    else:
        print('not defined')
        experim = ''

    experiment = (region,project,experim)
    print experiment
    
    return experiment

def Ksub_permafrost(T,phi,phi_i,phi_u,Kg25):

    """NOTE: Maybe not finished.
     
    Parameters
    ----------
    T : ndarray
    phi : scalar
    phi_i : scalar
    phi_u : scalar
    Kg25 : scalar
    Returns
    -------
	Result : ndarray
     
    Example
    -------
    >>> T = np.array([10.0,15.0,12.0,14.0])
	>>> phi = 0.025; phi_i = 0.3; phi_u = 0.15; Kg25 = 1.5
	>>> Ksub_permafrost(T,phi,phi_i,phi_u,Kg25)
    Out[]: array([ 0.57878726,  0.5888032 ,  0.58290002,  0.58686936])
    """

    dflag = 0
    Tk = T + 273.15
    nT = len(T)
    f = 0.25e-02
    fac = (1 - (T-25)*f)
    Kg = fac*Kg25
    a = 9.828
    b = 0.0057
    Ki = a * np.exp(-b*Tk)
    c = np.array([1.6630,-1.7781,1.1567,-0.432115])
    d = np.array([-1.15,-3.4,-6.0,-7.6])
    Tr = Tk/300
    Ku = np.zeros(len(T))
    
    for j in range(0,nT):
        sum1 = 0 
        for i in range(4):
            sum1 = sum1 + c[i] * Tr[j]**d[i]
        Ku[j] = sum1
        
    return Ku # looks like it should return Ku not K
    
def rho_permafrost(phi,phi_i,phi_u):

    """Density of permafrost?
     
    Parameters
    ----------
    phi : scalar
    phi_i : scalar
    phi_u : scalar
    Returns
    -------
	Result : ndarray
     
    Example
    -------
    >>> phi = 1.125; phi_i = 0.25; phi_u = 0.225
	>>> rho_permafrost(phi,phi_i,phi_u)
    Out[]: 123.0
    """

    rhog = 2650
    rhoi = 917
    rhow = 1000
    rho = (1-phi)*rhog + phi_i*rhoi + phi_u*rhow
    
    return rho
    
def Tprop_1(dummy_var1, materialC, dummy_var2):


    """Density of permafrost?
     
    Parameters
    ----------
    phi : scalar
    phi_i : scalar
    phi_u : scalar
    Returns
    -------
	Result : ndarray
     
    Example
    -------
    >>> phi = 1.125; phi_i = 0.25; phi_u = 0.225
	>>> rho_permafrost(phi,phi_i,phi_u)
    Out[]: 123.0
    """    
    
    typ = materialC[0]
    L=typ[0]
    print 'idx=',L

# Type 1 material -------
    if L == 1:
        rho = 2000*np.ones(typ.shape)
        Cp = 1000*np.ones(typ.shape)
        K = 2*np.ones(typ.shape)

# Type 2 material -------
    if L == 2:
        rho = 2500*np.ones(typ.shape)
        Cp = 1600*np.ones(typ.shape)
        K = 4*np.ones(typ.shape)

# bulk heat capacity
    C = rho * Cp

# temperature weighted specific heat
    Cpave = Cp

# store Cpave in diagnostics cell array
    diagC = [Cpave]

    return K,rho,Cp,C,diagC

