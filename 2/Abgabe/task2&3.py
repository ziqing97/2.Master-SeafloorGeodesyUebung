#%%
# pip install seawater
# https://pythonhosted.org/seawater/eos80.html

import seawater as sw
from seawater.library import T90conv
import numpy as np
import matplotlib.pyplot as plt

# seawater.eos80.salt(r, t, p)
#    Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.
#    Parameters:	
#    r : array_like
#        conductivity ratio R = \frac{C(S,T,P)}{C(35,15(IPTS-68),0)}
#    t : array_like
#        temperature 
#    C (ITS-90)]
#    p : array_like
#        pressure [db]
#    Returns:	
#    s : array_like
#        salinity [psu (PSS-78)]
#    References
#    Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for computation of fundamental properties of seawater. UNESCO Tech. Pap. in Mar. Sci., No. 44, 53 pp. Eqn.(31) p.39. http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

# p1 and p2 in dbar!  
# Conversion 1 KPa=0.1 dbar 
# Conversion done with
# https://www.asknumbers.com/PressureConversion.aspx
 
# Temperature deg. Celcius , Pressure in dbar  = 0.1 KPa , Salinity in PSU 



### Enter here the values which you expect for your deployment 

# Marmara 
t1 	= 14  			# In deg.
sal1 	= 38.6  		# In PSU
p1kPa 	= 8195.0  		# In kPa	
p1    	= ( p1kPa/float(10) )# In dbar



### Calculate changes in T°, P and SAL

# Change of Parameters used for calculation
deltat 		= 1 
deltas 		= 1 
deltapdkPa 	= 1
deltap 		= deltapdkPa / float(10)  

# Tempereture
t2 	= t1 + deltat  
sal2 	= sal1 
p2 	= p1  
p2kPa   = p2 * 10  


t = T90conv([  t1 ])
s = [ sal1  ]
p = [ p1 ]
sv1 = sw.svel(s, t, p)

t = T90conv([  t2 ])
s = [ sal2  ]
p = [ p2 ]
sv2 = sw.svel(s, t, p)

gradt= (t1-t2)/float((sv1-sv2))

# In kPa  --> Divide by 10 for dbar which is used by seawater 

# Pressure
t2 	= t1   
sal2 	= sal1 
p2 	= p1 + deltap  
p2kPa   = p2 * 10  

t = T90conv([  t1 ])
s = [ sal1  ]
p = [ p1 ]
sv1 = sw.svel(s, t, p)

t = T90conv([  t2 ])
s = [ sal2  ]
p = [ p2 ]
sv2 = sw.svel(s, t, p)

gradpdbar = (p1-p2)      /float((sv1-sv2))   
gradpkPa  = ((p1-p2)*10)/float((sv1-sv2))

# Salinity
t2 	= t1   
sal2 	= sal1 + deltas  
p2 	= p1   
p2kPa   = p2 * 10  

t = T90conv([  t1 ])
s = [ sal1  ]
p = [ p1 ]
sv1 = sw.svel(s, t, p)

t = T90conv([  t2 ])
s = [ sal2  ]
p = [ p2 ]
sv2 = sw.svel(s, t, p)

grads = (sal1-sal2)/float((sv1-sv2))

# print gradient 
print("dvdt" + "%8.3f" %(1/float(gradt)))
print("dvds" +  "%8.3f" %(1/float(grads)))
print("dvdp" +  "%8.3f" %(1/float(gradpkPa)))

# velocity 
t = T90conv([  t1 ])
s = [ sal1  ]
p = [ p1 ]
sv = sw.svel(s, t, p)
print("v= " + "%8.3f" %(sv) + "\n")

BSLlength = 1000

# for 1cm Baseline change
dl = 0.01
critical_dt = dl * float(gradt) * sv / BSLlength
critical_ds = dl * float(grads) * sv / BSLlength
critical_dp = dl * float(gradpkPa) * sv / BSLlength

print("Task 3: The parameter change for 1cm baseline change\n")
print("t: " + "%8.3f" %(critical_dt) + " deg")
print("sal: " + "%8.3f" %(critical_ds) + " PSU")
print("p: " + "%8.3f" %(critical_dp) + " kPa")


# %% Task 2 
# data
t 	= 14  			# In deg.
sal 	= 38.6  		# In PSU
pkPa 	= 8195.0  		# In kPa	
p    	= ( pkPa/float(10) )# In dbar

deltat 		= 1 
deltas 		= 0.1 
deltapdkPa 	= 1
deltap 		= deltapdkPa / float(10)  

# change in temperature
t1 = np.arange(0, 100, deltat)
t1 = T90conv(t1)
sal1 = sal
p1 = p

sv1 = sw.svel(sal1, t1, p1)
plt.plot(t1,sv1)
plt.xlabel("temperature (degree in ITS90)")
plt.ylabel("soudn velocity m/s")
plt.savefig("./tempereture")
plt.show()

# change in PSU 
sal2 = np.arange(28.6, 48.6, deltas)
t2 = T90conv(t)
p2 = p 
sv2 = sw.svel(sal2, t2, p2)
plt.plot(sal2,sv2)
plt.xlabel("salinity in PSU")
plt.ylabel("soudn velocity m/s")
plt.savefig("./salinity")
plt.show()

# change in pressure
sal3 = sal
t3 = T90conv(t)

p3 = np.arange(8095, 8295, deltap) 
sv3 = sw.svel(sal3, t3, p3)
plt.plot(p3,sv3)
plt.xlabel("pressure in dbar")
plt.ylabel("soudn velocity m/s")
plt.savefig("./pressure")
plt.show()


