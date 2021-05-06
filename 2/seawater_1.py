# pip install seawater
# https://pythonhosted.org/seawater/eos80.html

import seawater as sw
from seawater.library import T90conv

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

print ("########################")
print ("T Difference deg. C  ", deltat )
print ("########################")

t2 	= t1 + deltat  
sal2 	= sal1 
p2 	= p1  
p2kPa   = p2 * 10  

print ("t1,t2, sal1,sal2,  p1,p2 (dbar), p1,p2 (kPa): ")
print (t1,t2, sal1,sal2, p1,p2, p1kPa , p2kPa)

t = T90conv([  t1 ])
s = [ sal1  ]
p = [ p1 ]
sv1 = sw.svel(s, t, p)

t = T90conv([  t2 ])
s = [ sal2  ]
p = [ p2 ]
sv2 = sw.svel(s, t, p)

print ("sv1,sv2:")
print (sv1,sv2)
print ("Difference in SV") 
print (sv1-sv2)
#print "Difference T"
#print t1,t2,(t1-t2) 
print ("Gradient in deg. / (m/s)" )
gradt= (t1-t2)/float((sv1-sv2))
print (gradt)

print ("########################")
print ("For P in dbar/kPa ",   deltap, deltapdkPa )
print ("########################")

# In kPa  --> Divide by 10 for dbar which is used by seawater 

t2 	= t1   
sal2 	= sal1 
p2 	= p1 + deltap  
p2kPa   = p2 * 10  

print ("t1,t2, sal1,sal2,  p1,p2 (dbar), p1,p2 (kPa): ")
print (t1,t2, sal1,sal2, p1,p2, p1kPa , p2kPa)

t = T90conv([  t1 ])
s = [ sal1  ]
p = [ p1 ]
sv1 = sw.svel(s, t, p)

t = T90conv([  t2 ])
s = [ sal2  ]
p = [ p2 ]
sv2 = sw.svel(s, t, p)

print ("sv1,sv2:")
print (sv1,sv2)
print ("Differences sv1-sv2" )
print (sv1-sv2)
print ("Ratio (p1-p2)/(sv1-sv2)  in  dbar/ (km/s) " )
gradpdbar = (p1-p2)      /float((sv1-sv2))   
print (gradpdbar)
print ("Ratio (p1-p2)/(sv1-sv2)  in  kPa / (km/s) " )
gradpkPa  = ((p1-p2)*10)/float((sv1-sv2))
print (gradpkPa )

print ("########################")
print ("For SAL" )
print ("########################")

t2 	= t1   
sal2 	= sal1 + deltas  
p2 	= p1   
p2kPa   = p2 * 10  

print ("t1,t2, sal1,sal2,  p1,p2 (dbar), p1,p2 (kPa): ")
print (t1,t2, sal1,sal2, p1,p2, p1kPa , p2kPa)

t = T90conv([  t1 ])
s = [ sal1  ]
p = [ p1 ]
sv1 = sw.svel(s, t, p)

t = T90conv([  t2 ])
s = [ sal2  ]
p = [ p2 ]
sv2 = sw.svel(s, t, p)

print ("sv1,sv2:")
print (sv1,sv2)
print ("Differences" )
print (sv1-sv2)
print ("Difference in PSU/ ( km/s) " )
grads = (sal1-sal2)/float((sv1-sv2))
print( grads)

print ("")
print ("############################################")
print (" SV Gradients for t1,p1kPa,sal1" , t1 , p1kPa ,sal1)
print ("############################################")

print ('Temp. ', gradt       , ' deg./(m/s)  ' ,  (1/float(gradt ))	, '  (m/s)/deg.  '	)    
print ('Sal.  ', grads	    , ' PSU /(m/s)  ' ,  (1/float(grads))	, '  (m/s)/PSU   '	)   
print ('Pres. ', gradpkPa    , ' kPa /(m/s)  ' ,  (1/float(gradpkPa ))	, '  (m/s)/kPa  '	)   
print ('Pres. ', gradpdbar   , ' dbar/(m/s)  ' ,  (1/float(gradpdbar))	, '  (m/s)/dbar  ')	   



### You are looking at the data you collected during the mission

print( "")
print ("############################################")
print (" What are the Baseline changes for a given parameter change( 1 km long baseline)?")
print ("############################################")

#
#  delta SV
#  -------- = gradient [ m/s / kPa ] 
#  delta p
#
# --> delta SV [km/s] = delta p   [kPa ] * gradient [ m/s / kPa ]   	(1)
# --> delta p  [kbar] = delta SV [km/s] / gradient [ m/s / kbar ]   	(2) 
#
# With SV=BSL/t --> delta BSL = delta SV * t      			(3)
# Put in (1) into (3) 
#
# Dependency of BSL on P/SAL/T Change 
# 
# delta BSL [m] = delta SV * t 
# --> 
# delta BSL [m] = delta p    [kPa]  * gradient [ m/s / kPa ]   * t [s ]      for P 
# delta BSL [m] = delta SAL  [PSU ] * gradient [ m/s / PSU ]   * t [s ]      for SALinitiy
# delta BSL [m] = delta T    [deg.] * gradient [ m/s / deg ]   * t [s ]      for Temperature 
#

# TWT is two way travel time for a given Baseline 
BSLlength	= 750		        # in meter 
# TWT 	  	= float(float(2*BSLlength) /  1500 )	# TWT in seconds 	
TWT 	  	= float(float(2*BSLlength) /  sv1  )	# TWT in seconds 	
# TWT 		= 1 

print ("TWT for Baseline with length " , BSLlength , " m ")
print ("and SV ", sv1 , "k/s ")
print ("TWT:", TWT, " s")
print ("")

# delta BSL [m] = delta p    [kPa ] * gradient [ m/s / kPa ]   * t [s ]      for P 
# delta BSL [m] = delta SAL  [PSU ] * gradient [ m/s / PSU ]   * t [s ]      for SALinitiy
# delta BSL [m] = delta T    [deg.] * gradient [ m/s / deg ]   * t [s ]      for Temperature 

dt 	= 1    # in !c
dt 	= 0.01 # in !c

ds 	= 1    # in PSU
ds 	= 0.01 # in PSU

dp   	= 1 # In kPa
dpdbar 	= 1  # In dbar

deltaBSLt       =  "%8.3f" %(dt		* (1/float(gradt    ))	* TWT * 100)			
deltaBSLSAL     =  "%8.3f" %(ds		* (1/float(grads))    	* TWT * 100)	
deltaBSLp       =  "%8.3f" %(dp		* (1/float(gradpkPa ))	* TWT * 100)	
deltaBSLpdbar   =  "%8.3f" %(dpdbar	* (1/float(gradpdbar ))	* TWT * 100)	

print ('BSLchange' ,   deltaBSLt ,   ' (cm) from Temp. change (deg.) ', dt 	 ,	" Gradient "  , gradt	    , ' deg./(m/s)  ' ,  (1/float(gradt ))	, '  (m/s)/deg.  ')
print ('BSLchange' ,   deltaBSLSAL,  ' (cm) from Sal.  change (PSU)  ', ds 	 ,	" Gradient "  , grads	    , ' PSU /(m/s)  ' ,  (1/float(grads))	, '  (m/s)/PSU   ')
print( 'BSLchange' ,   deltaBSLp ,   ' (cm) from Pres. change (kPa)  ', dp        ,	" Gradient "  , gradpkPa    , ' kPa /(m/s)  ' ,  (1/float(gradpkPa ))	, '  (m/s)/kPa  ')
print( 'BSLchange' ,   deltaBSLpdbar,' (cm) from Pres. change (dbar) ', dpdbar    ,	" Gradient "  , gradpdbar   , ' dbar/(m/s)  ' ,  (1/float(gradpdbar))	, '  (m/s)/dbar  ')


