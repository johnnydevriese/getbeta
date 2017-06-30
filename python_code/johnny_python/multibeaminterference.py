import math as m 


#n1= float(raw_input("Define the refractive index n1(must be 1 or greater):")
#n2= float(raw_input("Define the refractive index n2(must be 1 or greater):") 

x = float(raw_input("Define the angle of incidence:"))

x = m.radians(x) 

z = m.sin(x)*(n1/n2)

t = asin(z) 
 

#t is the transmission angle in radians 
#x is the incident angle in radians

#now, to find the reflection and transmission coefficients. 

#reflection coefficient

R = m.tan(x-t)/m.tan(x+t) 

#transmission coefficient

T = 2*m.sin(t)*m.cos(x)/m.sin(x+t)*m.cos(x-t) 

#moving onto multiple beam interference (p.87 Fowles)
#finding the phase difference in a medium...

#defining variables
l = float(raw_input("What is the wavelength in the medium between two reflecting surfaces?")) 
d = float(raw_input("What is the distance between the two mirrors?")) 
Ei = float(raw_input("What is the initial Energy(Amplitude)")) 

delta =  (4*m.pi/l)*n2*d*m.cos(t) 

#taking 'delta' into account as the phase difference we can compute E_t (Energy Transmitted) with a geometric series

E_t = Ei*T**2/1-R**2*m.e**j*delta 

#we can also compute Intesity of the transmitted waves

I_t = abs(E_t)**2








