import math

#n1= float(raw_input("Define the refractive index n1(must be 1 or greater):")
#n2= float(raw_input("Define the refractive index n2(must be 1 or greater):") 

x = float(raw_input("Define the angle of incidence:"))

x = math.radians(x) 

z = math.sin(x)*(n1/n2)

t = asin(z) 
 

#t is the transmission angle in radians 
#x is the incident angle in radians

#now, to find the reflection and transmission coefficients. 

#reflection coefficient

R = math.tan(x-t)/math.tan(x-t) 

#transmission coefficient

T = 2*math.sin(t)*math.cos(x)/math.sin(x+t)*math.cos(x-t) 

#check to make sure conservation of energy is preserved we can check if this equals 1: 

(R**2+T**2)(n2*math.cos(t)/n1*math.cos(x)) 


