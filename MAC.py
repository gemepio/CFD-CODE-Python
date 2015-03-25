# 2D Marker and Cell Method for Solving Incompressible Flow
from pylab import *
#import matplotlib.pyplot as plt
#from scipy.optimize import fsolve
import math        

def GridGeneration(xlength, yhight):
	global node_x, node_y
	#print("This function generate the MAC grid")
	#Horizontal Length of domain in meters
	#length = 0.2
	#Vertical hight of domain in meters
	#hight = 0.1
	#grid size in x direction
	dx = 0.01
	#grid size in y direction
	dy = 0.002
	node_x = arange (0, xlength, dx)
	node_y = arange (0, yhight,  dy)
	
def main():
	x = 0.2
	y = 0.1
	GridGeneration(x, y)
	print node_x
	funcz=node_x**2+1
	#plt.plot(funcz)
	print("The Mac grid has been generated")

if __name__ == '__main__':
	main()
