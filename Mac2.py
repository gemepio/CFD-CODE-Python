# 2D Marker and Cell Method for Solving Incompressible Flow
from pylab import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math      
import numpy as np



def main():
	#domain_properties()
	lx=1.0
	ly=1.0
	gx=0.0
	gy=-100.0
	rho1=1.0
	rho2=2.0
	m0=0.01
	rro=rho1
	unorth=0.0
	usouth=0.0
	veast=0.0
	vwest=0.0
	time=0.0
	rad=0.15
	xc=0.5
	yc=0.7

	#grid_gen()
	nx=48
	ny=48
	dt =0.00125
	nstep=100
	maxit=300
	maxErr=0.001
	beta =1.2
	
	#=====================initialize zero arrarys================
	u=np.zeros([nx+1,ny+2])
	v=np.zeros([nx+2,ny+1])
	p=np.zeros([nx+2,ny+2])
	ut=np.zeros([nx+1,ny+2])
	vt=np.zeros([nx+2,ny+1])
	tmp1=np.zeros([nx+2,ny+2])
	tmp2=np.zeros([nx+2,ny+2])
	uu=np.zeros([nx+1,ny+1])
	vv=np.zeros([nx+1,ny+1]);
	x=np.zeros(nx+2)
	y=np.zeros(ny+2)
	#=====================Construct the grid=====================
	dx=lx/nx
	dy=ly/ny
	for i in range(0, nx+2):
		x[i]=dx*(i-0.5)
	for j in range(0, ny+2):
		y[j]=dy*(j-0.5)
	#Set Density
	r=np.zeros([nx+2, ny+2])+rho1
	for i in range (1,nx+1):
		for j in range (1,ny+1):
			if ( (x[i]-xc)**2+(y[j]-yc)**2 < rad**2 ):
				r[i,j]=rho2

	#f =open('data1','w')
	
	#=====================Time Loop=====================
	for step in range (1, 5):
		#=================Tangential Velocity=================
		u[0:nx+1,0]=2*usouth-u[0:nx+1,1]
		u[0:nx+1, ny+1]=2*unorth-u[0:nx+1,ny]
		v[0,0:ny+1]=2*vwest-v[1,0:ny+1]
		v[nx+1,0:ny+1]=2*veast-v[nx,0:ny+1]
		
		#=====================Temporary u=====================
		for i in range (1,nx):
			for j in range (1, ny+1):
				ut[i,j]=u[i,j]+dt*(-0.25*(((u[i+1,j]+u[i,j])**2 \
						-(u[i,j]+u[i-1,j])**2)/dx+((u[i,j+1]+u[i,j])* \
						(v[i+1,j]+v[i,j])-(u[i,j]+u[i,j-1])* \
						(v[i+1,j-1]+v[i,j-1]))/dy)+ \
						m0/(0.5*(r[i+1,j]+r[i,j]))*((u[i+1,j]-2*u[i,j] \
						+u[i-1,j])/dx**2+(u[i,j+1]-2*u[i,j]+u[i,j-1])/dy**2)+gx)
		
		#=====================Temporary v=====================
		for i in range (1,nx+1):
			for j in range (1,ny):
				vt[i,j]=v[i,j]+dt*(-0.25*(((u[i,j+1]+u[i,j])*(v[i+1,j] \
					+v[i,j])-(u[i-1,j+1]+u[i-1,j])*(v[i,j]+v[i-1,j]))/dx+ \
					((v[i,j+1]+v[i,j])**2-(v[i,j]+v[i,j-1])**2)/dy)+ \
					m0/(0.5*(r[i,j+1]+r[i,j]))*((v[i+1,j]-2*v[i,j]+v[i-1,j])/dx**2 \
					+(v[i,j+1]-2*v[i,j]+v[i,j-1])/dy**2)+gy)
		
		#=====================P(i,j) source term=====================
		
		rt=r
		lrg=1000
		rt[0:nx+1,0]=lrg
		rt[0:nx+1,ny+1]=lrg
		rt[0,0:ny+1]=lrg
		rt[nx+1,0:ny+1]=lrg
		
		for i in range (1,nx+1):
			for j in range (1,ny+1):
				tmp1[i,j]=(0.5/dt)*( (ut[i,j]-ut[i-1,j])/dx \
							+(vt[i,j]-vt[i,j-1])/dy)
				tmp2[i,j]=1.0/( (1/dx)*( 1/(dx*(rt[i+1,j]+rt[i,j]))+ \
					 	  1/(dx*(rt[i-1,j]+rt[i,j])))+ \
						  (1/dy)*(1/(dy*(rt[i,j+1]+rt[i,j]))+ \
						  1/(dy*(rt[i,j-1]+rt[i,j]))))
		
		#=====================Solve for pressure=====================
		for it in range (1,maxit+1):
			oldp=p
			for i in range (1,nx+1):
				for j in range (1,ny+1):
					p[i,j]=(1.0-beta)*p[i,j]+beta*tmp2[i,j]* \
							((1/dx)*(p[i+1,j]/(dx*(rt[i+1,j]+rt[i,j])) \
							+p[i-1,j]/(dx*(rt[i-1,j]+rt[i,j]))) \
							+(1/dy)*(p[i,j+1]/(dy*(rt[i,j+1]+rt[i,j])) \
							+p[i,j-1]/(dy*(rt[i,j-1]+rt[i,j]))) -tmp1[i,j])
			if np.amax(abs(oldp-p)) <maxErr:
				break
		
		#=====================Correc u-velocity=====================
		for i in range (1,nx):
			for j in range (1,ny+1):
				u[i,j]=ut[i,j]-dt*(2.0/dx)*(p[i+1,j]-p[i,j])/(r[i+1,j]+r[i,j])
		#=====================Correc v-velocity=====================
		for i in range (1,nx+1):
			for j in range (1,ny):
				v[i,j]=vt[i,j]-dt*(2.0/dy)*(p[i,j+1]-p[i,j])/(r[i,j+1]+r[i,j])
		
		#=====================Advect Density with Center difference=====================
		rold=r
		for i in range (1,nx+1):
			for j in range (1,ny+1):
				r[i,j]=rold[i,j]-(0.5*dt/dx)*(u[i,j]*(rold[i+1,j]+rold[i,j]) \
						-u[i-1,j]*(rold[i-1,j]+rold[i,j]))-(0.5*dt/dy)*(v[i,j]*(rold[i,j+1] \
						+rold[i,j])-v[i,j-1]*(rold[i,j-1]+rold[i,j])) \
						+(m0*dt/dx/dx)*(rold[i+1,j]-2.0*rold[i,j]+rold[i-1,j]) \
						+(m0*dt/dy/dy)*(rold[i,j+1]-2.0*rold[i,j]+rold[i,j-1])
		
		#=====================Results Visualization=====================
		time=time+dt
		print time
		uu[0:nx,0:ny]=0.5*(u[0:nx,1:ny+1]+u[0:nx,0:ny])
		vv[0:nx,0:ny]=0.5*(v[1:nx+1,0:ny]+v[0:nx,0:ny])
		#f.write(np.flipud(np.rot90(r)))
		np.savetxt('u'+str(step)+'.out', np.flipud(np.rot90(u)), fmt='%.4e',delimiter=' ')
		np.savetxt('v'+str(step)+'.out', np.flipud(np.rot90(v)), fmt='%.4e',delimiter=' ')
		np.savetxt('p'+str(step)+'.out', np.flipud(np.rot90(p)), fmt='%.4e',delimiter=' ')
		np.savetxt('rho'+str(step)+'.out', np.flipud(np.rot90(r)), fmt='%.4e',delimiter=' ')
		#for i in range (0,nx+1):
			#xh[i]=dx*(i-1)
		#for j in range (0,ny+1):
			#yh[i]=dy*(j-1)
		#Q=quiver(np.flipud(np.rot90(uu)),np.flipud(np.rot90(vv)))
		#plt.figure()
		#hold(False)
		#plt.contour(x,y,np.flipud(np.rot90(r)))
		#hold(True)
		#qk= quiverkey(Q,0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W',
        #       fontproperties={'weight': 'bold'})
		#show()
	#f.close()		
if __name__ == '__main__':
	main()
#------------------------------------------------------#
