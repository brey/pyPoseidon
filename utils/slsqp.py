import numpy as np
from scipy.optimize import minimize, fmin_slsqp
from pylab import plot,show,figure,title

def slsqp(R,V,maxR,rho,b,rmax,dp,k,vmax):
        dpmin=10.e2  # minimum value of  pressure drop P_central - P_env(=101kPa).
        dpmax=200.e2   # maximum value of  pressure drop P_central - P_env(=101kPa).
        kmin=0  # low limit of parameter k (=xn-.5-> k=0-> x=0.5)
        kmax=0.15 # upper limit for k (->xn=.65)  WHY?
        bmin=0.8
        bmax=1.8

        def func(p,x):
#                return (p[0]/rho*(p[1]/x)**p[0]*p[2]*np.exp(-(p[1]/x)**p[0]))**(0.5+(x-p[1])/(maxR-p[1])*p[3])
             #   return ((p[1]/x)**p[0]*vmax**2*np.exp(1.-(p[1]/x)**p[0]))**(0.5+(x-p[1])/(maxR-p[1])*p[2])
             #   return (p[0]/rho*(p[1]/x)**p[0]*dp*np.exp(-(p[1]/x)**p[0]))**(0.5+(x-p[1])/(maxR-p[1])*p[2])
                 return (p[0]/rho*(rmax/x)**p[0]*dp*np.exp(-(rmax/x)**p[0]))**(0.5+(x-rmax)/(maxR-rmax)*p[1])
        def errf(p,x,y):
            return np.sum((func(p,x)-y)**2)

#       def cf(p,x,y):
#	    return vmax-func(p,p[1])
        def cf(p,x,y):
                 return p[0]-vmax**2*rho*np.exp(1.)/p[1]

#	p0=[1.2,20000,400,.5]
	p0=[1.2,20000,.5]

#       bp=[(.8,1.8),(5000,R.min()*.99),(dpmin,dpmax),(kmin,kmax)]
        bp=[(.8,1.8),(5000,R.min()*.99),(kmin,kmax)]

        res1 = minimize(errf, p0, args=(R, V), method='L-BFGS-B', bounds=bp)
        print 'L-BFGS-B', res1.x
        res2 = minimize(errf, p0, args=(R, V), method='SLSQP', bounds=bp, tol=1e-1)
        print 'SLSQP', res2.x, res2.fun, res2.message
        res = fmin_slsqp(errf, p0, bounds=bp, args=(R, V),f_ieqcons=cf,acc=1e-1)
        print res

        

        print 'Rmax= ', res1.x[1],res2.x[1]
        print 'F(Rmax)= ', cf(res1.x,0.,0.), cf(res2.x,0.,0.)
        return func,res1.x,res2.x, res


        


