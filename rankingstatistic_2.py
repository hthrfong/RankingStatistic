from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import random
import time

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

start_time = time.time()

rho = 50.
t = [0.9, 1.0] # (tj \cdot tk) for [numerator, denominator]
#N_array = np.arange(10000,5010000,50000)
N_array = np.arange(1000,1010000,50000)
total = []

for N in N_array:
    print "N=%i" %N
    n = np.arange(1,N)
    log_y = [(lambda x: ( x*np.log( np.sqrt(2)/(rho*t[0]) )
                     + special.gammaln(N)
                     + special.gammaln((x+1)/2)
                     - special.gammaln(N-x)
                     - special.gammaln(x+1) ) ),
             (lambda x: ( x*np.log( np.sqrt(2)/(rho*t[1]) )
                     + special.gammaln(N)
                     + special.gammaln((x+1)/2)
                     - special.gammaln(N-x)
                     - special.gammaln(x+1) ) )]

    peak_index = []
    peak_deriv = []
    bound_index = []
    bound_deriv = []

    for m in range(len(t)):
        peak_index.append(log_y[m](n).argmax())

    y = [(lambda x: ( np.exp( x*np.log( np.sqrt(2)/(rho*t[0]) )
                     + special.gammaln(N)
                     + special.gammaln((x+1)/2)
                     - special.gammaln(N-x)
                     - special.gammaln(x+1) - log_y[0](n[peak_index[0]])) ) ),
         (lambda x: ( np.exp( x*np.log( np.sqrt(2)/(rho*t[1]) )
                     + special.gammaln(N)
                     + special.gammaln((x+1)/2)
                     - special.gammaln(N-x)
                     - special.gammaln(x+1) - log_y[1](n[peak_index[1]]))) )]

    for m in range(len(t)):
        n2 = n[:peak_index[m]]
        i = np.where(y[m](n2)==find_nearest(y[m](n2),0.01))[0]
        bound_index.append(i)

    cumsum = [ [],[] ]
    totalsum = [ [],[] ]

    for m in range(len(t)):
        cumsum[m] = np.cumsum(y[m](n[bound_index[m]:bound_index[m]+peak_index[m]]))
        totalsum[m] = np.cumsum(y[m](n))
        print "cumsum vs totalsum:", cumsum[m][-1],totalsum[m][-1]

    print "\t End time:", time.time()-start_time

    logterm = (N-1)*np.log(t[0]) - (0.5*rho**2*(1-t[0]**2)) + log_y[0](n[peak_index[0]])/log_y[1](n[peak_index[1]]) + np.log((cumsum[0][-1]/cumsum[1][-1]))

    total.append(logterm)

line=[]
for N in range(len(N_array)):
    line.append(total[N]-(N_array[N]-1)*np.log(t[0]))
plt.plot(N_array,line,marker="o",linestyle="-")
plt.xlabel('N')
#plt.ylabel('log(equation 24)')
plt.show()

#plt.plot(n,y[0](n))
#plt.plot(n[peak_index[0]],y[0](n[peak_index[0]]),marker="o")
#plt.plot(n[bound_index[0]],y[0](n[bound_index[0]]),marker="o")
#plt.figure()
#plt.plot(n,y_log[0](n))
#plt.plot(n[peak_index[0]],y_log[0](n[peak_index[0]]),marker="o")
