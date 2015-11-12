from sympy import *
from scipy import special
import math
import numpy as np
import matplotlib.pyplot as plt
import random

def incomplete_Gamma(a,x):
    # scipy.special.gammaincc is defined as
    # 1 / gamma(a) * integral(exp(-t) * t**(a-1), t=x..inf)
    # does gamma(a) * scipy.special.gammaincc to get the regular definition.
    gamma = special.gamma(a)
    incomplete_G = special.gammaincc(a,x) * gamma
    return incomplete_G

def list_sum(n, N, rho,tjtk):
    # finds summation of te first term in equation 23
    f = []
    for i in range(len(n)):
        #f.append( (rho*tjtk)**(-n[i])/math.factorial(n[i]) * 2.**(n[i]/2) * incomplete_Gamma((n[i]+1)/2,(rho*tjtk)**2/2) * (-(N-n[i]-2)/(math.factorial(N-n[i]-1))))
        f.append( (np.sqrt(2)/(rho*tjtk))**(n[i]) * incomplete_Gamma((n[i]+1)/2,(rho*tjtk)**2/2) * (1./special.gamma(N-n[i])) * (1./special.gamma(n[i]+1)) )
        #f.append( (np.sqrt(2)/rho*tjtk)**(n[i]) * incomplete_Gamma((n[i]+1)/2,(rho*tjtk)**2/2) * (1./special.gamma(N-n[i])) * (1./special.gamma(n[i]+1)) )
    f = np.array(f)
    Sum = 0
    SummedArray = []
    for i in range(len(f)):
        SummedArray.append(Sum)
        Sum += f[i]
    return Sum, SummedArray,f 

#def summation(n, N, rho, tjtk):
#    y = (np.sqrt(2)/(rho*tjtk))**(n[i]) * incomplete_Gamma((n[i]+1)/2,(rho*tjtk)**2/2) * (1./special.gamma(N-n[i])) * (1./special.gamma(n[i]+1))
#    return y

def derivative(x,y,i):
    diff1 = (y(x[i+1])-y(x[i]))/(x[i+1]-x[i])
    return diff1

def peak_finder(i, i_sample, n, x, y, option):
    if int(i_sample/n**2) < 1:
        if option == "-":
            i -= 1
        elif option == "+":
            i += 1
    if option == "-":
        if (i - int(i_sample/n**2)) < 0:
            i -= 1
        else:
            i -= int(i_sample/n**2)
    if option == "+":
        if (i + int(i_sample/n**2)) > len(x)-1:
            i += 1
        else:
            i += int(i_sample/n**2)
    diff = derivative(x,y,i)
    return i, diff, n

N = 171. # number of dimensions
rho = 10. # nominal SNR
n = np.arange(0,N) # range of n in summation
t = 0.9 # value of (t_j \cdot t_k)

y = lambda x: (np.sqrt(2)/(rho*t))**(x) * incomplete_Gamma((x+1)/2,(rho*t)**2/2) * (1./special.gamma(N-x)) * (1./special.gamma(x+1))

i_sample = int(N/4)
i = random.randint(0,len(n)-1)
diff = derivative(n,y,i)
print "random i", i, diff
if diff < 0:
    while diff < 0:
        i, diff, s = peak_finder(i,i_sample,1,n,y,"-")
        print i, diff
    if diff > 0:
        while diff > 0:
            i, diff, s = peak_finder(i,i_sample,2,n,y,"+")
            print i, diff
        if diff < 0:
            while diff < 0:
                i, diff, s = peak_finder(i,i_sample,4,n,y,"-")
                print i, diff
            if diff > 0:
                while diff > 0:
                    i, diff, s = peak_finder(i,i_sample,8,n,y,"+")
                    print i, diff
elif diff > 0:
    while diff > 0:
        i, diff, s = peak_finder(i,i_sample,1,n,y,"+")
        print i, diff
    if diff < 0:
        while diff < 0:
            i, diff, s = peak_finder(i,i_sample,2,n,y,"-")
            print i, diff
        if diff > 0:
            while diff > 0:
                i, diff, s = peak_finder(i,i_sample,4,n,y,"+")
                print i, diff
            if diff < 0:
                while diff < 0:
                    i, diff, s = peak_finder(i,i_sample,8,n,y,"-")
                    print i, diff

            
peak_index = i
peak_diff = diff

print "Peak:", peak_index, peak_diff
print "finding lower bound"

diff2 = derivative(n,y,peak_index-i_sample)
i = peak_index - i_sample
if diff2 < abs(peak_diff*0.8):
    while diff2 < abs(peak_diff*0.8):
        i, diff2, s = peak_finder(i,i_sample,2,n,y,"+")
        print i, diff2
    if diff2 > abs(peak_diff*0.8):
        while diff2 > abs(peak_diff*0.8):
            i, diff2, s = peak_finder(i,i_sample,4,n,y,"-")
            print i, diff2
elif diff2 > abs(peak_diff*0.8):
    while diff2 > abs(peak_diff*0.8):
        i, diff2, s = peak_finder(i,i_sample,2,n,y,"-")
        print i, diff2
    if diff2 < abs(peak_diff*0.8):
        while diff2 < abs(peak_diff*0.8):
            i, diff2, s = peak_finder(i,i_sample,4,n,y,"+")
            print i, diff2

lowerbound_index = i
lowerbound_diff = diff2

print "Lower bound:", lowerbound_index, lowerbound_diff



theSum0, b0, f0 = list_sum(n[lowerbound_index:lowerbound_index+peak_index],N,rho,tjtk=t)

theSum, b,f = list_sum(n,N,rho,tjtk=t)
theSum2, b2,f2 = list_sum(n,N,rho,tjtk=1.)

print "iterative vs total sum:", theSum0,theSum

plt.plot(n,y(n))
plt.plot(n[lowerbound_index],y(n[lowerbound_index]),marker="o")
plt.plot(n[peak_index],y(n[peak_index]),marker="o")
plt.figure()
plt.plot(n,b)
plt.plot(n[lowerbound_index:lowerbound_index+peak_index],b0,"--",color="red")
plt.show()

'''
# finds second term in equation 23
#c = 1./(math.factorial(N-1)) * (rho*t)**(-N+1) * 2.**((N-1)/2) * incomplete_Gamma(N/2,(rho*t)**2/2)
#c2 = 1./(math.factorial(N-1)) * (rho)**(-N+1) * 2.**((N-1)/2) * incomplete_Gamma(N/2,(rho)**2/2)

print "for tjtk=%.2f:" %t
print "\t sum: %.6e" %(theSum)
print "for tjtk=1:"
print "\t sum: %.6e" %(theSum2)

#print "total value:", t**(N-1)*np.exp(0.5*(rho*t)**2) * (theSum/theSum2)

# plots each summed term as a function of n
plt.plot(n,b,'-o',label="tjtk<1")
#plt.plot(n,b2,'-o',label="tjtk=1")
plt.xlabel("n = np.arange(0,N)")
plt.ylabel("sum at n")
plt.legend(loc="upper left")
#plt.figure()
#plt.plot(d)
plt.show()
#plt.savefig("summation.pdf")
'''
