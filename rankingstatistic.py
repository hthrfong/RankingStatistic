from sympy import *
from scipy import special
import math
import numpy as np
import matplotlib.pyplot as plt
import random

def log_incomplete_Gamma(a,x):
    # scipy.special.gammaincc is defined as
    # 1 / gamma(a) * integral(exp(-t) * t**(a-1), t=x..inf)
    # does gamma(a) * scipy.special.gammaincc to get the regular definition.

    # Gamma(a,x) = Gamma(a)-gamma(a,x)
    # ln(Gamma(a,x)) = ln(Gamma(a)-gamma(a,x))
    # ln(Gamma(a,x)) = ln(Gamma(a)*(1-gamma(a,x)/Gamma(a)))
    # ln(Gamma(a,x)) = ln(Gamma(a)+ln(1-gamma(a,x)/Gamma(a)) where gamma(a,x)/Gamma(a) = regularized Gamma function P(a,x)
    # ln(Gamma(a,x)) \approx ln(Gamma(a)) when a >> 1
    #
    #   Example: >>> special.gammainc(100,50)
    #            3.2000653245851495e-10
    #   Example: >>> special.gammainc(370,200)
    #            4.2466647446310334e-27
    
    log_incomplete_G = special.gammaln(a)
    return log_incomplete_G

def list_sum(n,N,rho,tjtk, y):
    # finds summation of te first term in equation 23
    f = []
    for i in range(len(n)):
        f.append( y(n[i]) )
        #if f[-1] < 1e-16*max(f):
        #    break
    f = np.array(f)
    #cumsum
    Sum = 0
    SummedArray = []
    for i in range(len(f)):
        SummedArray.append(Sum)
        Sum += f[i]
    return Sum, SummedArray,f 

def derivative(x,y,i):
    diff1 = (y(x[i+1])-y(x[i]))/(x[i+1]-x[i])
    return diff1

def peak_finder(i, i_sample, n, x, y, option):
    if int(i_sample/n**2) <= 1:
        if option == "-":
            i -= 1
        if option == "+":
            i += 1
    if int(i_sample/n**2) > 1 and option == "-":
        if (i - int(i_sample/n**2)) < 0:
            i -= 1
        if (i - int(i_sample/n**2)) > 0:
            i -= int(i_sample/n**2)
    if int(i_sample/n**2) > 1 and option == "+":
        if (i + int(i_sample/n**2)) > len(x)-1:
            i += 1
        if (i + int(i_sample/n**2)) < len(x)-1:
            i += int(i_sample/n**2)
    diff = derivative(x,y,i)
    print "\t \t diff and i", diff, i
    return i, diff, n

def sampling(i,i_sample,s,x,y,diff,number):
    if diff > number:
        while diff > number and i > 0:
            i,diff,s = peak_finder(i,i_sample,s,x,y,"+")
    if diff < number:
        while diff < number and i > 0:
            i,diff,s = peak_finder(i,i_sample,s,x,y,"-")
    s += 1
    return i, diff, s

def iterative_finding(N,x,y,i,i_sample,peak_diff,s):
    diff = derivative(x,y,i)
    print "\t i", i, diff, peak_diff
    i, diff, s = sampling(i,i_sample,s,x,y,diff,peak_diff)
    return i, diff, s

N = 100. # number of dimensions
rho = 10. # nominal SNR
n = np.arange(0,N) # range of n in summation
t = 0.9 # value of (t_j \cdot t_k)

y_num = (lambda x: ( x*np.log( np.sqrt(2)/(rho*t) )
                         + special.gammaln(N)
                         + special.gammaln((x+1)/2)
                         - special.gammaln(N-x)
                         - special.gammaln(x+1) ) )

#y_den = (lambda x: ( x*np.log( np.sqrt(2)/(rho) )
#                         + special.gammaln(N)
#                         + special.gammaln((x+1)/2)
#                         - special.gammaln(N-x)
#                         - special.gammaln(x+1) ) )

#y = (lambda x: ( x*np.log( np.sqrt(2)/(rho*t) )
#               + special.gammaln((x+1)/2)
#               - special.gammaln(N-x)
#               - special.gammaln(x+1) ) )

i_sample = int(N/4)
i = random.randint(0,N-1)
'''s = 1
while int(i_sample/s**2) >= 1 and i > 0:
    i,diff, s = iterative_finding(N,n,y,i,i_sample,0,s)
    if i < 0:
        break
peak_index = i
peak_diff = diff
print "Peak", peak_index,peak_diff,s
s = 1
while int(i_sample/s**2) >= 1 and i > 0:
    i,diff, s = iterative_finding(N,n,y,peak_index-i_sample,i_sample,peak_diff*0.8,s)
    if i < 0:
        break
lowerbound_index = i
lowerbound_diff = diff
print "Lower bound", lowerbound_index,lowerbound_diff,s'''
diff = derivative(n,y_num,i)
print "random i", i, diff
if diff < 0:
    while diff < 0:
        i, diff, s = peak_finder(i,i_sample,1,n,y_num,"-")
        print i, diff
    if diff > 0:
        while diff > 0:
            i, diff, s = peak_finder(i,i_sample,2,n,y_num,"+")
            print i, diff
        if diff < 0:
            while diff < 0:
                i, diff, s = peak_finder(i,i_sample,4,n,y_num,"-")
                print i, diff
            if diff > 0:
                while diff > 0:
                    i, diff, s = peak_finder(i,i_sample,8,n,y_num,"+")
                    print i, diff
elif diff > 0:
    while diff > 0:
        i, diff, s = peak_finder(i,i_sample,1,n,y_num,"+")
        print i, diff
    if diff < 0:
        while diff < 0:
            i, diff, s = peak_finder(i,i_sample,2,n,y_num,"-")
            print i, diff
        if diff > 0:
            while diff > 0:
                i, diff, s = peak_finder(i,i_sample,4,n,y_num,"+")
                print i, diff
            if diff < 0:
                while diff < 0:
                    i, diff, s = peak_finder(i,i_sample,8,n,y_num,"-")
                    print i, diff

            
peak_index_num = i
peak_diff_num = diff
'''
diff = derivative(n,y_den,i)
print "random i", i, diff
if diff < 0:
    while diff < 0:
        i, diff, s = peak_finder(i,i_sample,1,n,y_den,"-")
        print i, diff
    if diff > 0:
        while diff > 0:
            i, diff, s = peak_finder(i,i_sample,2,n,y_den,"+")
            print i, diff
        if diff < 0:
            while diff < 0:
                i, diff, s = peak_finder(i,i_sample,4,n,y_den,"-")
                print i, diff
            if diff > 0:
                while diff > 0:
                    i, diff, s = peak_finder(i,i_sample,8,n,y_den,"+")
                    print i, diff
elif diff > 0:
    while diff > 0:
        i, diff, s = peak_finder(i,i_sample,1,n,y_den,"+")
        print i, diff
    if diff < 0:
        while diff < 0:
            i, diff, s = peak_finder(i,i_sample,2,n,y_den,"-")
            print i, diff
        if diff > 0:
            while diff > 0:
                i, diff, s = peak_finder(i,i_sample,4,n,y_den,"+")
                print i, diff
            if diff < 0:
                while diff < 0:
                    i, diff, s = peak_finder(i,i_sample,8,n,y_den,"-")
                    print i, diff
'''
#peak_index_den = i
#peak_diff_den = diff

print "Peak:", peak_index_num, peak_diff_num#, peak_index_den, peak_diff_den

y_num2 = (lambda x: ( np.exp( x*np.log( np.sqrt(2)/(rho*t) )
                         + special.gammaln(N)
                         + special.gammaln((x+1)/2)
                         - special.gammaln(N-x)
                         - special.gammaln(x+1) - y_num(n[peak_index_num]) ) ) )

y_den2 = (lambda x: ( np.exp( x*np.log( np.sqrt(2)/(rho) )
                        + special.gammaln(N)
                         + special.gammaln((x+1)/2)
                         - special.gammaln(N-x)
                         - special.gammaln(x+1) - y_den(n[peak_index_den]) ) ) )

diff = derivative(n,y_num2,i)
i = peak_index_num - i_sample
number = abs(derivative(n,y_num2,peak_index_num)*0.1)
print "number:",number
if diff < number:
    while diff < number:
        i, diff, s = peak_finder(i,i_sample,2,n,y_num2,"+")
        print i, diff
    if diff > number:
        while diff > number:
            i, diff, s = peak_finder(i,i_sample,4,n,y_num2,"-")
            print i, diff
elif diff > number:
    while diff > number:
        i, diff, s = peak_finder(i,i_sample,2,n,y_num2,"-")
        print i, diff
    if diff < number:
        while diff < number:
            i, diff, s = peak_finder(i,i_sample,4,n,y_num2,"+")
            print i, diff

lowerbound_index_num = i
lowerbound_diff_num = diff

theSum0, b0, f0 = list_sum(n[lowerbound_index_num:lowerbound_index_num+peak_index_num],N,rho,t,y_num2)
theSum, b,f = list_sum(n,N,rho,t,y_num2)
theSum2, b2,f2 = list_sum(n,N,rho,1.,y_num2)

print np.exp(rho**2*(1-t**2))*(t**(N-1))*(theSum0/theSum2)*(y_num2(n[peak_index_num])/

'''
print "finding lower bound"
p = 0.001

diff2 = derivative(n,y,peak_index-i_sample)
i = peak_index - i_sample
if diff2 < abs(peak_diff*p):
    while diff2 < abs(peak_diff*p):
        i, diff2, s = peak_finder(i,i_sample,2,n,y,"+")
        print i, diff2
    if diff2 > abs(peak_diff*p):
        while diff2 > abs(peak_diff*p):
            i, diff2, s = peak_finder(i,i_sample,4,n,y,"-")
            print i, diff2
elif diff2 > abs(peak_diff*p):
    while diff2 > abs(peak_diff*p):
        i, diff2, s = peak_finder(i,i_sample,2,n,y,"-")
        print i, diff2
    if diff2 < abs(peak_diff*p):
        while diff2 < abs(peak_diff*p):
            i, diff2, s = peak_finder(i,i_sample,4,n,y,"+")
            print i, diff2

lowerbound_index = i
lowerbound_diff = diff2

print "Peak:", peak_index,peak_diff
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
#plt.show()

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
