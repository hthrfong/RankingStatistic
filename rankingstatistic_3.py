from scipy.special import *
import numpy as np
import matplotlib.pyplot as plt

def incomplete_Gamma(a,x):
    # scipy.special.gammaincc is defined as
    # 1 / gamma(a) * integral(exp(-t) * t**(a-1), t=x..inf)
    # this function does gamma(a) * scipy.special.gammaincc(a,x).

    # Gamma(a,x) = Gamma(a)-gamma(a,x)
    # ln(Gamma(a,x)) = ln(Gamma(a)-gamma(a,x))
    # ln(Gamma(a,x)) = ln(Gamma(a)*(1-gamma(a,x)/Gamma(a)))
    # ln(Gamma(a,x)) = ln(Gamma(a)+ln(1-gamma(a,x)/Gamma(a)) where gamma(a,x)/Gamma(a) = regularized Gamma function P(a,x)
    # ln(Gamma(a,x)) \approx ln(Gamma(a)) when a >> 1
    
    incomplete_G = gamma(a) * gammaincc(a,x)
    return incomplete_G

def log_y(rho,t,N,approx,alt):
    # Take the logarithm of the summation term
    if approx==False:
        if alt==True:
            return lambda x: (gammaln(N)-(gammaln(x+1)+gammaln(N-x))) + x*np.log(rho*t/np.sqrt(2)) + np.log(incomplete_Gamma((N-x)/2,0.5*(rho*t)**2))
        if alt==False:
            return lambda x: (gammaln(N)-(gammaln(x+1)+gammaln(N-x))) - x*np.log(rho*t/np.sqrt(2)) + np.log(incomplete_Gamma((x+1)/2,0.5*(rho*t)**2))
    if approx==True:
        if alt==True:
            return lambda x: (gammaln(N)-(gammaln(x+1)+gammaln(N-x))) + x*np.log(rho*t/np.sqrt(2)) + gammaln((N-x)/2)
        if alt==False:
            return lambda x: (gammaln(N)-(gammaln(x+1)+gammaln(N-x))) - x*np.log(rho*t/np.sqrt(2)) + gammaln((x+1)/2)

def exp_log_y(rho,t,N,log_peak,approx,alt):
    # Take the exponential of the NORMALISED logarithm of the summation term
    if approx==False:
        if alt==True:
            return lambda x: np.exp(gammaln(N)-(gammaln(x+1)+gammaln(N-x)) + x*np.log(rho*t/np.sqrt(2)) + np.log(incomplete_Gamma((N-x)/2,0.5*(rho*t)**2)) - log_peak)
        if alt==False:
            return lambda x: np.exp(gammaln(N)-(gammaln(x+1)+gammaln(N-x)) - x*np.log(rho*t/np.sqrt(2)) + np.log(incomplete_Gamma((x+1)/2,0.5*(rho*t)**2)) - log_peak)
    if approx==True:
        if alt==True:
            return lambda x: np.exp(gammaln(N)-(gammaln(x+1)+gammaln(N-x)) + x*np.log(rho*t/np.sqrt(2)) + gammaln((N-x)/2) - log_peak)
        if alt==False:
            return lambda x: np.exp(gammaln(N)-(gammaln(x+1)+gammaln(N-x)) - x*np.log(rho*t/np.sqrt(2)) + gammaln((x+1)/2) - log_peak)

NN = np.linspace(1000,10000,50) # number of dimensions
rho = 10. # nominal SNR
t = 0.9 # (t_j \cdot t_k)
alt = True # alt=True --> rho^n*x^(N-n-1), alt=False --> rho^(N-n-1)*x^n (binomial expansion option)

keep = [[],[]] # probability

for N in NN:
    for alt in [True,False]:
        n = np.arange(0,N-1) # range of n in summation

        if alt==True:
            idx_100 = np.where((N-n) <= 100) # indices where (N-n) <= 100 (i.e. range where Gamma(a,x) != Gamma(a) ).
            idx = np.where((N-n) > 100) # indices where (N-n) > 100 (i.e. range where Gamma(a,x) = Gamma(a) ).
        if alt==False:
            idx_100 = np.where((n+1) <= 100)
            idx = np.where((n+1) > 100)

        peak_num = max(log_y(rho,t,N,True,alt)(n[idx]))
        peak_den = max(log_y(rho,1,N,True,alt)(n[idx]))
    
        exp_log_y_num = [exp_log_y(rho,t,N,peak_num,False,alt),exp_log_y(rho,t,N,peak_num,True,alt)] # [exp(log(x0/xm)),...,exp(log(xn/xm))]
        exp_log_y_den = [exp_log_y(rho,1,N,peak_den,False,alt),exp_log_y(rho,1,N,peak_den,True,alt)] # = [x0/xm,...,xn/xm)]

        num = [np.cumsum(exp_log_y_num[0](n[idx_100])),np.cumsum(exp_log_y_num[1](n[idx]))] # sum{[x0/xm,...,xn/xm]}
        den = [np.cumsum(exp_log_y_den[0](n[idx_100])),np.cumsum(exp_log_y_den[1](n[idx]))] # sum{[x0,...,xn]} * (1/xm)
    
        if alt==True:
            keep[0].append( np.exp(-0.5*rho**2*(1-t**2)) * ((num[0][-1]+num[1][-1])/(den[0][-1]+den[1][-1])) * (np.exp(peak_num-peak_den)) )
        if alt==False:
            keep[1].append( np.exp(-0.5*rho**2*(1-t**2))  * ((num[0][-1]+num[1][-1])/(den[0][-1]+den[1][-1])) * (np.exp(peak_num-peak_den+(N-1)*np.log(t))) )

plt.plot(NN,keep[0],label="Eq53",linestyle="-",color="red")
plt.plot(NN,keep[1],label="Eq22",linestyle="--",color="blue")
plt.yscale('log')
plt.grid()
plt.xlabel('N (number of dimensions)')
plt.title(r"$\rho=$"+str(rho)+"$,(t_j\cdot t_k)=$"+str(t))
plt.legend()
plt.ylabel('Probability')
plt.savefig("rankingstatistic_3.pdf")
