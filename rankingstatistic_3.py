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

def log_y(rho,t,N,approx):
    if approx==False:
        return lambda x: -0.5*rho**2*(1-t**2) + (gammaln(N)-(gammaln(x+1)+gammaln(N-x))) + x*np.log(rho*t/np.sqrt(2)) + np.log(incomplete_Gamma((N-x)/2,0.5*(rho*t)**2))
    if approx==True:
        return lambda x: -0.5*rho**2*(1-t**2) + (gammaln(N)-(gammaln(x+1)+gammaln(N-x))) + x*np.log(rho*t/np.sqrt(2)) + gammaln((N-x)/2)

def exp_log_y(rho,t,N,log_peak,approx):
    if approx==False:
        return lambda x: np.exp(-0.5*rho**2*(1-t**2) + gammaln(N)-(gammaln(x+1)+gammaln(N-x)) + x*np.log(rho*t/np.sqrt(2)) + np.log(incomplete_Gamma((N-x)/2,0.5*(rho*t)**2)) - log_peak)
    if approx==True:
        return lambda x: np.exp(-0.5*rho**2*(1-t**2) + gammaln(N)-(gammaln(x+1)+gammaln(N-x)) + x*np.log(rho*t/np.sqrt(2)) + gammaln((N-x)/2) - log_peak)

NN = np.linspace(1000,10000,50) # number of dimensions
rho = 10. # nominal SNR
t = 0.9 # (t_j \cdot t_k)

keep = []

for N in NN:

    n = np.arange(0,N-1) # range of n in summation
    idx_100 = np.where((N-n) <= 100) # indices where (N-n) <= 100 (i.e. range where Gamma(a,x) != Gamma(a) ).
    idx = np.where((N-n) > 100) # indices where (N-n) > 100 (i.e. range where Gamma(a,x) = Gamma(a) ).

    peak_num = max(log_y(rho,t,N,True)(n[idx]))
    peak_den = max(log_y(rho,1,N,True)(n[idx]))
                   
    exp_log_y_num = [exp_log_y(rho,t,N,peak_num,False),exp_log_y(rho,t,N,peak_num,True)] # [exp(log(x0/xm)),...,exp(log(xn/xm))]
    exp_log_y_den = [exp_log_y(rho,1,N,peak_den,False),exp_log_y(rho,1,N,peak_den,True)] # = [x0/xm,...,xn/xm)]

    num = [np.cumsum(exp_log_y_num[0](n[idx_100])),np.cumsum(exp_log_y_num[1](n[idx]))] # sum{[x0/xm,...,xn/xm]}
    den = [np.cumsum(exp_log_y_den[0](n[idx_100])),np.cumsum(exp_log_y_den[1](n[idx]))] # sum{[x0,...,xn]} * (1/xm)
    
    keep.append( ((num[0][-1]+num[1][-1])/(den[0][-1]+den[1][-1])) * (np.exp(peak_num-peak_den)) )

plt.plot(NN,keep,label=r"$\rho=$"+str(rho)+"$,(t_j\cdot t_k)=$"+str(t))
plt.yscale('log')
plt.xlabel('N (number of dimensions)')
plt.legend()
plt.ylabel('Probability')
plt.savefig("rankingstatistic_3.pdf")
