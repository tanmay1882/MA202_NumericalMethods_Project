'''
s = stock price
tau = time to expiration
K = exercise price
'''

import math as m
import numpy as np
import matplotlib.pyplot as plt
global sigma, r, K, T, xL, xH

def BSPut(s, tau, K):
    d1 = []
    d2 = []
    p = []
    for i in range(len(s)):
        d1.append((m.log(s[i]/K)+(r+sigma**2/2)*tau)/(sigma*m.sqrt(tau)))
    for i in range(len(s)):
        d2.append(d1[i]-(sigma*m.sqrt(tau)))
    for i in range(len(s)):
        p.append((-s[i]*m.erfc(d1[i]/m.sqrt(2))+K*m.exp(-r*tau)*m.erfc(d2[i]/m.sqrt(2)))/2)
    return p

def fh(u,x,tau,T):
    mu = (((sigma/2)-(r/sigma))/(m.sqrt(2)))
    nu = r + mu**2
    t = T - tau
    s = []
    v = []
    for i in range(len(x)):
        s.append(m.exp((sigma*x[i])/m.sqrt(2)))
        v.append(m.exp(mu*x[i]-nu*tau)*(u))
    return s, v

def hf(v,s,t,T):
    mu = (((sigma/2)-(r/sigma))/(m.sqrt(2)))
    nu = r + mu**2
    tau = T - t
    x = []
    u = []
    if type(v)==list:
        for i in range(len(s)):
            x.append(m.sqrt(2)*m.log(s[i])/sigma)
            u.append(v[i]/m.exp(mu*x[i]-nu*tau))
    else:
        for i in range(len(s)):
            x.append(m.sqrt(2)*m.log(s[i])/sigma)
            u.append(v/m.exp(mu*x[i]-nu*tau))
    return x, u

def u0(x):
    s, trash = fh(0,x,0,0)
    v = vT(s)
    x, u = hf(v,s,0,0)
    return(u)

def uH(x, tau):
    s, v = fh(0,x,tau,tau)
    v = vH(s,0,tau)
    x, u = hf(v,s,0,tau)
    return u

def uL(x, tau):
    s, v = fh(0,x,tau,T)
    v = vL(s,0,tau)
    x, u = hf(v,s,0,tau)
    return u

def vH(s,t,T):
    v = 0
    return v

def vL(s,t,T):
    v = []
    for i in range(len(s)):
        v.append(K*m.exp(r*(T-t))-s[i])
    return v

def vT(s):
    v_mod = []
    for i in range(len(s)):
        v_mod.append((K-s[i]))
    for i in range(len(v_mod)):
        if(v_mod[i]<0):
            v_mod[i]=0
    return v_mod


Stock_lower = 50
Stock_higher = 150
K = 100
sigma = 0.2
r = 0.01
T = 1


x, u = hf(0,[Stock_lower, Stock_higher], 0, 0)
xL = x[0]
xH = x[1]

N = 100
M = 1000
dx = (xH-xL)/(N+1)
dt = T/M
alpha = dt/(dx)**2
x = np.linspace(xL, xH, N+2)
u = []
u = u0(x[1:N])
print(u)
tau = 0

u1 = u[1:N]
u2 = u[0:N-1]

for i in range(1,M):
    uH_value = alpha*uH([xH],tau)[0]
    uL_value = alpha*uL([xL],tau)[0]
    tau = i*dt


u_mod = []
u1_mod = []
u2_mod = []

for j in range(1,M):
        for i in range(len(u1)):
            if (0<i<M):
                u[i] = u1[i]*alpha+u2[i]*alpha+(1-(2*alpha))*u[i]
            
            elif (i==M):
                u[i] = uH([xH],tau)[0]+u2[i]*alpha+(1-(2*alpha))*u[i]
                

u_master = []
for i in range(len(u)):
    u_master.append(u[i]) 

u_master.insert(0,uL([xL],T)[0])
u_master.append(uH([xH],T)[0])

s, trash = fh(0,x,0,0)
true_values = hf(BSPut(s,T,K),s,0,T)

x_master = x[0:len(x)-1]
u.append(0)
u.append(0)

print("x:\n")
print(x_master)
print("u = \n")
print(u_master)

plt.plot(x_master,u_master,"*")
plt.xlabel('x - axis')
plt.ylabel('y - axis')
plt.show()
