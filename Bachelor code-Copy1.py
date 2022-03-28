#!/usr/bin/env python
# coding: utf-8

# In[295]:


import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as stats
from scipy.integrate import quad
from scipy.stats import logistic

#Parameters
beta=0.99
gamma=0.5
kappa=1.21
alpha=0.6
mu=0.43
bw=0.65
f=2.4
s=1.03
a_0=1
a_1 =0.5
c_f=0.23
phi = 0.03
phix = (2/3)*phi
#Steady state targets
eta_t=0.312
u_t=0.09
v_t=0.03971
theta_t = v_t / u_t


# In[299]:


#Defining wage level - equation 18)
def w_t(a_t,b_t):
    return gamma*a_t+(1-gamma)*b_t
w_0 = w_t(a_0, b_0)
w_1 = w_t(a_1, b_1)

#shock
def g(x):
    return logistic.pdf(x, 0, 2)

#Threshold for compensation - equation 3)
def v_k_t(a_t,w_t):
    return a_t-0.8*w_t
v_k_0 = v_k_t(a_0, w_0)
v_k_1 = v_k_t(a_1, w_1)


#Defining firing threshold in period 1
def v_f_period_1(w_t):
    return -(0.2*w_t-c_f)
v_f_1 = v_f_period_1(w_1)

#Equation (7)
#First line
def J_work(eps,a_t,w_t):
    return (a_t-w_t-eps)*g(eps)
#Integrating
def J_work_int(a_t,w_t,v_k_t):
    res, err = quad(J_work,-np.inf,v_k_t,args=(a_t,w_t))
    return res

#Second line
def J_comp(eps,w_t):
    return (-0.2*w_t-eps)*g(eps)
#Integrating
def J_comp_int(w_t, v_f_period_1, v_k_t):
    res, err = quad(J_comp,v_f_period_1, v_k_t,args=(w_t))
    return res

#Engogenous separation rate in t+1
def phie_t_1(v_f_period_1): 
    res, err = quad(g, v_f_period_1, np.inf)
    return res

phie_t_1 = phie_t_1(v_f_1)

phi_t_1= phix+(1-phix)*phie_t_1

first_line_integral = J_work_int(a_1,w_1,v_k_1)
second_line_integral = J_comp_int(w_1, v_f_1, v_k_1)
third_line_integral = (1-phi_t_1)*c_f-(1-phi)*phie_t_1*f

#Equation 7)
J_t_1 = (1-phi)* first_line_integral+ (1-phi)*second_line_integral-third_line_integral

#Equation 1)
#def J_eps(a_t, w_t, eps):
#    return a_t-w_t-eps-c_f+beta*J_t_1

#Threshold for firing period 0- Equation 4
def v_t_f(w_t):
    return -(0.2*w_t-c_f+beta*J_t_1)

v_f_0 = v_t_f(w_0)

#Unemployment benefits
def b_t(w_t):
    return 0.65*w_t
b_0= b_t(w_0)
b_1= b_t(w_1)

#Equation 5 - endogenous separation rate
def phie_t0(v_t_f):
    res, err = quad(g, v_t_f, np.inf)
    return res
phie_t_0 = phie_t0(v_f_0)

#Equation 6
chi_0 = res, err = quad(g, v_k_0, v_f_0)
chi_1 = res, err = quad(g, v_k_1, v_f_1)

#Equation 8
m_t = mu*(u_t**alpha)*(v_t**(1-alpha))

#Equatopm 9
q_t=mu*theta_t**(-alpha)

#Equatopm 10
eta_t2 = mu*theta_t**(1-alpha)

#Equation 12 (Tjek igennem)
kappa2 = beta*q_t*J_t_1

#Equation 13 (Tjek igennem. Skal gerne give 0.91)
n_t1 = 1-u_t
n_t = (1-phi)*n_t1+(1-phi)*eta_t*(1-n_t1)

#Equation 20
#available workers
abwk = n_t/(1-phie_t_0)
#First line
def abw(eps, a_t):
    return (a_t-eps)*g(eps)
#Integrating 
def abw_int(a_t, v_t_k):
    res, err = quad(abw, -np.inf , v_t_k ,args=(a_t))
    return res

#second Integrating 
def abw2_int(a_t, v_t_k, v_t_f):
    res, err = quad(abw2, v_t_k , v_t_f ,args=(a_t))
    return res

GDPfirst_t0 = abw_int(a_0, v_k_0)
GDPsecond_t0 = abw2_int(a_0, v_k_0, v_f_0)

#Equation 20
Y_0 = abwk*GDPfirst_t0+abwk*GDPsecond_t0-n_t*c_f-abwk*phi_t_1*f-v_t*kappa


# In[300]:


Y_0


# In[267]:


test

