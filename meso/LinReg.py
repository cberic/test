#!/usr/bin/env python
#from numpy import sqrt, pi, exp, linspace, loadtxt, savetxt
import numpy as np
#from sklearn.linear_model import LinearRegression
from scipy import stats
from numpy import sqrt, pi, exp, linspace, loadtxt, savetxt
data = loadtxt('./pGtot.dat')
#x = data[:, 0].reshape((-1, 1))
x = data[:, 0]
y = data[:, 1]
print x
print y
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print('Slope: ',slope,'\nIntercept: ',intercept)

