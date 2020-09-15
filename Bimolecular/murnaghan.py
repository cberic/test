#!/usr/bin/env python
#<examples/doc_model1.py>
from numpy import sqrt, pi, exp, linspace, loadtxt, savetxt
from lmfit import  Model

### import matplotlib.pyplot as plt

data = loadtxt('VGer.dat')
x = data[:, 0]/data[0,0]
y = data[:, 1]-data[0,1]

def murnaghan(x, a, b, c):
    "1-d gaussian: murnaghan(x, a,b, c)"
    return (a/b)*(1/x)**b+(a-c)*x

gmod = Model(murnaghan)
result = gmod.fit(y, x=x, a=0, b=5, c=0)

print(result.fit_report())


#<end examples/doc_model1.py>
