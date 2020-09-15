###!/usr/local/bin/WolframScript -script
t = Import["./pGtot.dat", "Table"]
model = LinearModelFit[t, x, x]
Print[model["BestFit"]]


