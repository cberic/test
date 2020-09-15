###!/usr/local/bin/WolframScript -script

(* generate high-precision samples of a mixed distribution *)


t = Import["./VGer.dat", "Table"]

Print[t]


<< NonlinearRegression`

Print[NonlinearRegress[t, 
 t[[1, 2]] + a*x ((1/b)*(t[[1, 1]]/x)^(b + 1) + 1) - c*x, {a, b, 
  c}, x, RegressionReport -> {BestFitParameters, ParameterCITable, 
   ANOVATable, BestFit, PredictedResponse}]]

