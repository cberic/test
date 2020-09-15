%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.118578  -0.436876   0.247024 
  6  -0.555477  -0.196368  -0.932206 
  1  -4.082851  -0.549917   0.728725 
  1   0.380936  -0.117306  -1.466881 
  6   1.090379   1.404463   0.689379 
  6   0.555477   0.196368   0.932206 
  1   0.521599   2.308823   0.864348 
  1  -0.380936   0.117306   1.466881 
  6   1.466810  -0.999509   0.989932 
  1   1.911942  -0.970218   1.991724 
  1   0.912395  -1.934620   0.924846 
  6   2.418004   1.516122   0.112667 
  1   2.827281   2.505496  -0.054615 
  6   3.118578   0.436876  -0.247024 
  1   4.082851   0.549917  -0.728725 
  6   2.591691  -0.957992  -0.044788 
  1   2.252888  -1.355863  -1.009385 
  1   3.400558  -1.620027   0.271880 
  6  -1.090379  -1.404463  -0.689379 
  1  -0.521599  -2.308823  -0.864348 
  6  -2.418004  -1.516122  -0.112667 
  1  -2.827281  -2.505496   0.054615 
  6  -1.466810   0.999509  -0.989932 
  1  -0.912395   1.934620  -0.924846 
  1  -1.911942   0.970219  -1.991724 
  6  -2.591691   0.957992   0.044788 
  1  -3.400558   1.620027  -0.271880 
  1  -2.252888   1.355863   1.009385 

norep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

    1     1.700000     1.200000
    2     1.700000     1.200000
    3     1.200000     1.200000
    4     1.200000     1.200000
    5     1.700000     1.200000
    6     1.700000     1.200000
    7     1.200000     1.200000
    8     1.200000     1.200000
    9     1.700000     1.200000
   10     1.200000     1.200000
   11     1.200000     1.200000
   12     1.700000     1.200000
   13     1.200000     1.200000
   14     1.700000     1.200000
   15     1.200000     1.200000
   16     1.700000     1.200000
   17     1.200000     1.200000
   18     1.200000     1.200000
   19     1.700000     1.200000
   20     1.200000     1.200000
   21     1.700000     1.200000
   22     1.200000     1.200000
   23     1.700000     1.200000
   24     1.200000     1.200000
   25     1.200000     1.200000
   26     1.700000     1.200000
   27     1.200000     1.200000
   28     1.200000     1.200000

