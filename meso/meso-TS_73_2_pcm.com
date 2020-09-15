%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120165  -0.439520   0.280190 
  6  -0.807998  -0.270483  -1.349166 
  1  -4.039260  -0.519446   0.849052 
  1   0.117969  -0.196755  -1.906627 
  6   1.221266   1.438778   0.865571 
  6   0.807998   0.270483   1.349166 
  1   0.634256   2.337167   1.015854 
  1  -0.117969   0.196755   1.906627 
  6   1.635901  -0.972988   1.195334 
  1   2.239675  -1.089626   2.104856 
  1   0.987632  -1.849296   1.145396 
  6   2.482134   1.535015   0.127997 
  1   2.882920   2.517315  -0.091872 
  6   3.120165   0.439520  -0.280190 
  1   4.039260   0.519446  -0.849052 
  6   2.551272  -0.930669  -0.031151 
  1   1.991682  -1.236938  -0.923478 
  1   3.354963  -1.660578   0.081874 
  6  -1.221266  -1.438778  -0.865571 
  1  -0.634256  -2.337167  -1.015854 
  6  -2.482134  -1.535015  -0.127997 
  1  -2.882920  -2.517315   0.091872 
  6  -1.635901   0.972988  -1.195334 
  1  -0.987632   1.849296  -1.145396 
  1  -2.239675   1.089626  -2.104856 
  6  -2.551272   0.930669   0.031151 
  1  -3.354963   1.660578  -0.081874 
  1  -1.991682   1.236938   0.923479 

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

