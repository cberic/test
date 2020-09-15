%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.126543  -0.441233   0.268190 
  6  -0.754024  -0.256533  -1.270411 
  1  -4.062711  -0.529769   0.807219 
  1   0.190409  -0.176659  -1.794641 
  6   1.183128   1.429264   0.813376 
  6   0.754024   0.256533   1.270411 
  1   0.591290   2.326044   0.953170 
  1  -0.190409   0.176659   1.794641 
  6   1.584536  -0.987345   1.138633 
  1   2.135575  -1.119037   2.078702 
  1   0.936296  -1.859379   1.043739 
  6   2.464520   1.531740   0.114055 
  1   2.863784   2.516400  -0.097950 
  6   3.126543   0.441233  -0.268190 
  1   4.062711   0.529769  -0.807219 
  6   2.569963  -0.936686  -0.033648 
  1   2.075826  -1.268120  -0.955200 
  1   3.381182  -1.646980   0.137036 
  6  -1.183128  -1.429264  -0.813376 
  1  -0.591290  -2.326044  -0.953170 
  6  -2.464520  -1.531740  -0.114055 
  1  -2.863784  -2.516400   0.097950 
  6  -1.584536   0.987345  -1.138633 
  1  -0.936296   1.859379  -1.043739 
  1  -2.135575   1.119037  -2.078702 
  6  -2.569963   0.936686   0.033648 
  1  -3.381182   1.646980  -0.137036 
  1  -2.075826   1.268120   0.955200 

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

