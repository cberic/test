%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.119691  -0.425826   0.248428 
  6  -0.438059  -0.127566  -0.740444 
  1  -4.086791  -0.546863   0.720948 
  1   0.329828  -0.109733  -1.515278 
  6   1.116208   1.432588   0.668787 
  6   0.438059   0.127566   0.740444 
  1   0.567375   2.333262   0.913532 
  1  -0.329828   0.109733   1.515278 
  6   1.453963  -0.997822   0.978822 
  1   1.872621  -0.887919   1.983401 
  1   0.950096  -1.965338   0.950002 
  6   2.395835   1.529761   0.125114 
  1   2.827943   2.513325  -0.022973 
  6   3.119691   0.425826  -0.248428 
  1   4.086791   0.546863  -0.720948 
  6   2.592903  -0.959744  -0.042391 
  1   2.246941  -1.360164  -1.004660 
  1   3.401404  -1.625003   0.272285 
  6  -1.116208  -1.432588  -0.668787 
  1  -0.567375  -2.333262  -0.913532 
  6  -2.395835  -1.529761  -0.125114 
  1  -2.827943  -2.513325   0.022973 
  6  -1.453963   0.997822  -0.978822 
  1  -0.950096   1.965338  -0.950002 
  1  -1.872621   0.887919  -1.983401 
  6  -2.592903   0.959744   0.042391 
  1  -3.401404   1.625003  -0.272284 
  1  -2.246941   1.360164   1.004660 

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

