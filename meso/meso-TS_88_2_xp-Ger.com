%chk=meso-TS_88_2_xp-12.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.100740  -0.436725   0.295815 
  6  -0.857979  -0.286324  -1.428905 
  1  -3.998989  -0.507198   0.898194 
  1   0.044923  -0.220426  -2.024218 
  6   1.256794   1.450448   0.922302 
  6   0.857979   0.286324   1.428905 
  1   0.679926   2.352386   1.089821 
  1  -0.044923   0.220426   2.024218 
  6   1.678016  -0.959021   1.247210 
  1   2.333558  -1.064055   2.121932 
  1   1.030749  -1.837588   1.242211 
  6   2.491654   1.537925   0.140607 
  1   2.894748   2.517041  -0.088961 
  6   3.100740   0.436725  -0.295815 
  1   3.998989   0.507198  -0.898194 
  6   2.520053  -0.925099  -0.029697 
  1   1.898127  -1.206584  -0.888225 
  1   3.312455  -1.673451   0.027707 
  6  -1.256794  -1.450448  -0.922302 
  1  -0.679926  -2.352386  -1.089821 
  6  -2.491654  -1.537925  -0.140607 
  1  -2.894748  -2.517041   0.088961 
  6  -1.678016   0.959021  -1.247210 
  1  -1.030749   1.837588  -1.242211 
  1  -2.333558   1.064055  -2.121932 
  6  -2.520053   0.925099   0.029697 
  1  -3.312455   1.673451  -0.027707 
  1  -1.898127   1.206584   0.888225 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.778099999999998

   -3.100740    -0.436725     0.295815     1.700000     1.200000
   -0.857979    -0.286324    -1.428905     1.700000     1.200000
   -3.998989    -0.507198     0.898194     1.200000     1.200000
    0.044923    -0.220426    -2.024218     1.200000     1.200000
    1.256794     1.450448     0.922302     1.700000     1.200000
    0.857979     0.286324     1.428905     1.700000     1.200000
    0.679926     2.352386     1.089821     1.200000     1.200000
   -0.044923     0.220426     2.024218     1.200000     1.200000
    1.678016    -0.959021     1.247210     1.700000     1.200000
    2.333558    -1.064055     2.121932     1.200000     1.200000
    1.030749    -1.837588     1.242211     1.200000     1.200000
    2.491654     1.537925     0.140607     1.700000     1.200000
    2.894748     2.517041    -0.088961     1.200000     1.200000
    3.100740     0.436725    -0.295815     1.700000     1.200000
    3.998989     0.507198    -0.898194     1.200000     1.200000
    2.520053    -0.925099    -0.029697     1.700000     1.200000
    1.898127    -1.206584    -0.888225     1.200000     1.200000
    3.312455    -1.673451     0.027707     1.200000     1.200000
   -1.256794    -1.450448    -0.922302     1.700000     1.200000
   -0.679926    -2.352386    -1.089821     1.200000     1.200000
   -2.491654    -1.537925    -0.140607     1.700000     1.200000
   -2.894748    -2.517041     0.088961     1.200000     1.200000
   -1.678016     0.959021    -1.247210     1.700000     1.200000
   -1.030749     1.837588    -1.242211     1.200000     1.200000
   -2.333558     1.064055    -2.121932     1.200000     1.200000
   -2.520053     0.925099     0.029697     1.700000     1.200000
   -3.312455     1.673451    -0.027707     1.200000     1.200000
   -1.898127     1.206584     0.888225     1.200000     1.200000

--link1--
%chk=meso-TS_88_2_xp-115.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.100740  -0.436725   0.295815 
  6  -0.857979  -0.286324  -1.428905 
  1  -3.998989  -0.507198   0.898194 
  1   0.044923  -0.220426  -2.024218 
  6   1.256794   1.450448   0.922302 
  6   0.857979   0.286324   1.428905 
  1   0.679926   2.352386   1.089821 
  1  -0.044923   0.220426   2.024218 
  6   1.678016  -0.959021   1.247210 
  1   2.333558  -1.064055   2.121932 
  1   1.030749  -1.837588   1.242211 
  6   2.491654   1.537925   0.140607 
  1   2.894748   2.517041  -0.088961 
  6   3.100740   0.436725  -0.295815 
  1   3.998989   0.507198  -0.898194 
  6   2.520053  -0.925099  -0.029697 
  1   1.898127  -1.206584  -0.888225 
  1   3.312455  -1.673451   0.027707 
  6  -1.256794  -1.450448  -0.922302 
  1  -0.679926  -2.352386  -1.089821 
  6  -2.491654  -1.537925  -0.140607 
  1  -2.894748  -2.517041   0.088961 
  6  -1.678016   0.959021  -1.247210 
  1  -1.030749   1.837588  -1.242211 
  1  -2.333558   1.064055  -2.121932 
  6  -2.520053   0.925099   0.029697 
  1  -3.312455   1.673451  -0.027707 
  1  -1.898127   1.206584   0.888225 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.08074995319988 RhoS=0.87957140140931

   -3.100740    -0.436725     0.295815     1.700000     1.150000
   -0.857979    -0.286324    -1.428905     1.700000     1.150000
   -3.998989    -0.507198     0.898194     1.200000     1.150000
    0.044923    -0.220426    -2.024218     1.200000     1.150000
    1.256794     1.450448     0.922302     1.700000     1.150000
    0.857979     0.286324     1.428905     1.700000     1.150000
    0.679926     2.352386     1.089821     1.200000     1.150000
   -0.044923     0.220426     2.024218     1.200000     1.150000
    1.678016    -0.959021     1.247210     1.700000     1.150000
    2.333558    -1.064055     2.121932     1.200000     1.150000
    1.030749    -1.837588     1.242211     1.200000     1.150000
    2.491654     1.537925     0.140607     1.700000     1.150000
    2.894748     2.517041    -0.088961     1.200000     1.150000
    3.100740     0.436725    -0.295815     1.700000     1.150000
    3.998989     0.507198    -0.898194     1.200000     1.150000
    2.520053    -0.925099    -0.029697     1.700000     1.150000
    1.898127    -1.206584    -0.888225     1.200000     1.150000
    3.312455    -1.673451     0.027707     1.200000     1.150000
   -1.256794    -1.450448    -0.922302     1.700000     1.150000
   -0.679926    -2.352386    -1.089821     1.200000     1.150000
   -2.491654    -1.537925    -0.140607     1.700000     1.150000
   -2.894748    -2.517041     0.088961     1.200000     1.150000
   -1.678016     0.959021    -1.247210     1.700000     1.150000
   -1.030749     1.837588    -1.242211     1.200000     1.150000
   -2.333558     1.064055    -2.121932     1.200000     1.150000
   -2.520053     0.925099     0.029697     1.700000     1.150000
   -3.312455     1.673451    -0.027707     1.200000     1.150000
   -1.898127     1.206584     0.888225     1.200000     1.150000

--link1--
%chk=meso-TS_88_2_xp-110.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.100740  -0.436725   0.295815 
  6  -0.857979  -0.286324  -1.428905 
  1  -3.998989  -0.507198   0.898194 
  1   0.044923  -0.220426  -2.024218 
  6   1.256794   1.450448   0.922302 
  6   0.857979   0.286324   1.428905 
  1   0.679926   2.352386   1.089821 
  1  -0.044923   0.220426   2.024218 
  6   1.678016  -0.959021   1.247210 
  1   2.333558  -1.064055   2.121932 
  1   1.030749  -1.837588   1.242211 
  6   2.491654   1.537925   0.140607 
  1   2.894748   2.517041  -0.088961 
  6   3.100740   0.436725  -0.295815 
  1   3.998989   0.507198  -0.898194 
  6   2.520053  -0.925099  -0.029697 
  1   1.898127  -1.206584  -0.888225 
  1   3.312455  -1.673451   0.027707 
  6  -1.256794  -1.450448  -0.922302 
  1  -0.679926  -2.352386  -1.089821 
  6  -2.491654  -1.537925  -0.140607 
  1  -2.894748  -2.517041   0.088961 
  6  -1.678016   0.959021  -1.247210 
  1  -1.030749   1.837588  -1.242211 
  1  -2.333558   1.064055  -2.121932 
  6  -2.520053   0.925099   0.029697 
  1  -3.312455   1.673451  -0.027707 
  1  -1.898127   1.206584   0.888225 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.14797859684703 RhoS=0.992403372535132

   -3.100740    -0.436725     0.295815     1.700000     1.100000
   -0.857979    -0.286324    -1.428905     1.700000     1.100000
   -3.998989    -0.507198     0.898194     1.200000     1.100000
    0.044923    -0.220426    -2.024218     1.200000     1.100000
    1.256794     1.450448     0.922302     1.700000     1.100000
    0.857979     0.286324     1.428905     1.700000     1.100000
    0.679926     2.352386     1.089821     1.200000     1.100000
   -0.044923     0.220426     2.024218     1.200000     1.100000
    1.678016    -0.959021     1.247210     1.700000     1.100000
    2.333558    -1.064055     2.121932     1.200000     1.100000
    1.030749    -1.837588     1.242211     1.200000     1.100000
    2.491654     1.537925     0.140607     1.700000     1.100000
    2.894748     2.517041    -0.088961     1.200000     1.100000
    3.100740     0.436725    -0.295815     1.700000     1.100000
    3.998989     0.507198    -0.898194     1.200000     1.100000
    2.520053    -0.925099    -0.029697     1.700000     1.100000
    1.898127    -1.206584    -0.888225     1.200000     1.100000
    3.312455    -1.673451     0.027707     1.200000     1.100000
   -1.256794    -1.450448    -0.922302     1.700000     1.100000
   -0.679926    -2.352386    -1.089821     1.200000     1.100000
   -2.491654    -1.537925    -0.140607     1.700000     1.100000
   -2.894748    -2.517041     0.088961     1.200000     1.100000
   -1.678016     0.959021    -1.247210     1.700000     1.100000
   -1.030749     1.837588    -1.242211     1.200000     1.100000
   -2.333558     1.064055    -2.121932     1.200000     1.100000
   -2.520053     0.925099     0.029697     1.700000     1.100000
   -3.312455     1.673451    -0.027707     1.200000     1.100000
   -1.898127     1.206584     0.888225     1.200000     1.100000

--link1--
%chk=meso-TS_88_2_xp-105.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.100740  -0.436725   0.295815 
  6  -0.857979  -0.286324  -1.428905 
  1  -3.998989  -0.507198   0.898194 
  1   0.044923  -0.220426  -2.024218 
  6   1.256794   1.450448   0.922302 
  6   0.857979   0.286324   1.428905 
  1   0.679926   2.352386   1.089821 
  1  -0.044923   0.220426   2.024218 
  6   1.678016  -0.959021   1.247210 
  1   2.333558  -1.064055   2.121932 
  1   1.030749  -1.837588   1.242211 
  6   2.491654   1.537925   0.140607 
  1   2.894748   2.517041  -0.088961 
  6   3.100740   0.436725  -0.295815 
  1   3.998989   0.507198  -0.898194 
  6   2.520053  -0.925099  -0.029697 
  1   1.898127  -1.206584  -0.888225 
  1   3.312455  -1.673451   0.027707 
  6  -1.256794  -1.450448  -0.922302 
  1  -0.679926  -2.352386  -1.089821 
  6  -2.491654  -1.537925  -0.140607 
  1  -2.894748  -2.517041   0.088961 
  6  -1.678016   0.959021  -1.247210 
  1  -1.030749   1.837588  -1.242211 
  1  -2.333558   1.064055  -2.121932 
  6  -2.520053   0.925099   0.029697 
  1  -3.312455   1.673451  -0.027707 
  1  -1.898127   1.206584   0.888225 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.22376880117094 RhoS=1.12776705486046

   -3.100740    -0.436725     0.295815     1.700000     1.050000
   -0.857979    -0.286324    -1.428905     1.700000     1.050000
   -3.998989    -0.507198     0.898194     1.200000     1.050000
    0.044923    -0.220426    -2.024218     1.200000     1.050000
    1.256794     1.450448     0.922302     1.700000     1.050000
    0.857979     0.286324     1.428905     1.700000     1.050000
    0.679926     2.352386     1.089821     1.200000     1.050000
   -0.044923     0.220426     2.024218     1.200000     1.050000
    1.678016    -0.959021     1.247210     1.700000     1.050000
    2.333558    -1.064055     2.121932     1.200000     1.050000
    1.030749    -1.837588     1.242211     1.200000     1.050000
    2.491654     1.537925     0.140607     1.700000     1.050000
    2.894748     2.517041    -0.088961     1.200000     1.050000
    3.100740     0.436725    -0.295815     1.700000     1.050000
    3.998989     0.507198    -0.898194     1.200000     1.050000
    2.520053    -0.925099    -0.029697     1.700000     1.050000
    1.898127    -1.206584    -0.888225     1.200000     1.050000
    3.312455    -1.673451     0.027707     1.200000     1.050000
   -1.256794    -1.450448    -0.922302     1.700000     1.050000
   -0.679926    -2.352386    -1.089821     1.200000     1.050000
   -2.491654    -1.537925    -0.140607     1.700000     1.050000
   -2.894748    -2.517041     0.088961     1.200000     1.050000
   -1.678016     0.959021    -1.247210     1.700000     1.050000
   -1.030749     1.837588    -1.242211     1.200000     1.050000
   -2.333558     1.064055    -2.121932     1.200000     1.050000
   -2.520053     0.925099     0.029697     1.700000     1.050000
   -3.312455     1.673451    -0.027707     1.200000     1.050000
   -1.898127     1.206584     0.888225     1.200000     1.050000

--link1--
%chk=meso-TS_88_2_xp-100.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.100740  -0.436725   0.295815 
  6  -0.857979  -0.286324  -1.428905 
  1  -3.998989  -0.507198   0.898194 
  1   0.044923  -0.220426  -2.024218 
  6   1.256794   1.450448   0.922302 
  6   0.857979   0.286324   1.428905 
  1   0.679926   2.352386   1.089821 
  1  -0.044923   0.220426   2.024218 
  6   1.678016  -0.959021   1.247210 
  1   2.333558  -1.064055   2.121932 
  1   1.030749  -1.837588   1.242211 
  6   2.491654   1.537925   0.140607 
  1   2.894748   2.517041  -0.088961 
  6   3.100740   0.436725  -0.295815 
  1   3.998989   0.507198  -0.898194 
  6   2.520053  -0.925099  -0.029697 
  1   1.898127  -1.206584  -0.888225 
  1   3.312455  -1.673451   0.027707 
  6  -1.256794  -1.450448  -0.922302 
  1  -0.679926  -2.352386  -1.089821 
  6  -2.491654  -1.537925  -0.140607 
  1  -2.894748  -2.517041   0.088961 
  6  -1.678016   0.959021  -1.247210 
  1  -1.030749   1.837588  -1.242211 
  1  -2.333558   1.064055  -2.121932 
  6  -2.520053   0.925099   0.029697 
  1  -3.312455   1.673451  -0.027707 
  1  -1.898127   1.206584   0.888225 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.30480546666792 RhoS=1.28207131826017

   -3.100740    -0.436725     0.295815     1.700000     1.000000
   -0.857979    -0.286324    -1.428905     1.700000     1.000000
   -3.998989    -0.507198     0.898194     1.200000     1.000000
    0.044923    -0.220426    -2.024218     1.200000     1.000000
    1.256794     1.450448     0.922302     1.700000     1.000000
    0.857979     0.286324     1.428905     1.700000     1.000000
    0.679926     2.352386     1.089821     1.200000     1.000000
   -0.044923     0.220426     2.024218     1.200000     1.000000
    1.678016    -0.959021     1.247210     1.700000     1.000000
    2.333558    -1.064055     2.121932     1.200000     1.000000
    1.030749    -1.837588     1.242211     1.200000     1.000000
    2.491654     1.537925     0.140607     1.700000     1.000000
    2.894748     2.517041    -0.088961     1.200000     1.000000
    3.100740     0.436725    -0.295815     1.700000     1.000000
    3.998989     0.507198    -0.898194     1.200000     1.000000
    2.520053    -0.925099    -0.029697     1.700000     1.000000
    1.898127    -1.206584    -0.888225     1.200000     1.000000
    3.312455    -1.673451     0.027707     1.200000     1.000000
   -1.256794    -1.450448    -0.922302     1.700000     1.000000
   -0.679926    -2.352386    -1.089821     1.200000     1.000000
   -2.491654    -1.537925    -0.140607     1.700000     1.000000
   -2.894748    -2.517041     0.088961     1.200000     1.000000
   -1.678016     0.959021    -1.247210     1.700000     1.000000
   -1.030749     1.837588    -1.242211     1.200000     1.000000
   -2.333558     1.064055    -2.121932     1.200000     1.000000
   -2.520053     0.925099     0.029697     1.700000     1.000000
   -3.312455     1.673451    -0.027707     1.200000     1.000000
   -1.898127     1.206584     0.888225     1.200000     1.000000

--link1--
%chk=meso-TS_88_2_xp-0975.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.100740  -0.436725   0.295815 
  6  -0.857979  -0.286324  -1.428905 
  1  -3.998989  -0.507198   0.898194 
  1   0.044923  -0.220426  -2.024218 
  6   1.256794   1.450448   0.922302 
  6   0.857979   0.286324   1.428905 
  1   0.679926   2.352386   1.089821 
  1  -0.044923   0.220426   2.024218 
  6   1.678016  -0.959021   1.247210 
  1   2.333558  -1.064055   2.121932 
  1   1.030749  -1.837588   1.242211 
  6   2.491654   1.537925   0.140607 
  1   2.894748   2.517041  -0.088961 
  6   3.100740   0.436725  -0.295815 
  1   3.998989   0.507198  -0.898194 
  6   2.520053  -0.925099  -0.029697 
  1   1.898127  -1.206584  -0.888225 
  1   3.312455  -1.673451   0.027707 
  6  -1.256794  -1.450448  -0.922302 
  1  -0.679926  -2.352386  -1.089821 
  6  -2.491654  -1.537925  -0.140607 
  1  -2.894748  -2.517041   0.088961 
  6  -1.678016   0.959021  -1.247210 
  1  -1.030749   1.837588  -1.242211 
  1  -2.333558   1.064055  -2.121932 
  6  -2.520053   0.925099   0.029697 
  1  -3.312455   1.673451  -0.027707 
  1  -1.898127   1.206584   0.888225 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.34745273026226 RhoS=1.36724935387803

   -3.100740    -0.436725     0.295815     1.700000     0.975000
   -0.857979    -0.286324    -1.428905     1.700000     0.975000
   -3.998989    -0.507198     0.898194     1.200000     0.975000
    0.044923    -0.220426    -2.024218     1.200000     0.975000
    1.256794     1.450448     0.922302     1.700000     0.975000
    0.857979     0.286324     1.428905     1.700000     0.975000
    0.679926     2.352386     1.089821     1.200000     0.975000
   -0.044923     0.220426     2.024218     1.200000     0.975000
    1.678016    -0.959021     1.247210     1.700000     0.975000
    2.333558    -1.064055     2.121932     1.200000     0.975000
    1.030749    -1.837588     1.242211     1.200000     0.975000
    2.491654     1.537925     0.140607     1.700000     0.975000
    2.894748     2.517041    -0.088961     1.200000     0.975000
    3.100740     0.436725    -0.295815     1.700000     0.975000
    3.998989     0.507198    -0.898194     1.200000     0.975000
    2.520053    -0.925099    -0.029697     1.700000     0.975000
    1.898127    -1.206584    -0.888225     1.200000     0.975000
    3.312455    -1.673451     0.027707     1.200000     0.975000
   -1.256794    -1.450448    -0.922302     1.700000     0.975000
   -0.679926    -2.352386    -1.089821     1.200000     0.975000
   -2.491654    -1.537925    -0.140607     1.700000     0.975000
   -2.894748    -2.517041     0.088961     1.200000     0.975000
   -1.678016     0.959021    -1.247210     1.700000     0.975000
   -1.030749     1.837588    -1.242211     1.200000     0.975000
   -2.333558     1.064055    -2.121932     1.200000     0.975000
   -2.520053     0.925099     0.029697     1.700000     0.975000
   -3.312455     1.673451    -0.027707     1.200000     0.975000
   -1.898127     1.206584     0.888225     1.200000     0.975000

--link1--
%chk=meso-TS_88_2_xp-095.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.100740  -0.436725   0.295815 
  6  -0.857979  -0.286324  -1.428905 
  1  -3.998989  -0.507198   0.898194 
  1   0.044923  -0.220426  -2.024218 
  6   1.256794   1.450448   0.922302 
  6   0.857979   0.286324   1.428905 
  1   0.679926   2.352386   1.089821 
  1  -0.044923   0.220426   2.024218 
  6   1.678016  -0.959021   1.247210 
  1   2.333558  -1.064055   2.121932 
  1   1.030749  -1.837588   1.242211 
  6   2.491654   1.537925   0.140607 
  1   2.894748   2.517041  -0.088961 
  6   3.100740   0.436725  -0.295815 
  1   3.998989   0.507198  -0.898194 
  6   2.520053  -0.925099  -0.029697 
  1   1.898127  -1.206584  -0.888225 
  1   3.312455  -1.673451   0.027707 
  6  -1.256794  -1.450448  -0.922302 
  1  -0.679926  -2.352386  -1.089821 
  6  -2.491654  -1.537925  -0.140607 
  1  -2.894748  -2.517041   0.088961 
  6  -1.678016   0.959021  -1.247210 
  1  -1.030749   1.837588  -1.242211 
  1  -2.333558   1.064055  -2.121932 
  6  -2.520053   0.925099   0.029697 
  1  -3.312455   1.673451  -0.027707 
  1  -1.898127   1.206584   0.888225 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.39339884323873 RhoS=1.4620813611373

   -3.100740    -0.436725     0.295815     1.700000     0.950000
   -0.857979    -0.286324    -1.428905     1.700000     0.950000
   -3.998989    -0.507198     0.898194     1.200000     0.950000
    0.044923    -0.220426    -2.024218     1.200000     0.950000
    1.256794     1.450448     0.922302     1.700000     0.950000
    0.857979     0.286324     1.428905     1.700000     0.950000
    0.679926     2.352386     1.089821     1.200000     0.950000
   -0.044923     0.220426     2.024218     1.200000     0.950000
    1.678016    -0.959021     1.247210     1.700000     0.950000
    2.333558    -1.064055     2.121932     1.200000     0.950000
    1.030749    -1.837588     1.242211     1.200000     0.950000
    2.491654     1.537925     0.140607     1.700000     0.950000
    2.894748     2.517041    -0.088961     1.200000     0.950000
    3.100740     0.436725    -0.295815     1.700000     0.950000
    3.998989     0.507198    -0.898194     1.200000     0.950000
    2.520053    -0.925099    -0.029697     1.700000     0.950000
    1.898127    -1.206584    -0.888225     1.200000     0.950000
    3.312455    -1.673451     0.027707     1.200000     0.950000
   -1.256794    -1.450448    -0.922302     1.700000     0.950000
   -0.679926    -2.352386    -1.089821     1.200000     0.950000
   -2.491654    -1.537925    -0.140607     1.700000     0.950000
   -2.894748    -2.517041     0.088961     1.200000     0.950000
   -1.678016     0.959021    -1.247210     1.700000     0.950000
   -1.030749     1.837588    -1.242211     1.200000     0.950000
   -2.333558     1.064055    -2.121932     1.200000     0.950000
   -2.520053     0.925099     0.029697     1.700000     0.950000
   -3.312455     1.673451    -0.027707     1.200000     0.950000
   -1.898127     1.206584     0.888225     1.200000     0.950000

