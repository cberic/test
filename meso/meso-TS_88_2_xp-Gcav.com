%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=108.160262177098 cm^3/mol 

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

norep  pcmdoc geomview  nodis cav g03defaults noaddsph  tsare=0.075 
nsfe=28
Vmol=108.160262177098 Rsolv=2.815 

    1     2.040000     1.000000
    2     2.040000     1.000000
    3     1.440000     1.000000
    4     1.440000     1.000000
    5     2.040000     1.000000
    6     2.040000     1.000000
    7     1.440000     1.000000
    8     1.440000     1.000000
    9     2.040000     1.000000
   10     1.440000     1.000000
   11     1.440000     1.000000
   12     2.040000     1.000000
   13     1.440000     1.000000
   14     2.040000     1.000000
   15     1.440000     1.000000
   16     2.040000     1.000000
   17     1.440000     1.000000
   18     1.440000     1.000000
   19     2.040000     1.000000
   20     1.440000     1.000000
   21     2.040000     1.000000
   22     1.440000     1.000000
   23     2.040000     1.000000
   24     1.440000     1.000000
   25     1.440000     1.000000
   26     2.040000     1.000000
   27     1.440000     1.000000
   28     1.440000     1.000000

--link1--
%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=101.730197792279 cm^3/mol 

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

norep  pcmdoc geomview  nodis cav g03defaults noaddsph  tsare=0.075 
nsfe=28
Vmol=101.730197792279 Rsolv=2.815 

    1     2.040000     1.000000
    2     2.040000     1.000000
    3     1.440000     1.000000
    4     1.440000     1.000000
    5     2.040000     1.000000
    6     2.040000     1.000000
    7     1.440000     1.000000
    8     1.440000     1.000000
    9     2.040000     1.000000
   10     1.440000     1.000000
   11     1.440000     1.000000
   12     2.040000     1.000000
   13     1.440000     1.000000
   14     2.040000     1.000000
   15     1.440000     1.000000
   16     2.040000     1.000000
   17     1.440000     1.000000
   18     1.440000     1.000000
   19     2.040000     1.000000
   20     1.440000     1.000000
   21     2.040000     1.000000
   22     1.440000     1.000000
   23     2.040000     1.000000
   24     1.440000     1.000000
   25     1.440000     1.000000
   26     2.040000     1.000000
   27     1.440000     1.000000
   28     1.440000     1.000000

--link1--
%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=95.7726100512575 cm^3/mol 

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

norep  pcmdoc geomview  nodis cav g03defaults noaddsph  tsare=0.075 
nsfe=28
Vmol=95.7726100512575 Rsolv=2.815 

    1     2.040000     1.000000
    2     2.040000     1.000000
    3     1.440000     1.000000
    4     1.440000     1.000000
    5     2.040000     1.000000
    6     2.040000     1.000000
    7     1.440000     1.000000
    8     1.440000     1.000000
    9     2.040000     1.000000
   10     1.440000     1.000000
   11     1.440000     1.000000
   12     2.040000     1.000000
   13     1.440000     1.000000
   14     2.040000     1.000000
   15     1.440000     1.000000
   16     2.040000     1.000000
   17     1.440000     1.000000
   18     1.440000     1.000000
   19     2.040000     1.000000
   20     1.440000     1.000000
   21     2.040000     1.000000
   22     1.440000     1.000000
   23     2.040000     1.000000
   24     1.440000     1.000000
   25     1.440000     1.000000
   26     2.040000     1.000000
   27     1.440000     1.000000
   28     1.440000     1.000000

--link1--
%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=89.84123994485 cm^3/mol 

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

norep  pcmdoc geomview  nodis cav g03defaults noaddsph  tsare=0.075 
nsfe=28
Vmol=89.84123994485 Rsolv=2.815 

    1     2.040000     1.000000
    2     2.040000     1.000000
    3     1.440000     1.000000
    4     1.440000     1.000000
    5     2.040000     1.000000
    6     2.040000     1.000000
    7     1.440000     1.000000
    8     1.440000     1.000000
    9     2.040000     1.000000
   10     1.440000     1.000000
   11     1.440000     1.000000
   12     2.040000     1.000000
   13     1.440000     1.000000
   14     2.040000     1.000000
   15     1.440000     1.000000
   16     2.040000     1.000000
   17     1.440000     1.000000
   18     1.440000     1.000000
   19     2.040000     1.000000
   20     1.440000     1.000000
   21     2.040000     1.000000
   22     1.440000     1.000000
   23     2.040000     1.000000
   24     1.440000     1.000000
   25     1.440000     1.000000
   26     2.040000     1.000000
   27     1.440000     1.000000
   28     1.440000     1.000000

--link1--
%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=84.2615311719887 cm^3/mol 

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

norep  pcmdoc geomview  nodis cav g03defaults noaddsph  tsare=0.075 
nsfe=28
Vmol=84.2615311719887 Rsolv=2.815 

    1     2.040000     1.000000
    2     2.040000     1.000000
    3     1.440000     1.000000
    4     1.440000     1.000000
    5     2.040000     1.000000
    6     2.040000     1.000000
    7     1.440000     1.000000
    8     1.440000     1.000000
    9     2.040000     1.000000
   10     1.440000     1.000000
   11     1.440000     1.000000
   12     2.040000     1.000000
   13     1.440000     1.000000
   14     2.040000     1.000000
   15     1.440000     1.000000
   16     2.040000     1.000000
   17     1.440000     1.000000
   18     1.440000     1.000000
   19     2.040000     1.000000
   20     1.440000     1.000000
   21     2.040000     1.000000
   22     1.440000     1.000000
   23     2.040000     1.000000
   24     1.440000     1.000000
   25     1.440000     1.000000
   26     2.040000     1.000000
   27     1.440000     1.000000
   28     1.440000     1.000000

--link1--
%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=81.5946296547421 cm^3/mol 

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

norep  pcmdoc geomview  nodis cav g03defaults noaddsph  tsare=0.075 
nsfe=28
Vmol=81.5946296547421 Rsolv=2.815 

    1     2.040000     1.000000
    2     2.040000     1.000000
    3     1.440000     1.000000
    4     1.440000     1.000000
    5     2.040000     1.000000
    6     2.040000     1.000000
    7     1.440000     1.000000
    8     1.440000     1.000000
    9     2.040000     1.000000
   10     1.440000     1.000000
   11     1.440000     1.000000
   12     2.040000     1.000000
   13     1.440000     1.000000
   14     2.040000     1.000000
   15     1.440000     1.000000
   16     2.040000     1.000000
   17     1.440000     1.000000
   18     1.440000     1.000000
   19     2.040000     1.000000
   20     1.440000     1.000000
   21     2.040000     1.000000
   22     1.440000     1.000000
   23     2.040000     1.000000
   24     1.440000     1.000000
   25     1.440000     1.000000
   26     2.040000     1.000000
   27     1.440000     1.000000
   28     1.440000     1.000000

--link1--
%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=78.9041178242053 cm^3/mol 

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

norep  pcmdoc geomview  nodis cav g03defaults noaddsph  tsare=0.075 
nsfe=28
Vmol=78.9041178242053 Rsolv=2.815 

    1     2.040000     1.000000
    2     2.040000     1.000000
    3     1.440000     1.000000
    4     1.440000     1.000000
    5     2.040000     1.000000
    6     2.040000     1.000000
    7     1.440000     1.000000
    8     1.440000     1.000000
    9     2.040000     1.000000
   10     1.440000     1.000000
   11     1.440000     1.000000
   12     2.040000     1.000000
   13     1.440000     1.000000
   14     2.040000     1.000000
   15     1.440000     1.000000
   16     2.040000     1.000000
   17     1.440000     1.000000
   18     1.440000     1.000000
   19     2.040000     1.000000
   20     1.440000     1.000000
   21     2.040000     1.000000
   22     1.440000     1.000000
   23     2.040000     1.000000
   24     1.440000     1.000000
   25     1.440000     1.000000
   26     2.040000     1.000000
   27     1.440000     1.000000
   28     1.440000     1.000000

