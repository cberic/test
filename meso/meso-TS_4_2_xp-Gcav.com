%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=108.160262177098 cm^3/mol 

0 1
  6  -3.118154  -0.434970   0.246911 
  6  -0.525244  -0.179089  -0.881108 
  1  -4.082876  -0.549545   0.727275 
  1   0.369421  -0.116466  -1.488305 
  6   1.095749   1.409502   0.684953 
  6   0.525244   0.179089   0.881108 
  1   0.534376   2.315268   0.877539 
  1  -0.369421   0.116466   1.488305 
  6   1.464053  -0.999907   0.987811 
  1   1.902396  -0.948136   1.990906 
  1   0.922907  -1.944294   0.931362 
  6   2.412752   1.518768   0.114990 
  1   2.825512   2.507319  -0.049863 
  6   3.118154   0.434970  -0.246911 
  1   4.082876   0.549545  -0.727275 
  6   2.591960  -0.958481  -0.043904 
  1   2.251216  -1.356631  -1.007828 
  1   3.400803  -1.621376   0.271761 
  6  -1.095749  -1.409502  -0.684954 
  1  -0.534376  -2.315268  -0.877539 
  6  -2.412752  -1.518768  -0.114990 
  1  -2.825512  -2.507319   0.049863 
  6  -1.464053   0.999907  -0.987811 
  1  -0.922907   1.944294  -0.931362 
  1  -1.902396   0.948136  -1.990906 
  6  -2.591960   0.958481   0.043904 
  1  -3.400803   1.621376  -0.271761 
  1  -2.251216   1.356631   1.007828 

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
  6  -3.118154  -0.434970   0.246911 
  6  -0.525244  -0.179089  -0.881108 
  1  -4.082876  -0.549545   0.727275 
  1   0.369421  -0.116466  -1.488305 
  6   1.095749   1.409502   0.684953 
  6   0.525244   0.179089   0.881108 
  1   0.534376   2.315268   0.877539 
  1  -0.369421   0.116466   1.488305 
  6   1.464053  -0.999907   0.987811 
  1   1.902396  -0.948136   1.990906 
  1   0.922907  -1.944294   0.931362 
  6   2.412752   1.518768   0.114990 
  1   2.825512   2.507319  -0.049863 
  6   3.118154   0.434970  -0.246911 
  1   4.082876   0.549545  -0.727275 
  6   2.591960  -0.958481  -0.043904 
  1   2.251216  -1.356631  -1.007828 
  1   3.400803  -1.621376   0.271761 
  6  -1.095749  -1.409502  -0.684954 
  1  -0.534376  -2.315268  -0.877539 
  6  -2.412752  -1.518768  -0.114990 
  1  -2.825512  -2.507319   0.049863 
  6  -1.464053   0.999907  -0.987811 
  1  -0.922907   1.944294  -0.931362 
  1  -1.902396   0.948136  -1.990906 
  6  -2.591960   0.958481   0.043904 
  1  -3.400803   1.621376  -0.271761 
  1  -2.251216   1.356631   1.007828 

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
  6  -3.118154  -0.434970   0.246911 
  6  -0.525244  -0.179089  -0.881108 
  1  -4.082876  -0.549545   0.727275 
  1   0.369421  -0.116466  -1.488305 
  6   1.095749   1.409502   0.684953 
  6   0.525244   0.179089   0.881108 
  1   0.534376   2.315268   0.877539 
  1  -0.369421   0.116466   1.488305 
  6   1.464053  -0.999907   0.987811 
  1   1.902396  -0.948136   1.990906 
  1   0.922907  -1.944294   0.931362 
  6   2.412752   1.518768   0.114990 
  1   2.825512   2.507319  -0.049863 
  6   3.118154   0.434970  -0.246911 
  1   4.082876   0.549545  -0.727275 
  6   2.591960  -0.958481  -0.043904 
  1   2.251216  -1.356631  -1.007828 
  1   3.400803  -1.621376   0.271761 
  6  -1.095749  -1.409502  -0.684954 
  1  -0.534376  -2.315268  -0.877539 
  6  -2.412752  -1.518768  -0.114990 
  1  -2.825512  -2.507319   0.049863 
  6  -1.464053   0.999907  -0.987811 
  1  -0.922907   1.944294  -0.931362 
  1  -1.902396   0.948136  -1.990906 
  6  -2.591960   0.958481   0.043904 
  1  -3.400803   1.621376  -0.271761 
  1  -2.251216   1.356631   1.007828 

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
  6  -3.118154  -0.434970   0.246911 
  6  -0.525244  -0.179089  -0.881108 
  1  -4.082876  -0.549545   0.727275 
  1   0.369421  -0.116466  -1.488305 
  6   1.095749   1.409502   0.684953 
  6   0.525244   0.179089   0.881108 
  1   0.534376   2.315268   0.877539 
  1  -0.369421   0.116466   1.488305 
  6   1.464053  -0.999907   0.987811 
  1   1.902396  -0.948136   1.990906 
  1   0.922907  -1.944294   0.931362 
  6   2.412752   1.518768   0.114990 
  1   2.825512   2.507319  -0.049863 
  6   3.118154   0.434970  -0.246911 
  1   4.082876   0.549545  -0.727275 
  6   2.591960  -0.958481  -0.043904 
  1   2.251216  -1.356631  -1.007828 
  1   3.400803  -1.621376   0.271761 
  6  -1.095749  -1.409502  -0.684954 
  1  -0.534376  -2.315268  -0.877539 
  6  -2.412752  -1.518768  -0.114990 
  1  -2.825512  -2.507319   0.049863 
  6  -1.464053   0.999907  -0.987811 
  1  -0.922907   1.944294  -0.931362 
  1  -1.902396   0.948136  -1.990906 
  6  -2.591960   0.958481   0.043904 
  1  -3.400803   1.621376  -0.271761 
  1  -2.251216   1.356631   1.007828 

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
  6  -3.118154  -0.434970   0.246911 
  6  -0.525244  -0.179089  -0.881108 
  1  -4.082876  -0.549545   0.727275 
  1   0.369421  -0.116466  -1.488305 
  6   1.095749   1.409502   0.684953 
  6   0.525244   0.179089   0.881108 
  1   0.534376   2.315268   0.877539 
  1  -0.369421   0.116466   1.488305 
  6   1.464053  -0.999907   0.987811 
  1   1.902396  -0.948136   1.990906 
  1   0.922907  -1.944294   0.931362 
  6   2.412752   1.518768   0.114990 
  1   2.825512   2.507319  -0.049863 
  6   3.118154   0.434970  -0.246911 
  1   4.082876   0.549545  -0.727275 
  6   2.591960  -0.958481  -0.043904 
  1   2.251216  -1.356631  -1.007828 
  1   3.400803  -1.621376   0.271761 
  6  -1.095749  -1.409502  -0.684954 
  1  -0.534376  -2.315268  -0.877539 
  6  -2.412752  -1.518768  -0.114990 
  1  -2.825512  -2.507319   0.049863 
  6  -1.464053   0.999907  -0.987811 
  1  -0.922907   1.944294  -0.931362 
  1  -1.902396   0.948136  -1.990906 
  6  -2.591960   0.958481   0.043904 
  1  -3.400803   1.621376  -0.271761 
  1  -2.251216   1.356631   1.007828 

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
  6  -3.118154  -0.434970   0.246911 
  6  -0.525244  -0.179089  -0.881108 
  1  -4.082876  -0.549545   0.727275 
  1   0.369421  -0.116466  -1.488305 
  6   1.095749   1.409502   0.684953 
  6   0.525244   0.179089   0.881108 
  1   0.534376   2.315268   0.877539 
  1  -0.369421   0.116466   1.488305 
  6   1.464053  -0.999907   0.987811 
  1   1.902396  -0.948136   1.990906 
  1   0.922907  -1.944294   0.931362 
  6   2.412752   1.518768   0.114990 
  1   2.825512   2.507319  -0.049863 
  6   3.118154   0.434970  -0.246911 
  1   4.082876   0.549545  -0.727275 
  6   2.591960  -0.958481  -0.043904 
  1   2.251216  -1.356631  -1.007828 
  1   3.400803  -1.621376   0.271761 
  6  -1.095749  -1.409502  -0.684954 
  1  -0.534376  -2.315268  -0.877539 
  6  -2.412752  -1.518768  -0.114990 
  1  -2.825512  -2.507319   0.049863 
  6  -1.464053   0.999907  -0.987811 
  1  -0.922907   1.944294  -0.931362 
  1  -1.902396   0.948136  -1.990906 
  6  -2.591960   0.958481   0.043904 
  1  -3.400803   1.621376  -0.271761 
  1  -2.251216   1.356631   1.007828 

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
  6  -3.118154  -0.434970   0.246911 
  6  -0.525244  -0.179089  -0.881108 
  1  -4.082876  -0.549545   0.727275 
  1   0.369421  -0.116466  -1.488305 
  6   1.095749   1.409502   0.684953 
  6   0.525244   0.179089   0.881108 
  1   0.534376   2.315268   0.877539 
  1  -0.369421   0.116466   1.488305 
  6   1.464053  -0.999907   0.987811 
  1   1.902396  -0.948136   1.990906 
  1   0.922907  -1.944294   0.931362 
  6   2.412752   1.518768   0.114990 
  1   2.825512   2.507319  -0.049863 
  6   3.118154   0.434970  -0.246911 
  1   4.082876   0.549545  -0.727275 
  6   2.591960  -0.958481  -0.043904 
  1   2.251216  -1.356631  -1.007828 
  1   3.400803  -1.621376   0.271761 
  6  -1.095749  -1.409502  -0.684954 
  1  -0.534376  -2.315268  -0.877539 
  6  -2.412752  -1.518768  -0.114990 
  1  -2.825512  -2.507319   0.049863 
  6  -1.464053   0.999907  -0.987811 
  1  -0.922907   1.944294  -0.931362 
  1  -1.902396   0.948136  -1.990906 
  6  -2.591960   0.958481   0.043904 
  1  -3.400803   1.621376  -0.271761 
  1  -2.251216   1.356631   1.007828 

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

