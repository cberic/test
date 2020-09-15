%chk=scr.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm with Vmol=108.160262177098 cm^3/mol 

0 1
  6  -3.119734  -0.426930   0.248214 
  6  -0.452630  -0.134377  -0.766522 
  1  -4.086028  -0.547352   0.722629 
  1   0.332703  -0.111776  -1.521768 
  6   1.115166   1.430657   0.672626 
  6   0.452630   0.134377   0.766522 
  1   0.561998   2.331699   0.905522 
  1  -0.332703   0.111776   1.521768 
  6   1.457223  -0.998723   0.982555 
  1   1.882127  -0.904241   1.986471 
  1   0.944658  -1.961243   0.946268 
  6   2.397689   1.528359   0.123394 
  1   2.825979   2.512496  -0.031045 
  6   3.119734   0.426930  -0.248214 
  1   4.086028   0.547352  -0.722629 
  6   2.592798  -0.959299  -0.042549 
  1   2.246207  -1.358828  -1.004791 
  1   3.401603  -1.624548   0.271052 
  6  -1.115166  -1.430657  -0.672626 
  1  -0.561998  -2.331699  -0.905522 
  6  -2.397689  -1.528359  -0.123394 
  1  -2.825979  -2.512496   0.031045 
  6  -1.457223   0.998723  -0.982556 
  1  -0.944658   1.961243  -0.946268 
  1  -1.882127   0.904241  -1.986471 
  6  -2.592798   0.959299   0.042549 
  1  -3.401603   1.624548  -0.271052 
  1  -2.246207   1.358828   1.004791 

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
  6  -3.119734  -0.426930   0.248214 
  6  -0.452630  -0.134377  -0.766522 
  1  -4.086028  -0.547352   0.722629 
  1   0.332703  -0.111776  -1.521768 
  6   1.115166   1.430657   0.672626 
  6   0.452630   0.134377   0.766522 
  1   0.561998   2.331699   0.905522 
  1  -0.332703   0.111776   1.521768 
  6   1.457223  -0.998723   0.982555 
  1   1.882127  -0.904241   1.986471 
  1   0.944658  -1.961243   0.946268 
  6   2.397689   1.528359   0.123394 
  1   2.825979   2.512496  -0.031045 
  6   3.119734   0.426930  -0.248214 
  1   4.086028   0.547352  -0.722629 
  6   2.592798  -0.959299  -0.042549 
  1   2.246207  -1.358828  -1.004791 
  1   3.401603  -1.624548   0.271052 
  6  -1.115166  -1.430657  -0.672626 
  1  -0.561998  -2.331699  -0.905522 
  6  -2.397689  -1.528359  -0.123394 
  1  -2.825979  -2.512496   0.031045 
  6  -1.457223   0.998723  -0.982556 
  1  -0.944658   1.961243  -0.946268 
  1  -1.882127   0.904241  -1.986471 
  6  -2.592798   0.959299   0.042549 
  1  -3.401603   1.624548  -0.271052 
  1  -2.246207   1.358828   1.004791 

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
  6  -3.119734  -0.426930   0.248214 
  6  -0.452630  -0.134377  -0.766522 
  1  -4.086028  -0.547352   0.722629 
  1   0.332703  -0.111776  -1.521768 
  6   1.115166   1.430657   0.672626 
  6   0.452630   0.134377   0.766522 
  1   0.561998   2.331699   0.905522 
  1  -0.332703   0.111776   1.521768 
  6   1.457223  -0.998723   0.982555 
  1   1.882127  -0.904241   1.986471 
  1   0.944658  -1.961243   0.946268 
  6   2.397689   1.528359   0.123394 
  1   2.825979   2.512496  -0.031045 
  6   3.119734   0.426930  -0.248214 
  1   4.086028   0.547352  -0.722629 
  6   2.592798  -0.959299  -0.042549 
  1   2.246207  -1.358828  -1.004791 
  1   3.401603  -1.624548   0.271052 
  6  -1.115166  -1.430657  -0.672626 
  1  -0.561998  -2.331699  -0.905522 
  6  -2.397689  -1.528359  -0.123394 
  1  -2.825979  -2.512496   0.031045 
  6  -1.457223   0.998723  -0.982556 
  1  -0.944658   1.961243  -0.946268 
  1  -1.882127   0.904241  -1.986471 
  6  -2.592798   0.959299   0.042549 
  1  -3.401603   1.624548  -0.271052 
  1  -2.246207   1.358828   1.004791 

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
  6  -3.119734  -0.426930   0.248214 
  6  -0.452630  -0.134377  -0.766522 
  1  -4.086028  -0.547352   0.722629 
  1   0.332703  -0.111776  -1.521768 
  6   1.115166   1.430657   0.672626 
  6   0.452630   0.134377   0.766522 
  1   0.561998   2.331699   0.905522 
  1  -0.332703   0.111776   1.521768 
  6   1.457223  -0.998723   0.982555 
  1   1.882127  -0.904241   1.986471 
  1   0.944658  -1.961243   0.946268 
  6   2.397689   1.528359   0.123394 
  1   2.825979   2.512496  -0.031045 
  6   3.119734   0.426930  -0.248214 
  1   4.086028   0.547352  -0.722629 
  6   2.592798  -0.959299  -0.042549 
  1   2.246207  -1.358828  -1.004791 
  1   3.401603  -1.624548   0.271052 
  6  -1.115166  -1.430657  -0.672626 
  1  -0.561998  -2.331699  -0.905522 
  6  -2.397689  -1.528359  -0.123394 
  1  -2.825979  -2.512496   0.031045 
  6  -1.457223   0.998723  -0.982556 
  1  -0.944658   1.961243  -0.946268 
  1  -1.882127   0.904241  -1.986471 
  6  -2.592798   0.959299   0.042549 
  1  -3.401603   1.624548  -0.271052 
  1  -2.246207   1.358828   1.004791 

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
  6  -3.119734  -0.426930   0.248214 
  6  -0.452630  -0.134377  -0.766522 
  1  -4.086028  -0.547352   0.722629 
  1   0.332703  -0.111776  -1.521768 
  6   1.115166   1.430657   0.672626 
  6   0.452630   0.134377   0.766522 
  1   0.561998   2.331699   0.905522 
  1  -0.332703   0.111776   1.521768 
  6   1.457223  -0.998723   0.982555 
  1   1.882127  -0.904241   1.986471 
  1   0.944658  -1.961243   0.946268 
  6   2.397689   1.528359   0.123394 
  1   2.825979   2.512496  -0.031045 
  6   3.119734   0.426930  -0.248214 
  1   4.086028   0.547352  -0.722629 
  6   2.592798  -0.959299  -0.042549 
  1   2.246207  -1.358828  -1.004791 
  1   3.401603  -1.624548   0.271052 
  6  -1.115166  -1.430657  -0.672626 
  1  -0.561998  -2.331699  -0.905522 
  6  -2.397689  -1.528359  -0.123394 
  1  -2.825979  -2.512496   0.031045 
  6  -1.457223   0.998723  -0.982556 
  1  -0.944658   1.961243  -0.946268 
  1  -1.882127   0.904241  -1.986471 
  6  -2.592798   0.959299   0.042549 
  1  -3.401603   1.624548  -0.271052 
  1  -2.246207   1.358828   1.004791 

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
  6  -3.119734  -0.426930   0.248214 
  6  -0.452630  -0.134377  -0.766522 
  1  -4.086028  -0.547352   0.722629 
  1   0.332703  -0.111776  -1.521768 
  6   1.115166   1.430657   0.672626 
  6   0.452630   0.134377   0.766522 
  1   0.561998   2.331699   0.905522 
  1  -0.332703   0.111776   1.521768 
  6   1.457223  -0.998723   0.982555 
  1   1.882127  -0.904241   1.986471 
  1   0.944658  -1.961243   0.946268 
  6   2.397689   1.528359   0.123394 
  1   2.825979   2.512496  -0.031045 
  6   3.119734   0.426930  -0.248214 
  1   4.086028   0.547352  -0.722629 
  6   2.592798  -0.959299  -0.042549 
  1   2.246207  -1.358828  -1.004791 
  1   3.401603  -1.624548   0.271052 
  6  -1.115166  -1.430657  -0.672626 
  1  -0.561998  -2.331699  -0.905522 
  6  -2.397689  -1.528359  -0.123394 
  1  -2.825979  -2.512496   0.031045 
  6  -1.457223   0.998723  -0.982556 
  1  -0.944658   1.961243  -0.946268 
  1  -1.882127   0.904241  -1.986471 
  6  -2.592798   0.959299   0.042549 
  1  -3.401603   1.624548  -0.271052 
  1  -2.246207   1.358828   1.004791 

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
  6  -3.119734  -0.426930   0.248214 
  6  -0.452630  -0.134377  -0.766522 
  1  -4.086028  -0.547352   0.722629 
  1   0.332703  -0.111776  -1.521768 
  6   1.115166   1.430657   0.672626 
  6   0.452630   0.134377   0.766522 
  1   0.561998   2.331699   0.905522 
  1  -0.332703   0.111776   1.521768 
  6   1.457223  -0.998723   0.982555 
  1   1.882127  -0.904241   1.986471 
  1   0.944658  -1.961243   0.946268 
  6   2.397689   1.528359   0.123394 
  1   2.825979   2.512496  -0.031045 
  6   3.119734   0.426930  -0.248214 
  1   4.086028   0.547352  -0.722629 
  6   2.592798  -0.959299  -0.042549 
  1   2.246207  -1.358828  -1.004791 
  1   3.401603  -1.624548   0.271052 
  6  -1.115166  -1.430657  -0.672626 
  1  -0.561998  -2.331699  -0.905522 
  6  -2.397689  -1.528359  -0.123394 
  1  -2.825979  -2.512496   0.031045 
  6  -1.457223   0.998723  -0.982556 
  1  -0.944658   1.961243  -0.946268 
  1  -1.882127   0.904241  -1.986471 
  6  -2.592798   0.959299   0.042549 
  1  -3.401603   1.624548  -0.271052 
  1  -2.246207   1.358828   1.004791 

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

