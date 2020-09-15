%chk=meso-TS_28_1_xp-12.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.130926  -0.423925   0.230614 
  6  -0.402668  -0.125188  -0.661714 
  1  -4.121721  -0.544101   0.651020 
  1   0.324925  -0.126026  -1.477873 
  6   1.104701   1.443859   0.614659 
  6   0.402668   0.125188   0.661714 
  1   0.570306   2.342088   0.897018 
  1  -0.324925   0.126026   1.477873 
  6   1.426074  -0.983719   0.936566 
  1   1.827679  -0.830304   1.941515 
  1   0.940037  -1.960429   0.934664 
  6   2.402045   1.533898   0.142597 
  1   2.860563   2.512891   0.054740 
  6   3.130926   0.423925  -0.230614 
  1   4.121721   0.544101  -0.651020 
  6   2.584558  -0.957939  -0.060346 
  1   2.258123  -1.348285  -1.033742 
  1   3.376997  -1.634728   0.270453 
  6  -1.104701  -1.443859  -0.614659 
  1  -0.570306  -2.342088  -0.897018 
  6  -2.402045  -1.533898  -0.142597 
  1  -2.860563  -2.512891  -0.054740 
  6  -1.426074   0.983719  -0.936566 
  1  -0.940037   1.960429  -0.934664 
  1  -1.827679   0.830304  -1.941515 
  6  -2.584558   0.957939   0.060346 
  1  -3.376997   1.634728  -0.270453 
  1  -2.258122   1.348285   1.033742 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.130926    -0.423925     0.230614     1.700000     1.200000
   -0.402668    -0.125188    -0.661714     1.700000     1.200000
   -4.121721    -0.544101     0.651020     1.200000     1.200000
    0.324925    -0.126026    -1.477873     1.200000     1.200000
    1.104701     1.443859     0.614659     1.700000     1.200000
    0.402668     0.125188     0.661714     1.700000     1.200000
    0.570306     2.342088     0.897018     1.200000     1.200000
   -0.324925     0.126026     1.477873     1.200000     1.200000
    1.426074    -0.983719     0.936566     1.700000     1.200000
    1.827679    -0.830304     1.941515     1.200000     1.200000
    0.940037    -1.960429     0.934664     1.200000     1.200000
    2.402045     1.533898     0.142597     1.700000     1.200000
    2.860563     2.512891     0.054740     1.200000     1.200000
    3.130926     0.423925    -0.230614     1.700000     1.200000
    4.121721     0.544101    -0.651020     1.200000     1.200000
    2.584558    -0.957939    -0.060346     1.700000     1.200000
    2.258123    -1.348285    -1.033742     1.200000     1.200000
    3.376997    -1.634728     0.270453     1.200000     1.200000
   -1.104701    -1.443859    -0.614659     1.700000     1.200000
   -0.570306    -2.342088    -0.897018     1.200000     1.200000
   -2.402045    -1.533898    -0.142597     1.700000     1.200000
   -2.860563    -2.512891    -0.054740     1.200000     1.200000
   -1.426074     0.983719    -0.936566     1.700000     1.200000
   -0.940037     1.960429    -0.934664     1.200000     1.200000
   -1.827679     0.830304    -1.941515     1.200000     1.200000
   -2.584558     0.957939     0.060346     1.700000     1.200000
   -3.376997     1.634728    -0.270453     1.200000     1.200000
   -2.258122     1.348285     1.033742     1.200000     1.200000

--link1--
%chk=meso-TS_28_1_xp-115.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.130926  -0.423925   0.230614 
  6  -0.402668  -0.125188  -0.661714 
  1  -4.121721  -0.544101   0.651020 
  1   0.324925  -0.126026  -1.477873 
  6   1.104701   1.443859   0.614659 
  6   0.402668   0.125188   0.661714 
  1   0.570306   2.342088   0.897018 
  1  -0.324925   0.126026   1.477873 
  6   1.426074  -0.983719   0.936566 
  1   1.827679  -0.830304   1.941515 
  1   0.940037  -1.960429   0.934664 
  6   2.402045   1.533898   0.142597 
  1   2.860563   2.512891   0.054740 
  6   3.130926   0.423925  -0.230614 
  1   4.121721   0.544101  -0.651020 
  6   2.584558  -0.957939  -0.060346 
  1   2.258123  -1.348285  -1.033742 
  1   3.376997  -1.634728   0.270453 
  6  -1.104701  -1.443859  -0.614659 
  1  -0.570306  -2.342088  -0.897018 
  6  -2.402045  -1.533898  -0.142597 
  1  -2.860563  -2.512891  -0.054740 
  6  -1.426074   0.983719  -0.936566 
  1  -0.940037   1.960429  -0.934664 
  1  -1.827679   0.830304  -1.941515 
  6  -2.584558   0.957939   0.060346 
  1  -3.376997   1.634728  -0.270453 
  1  -2.258122   1.348285   1.033742 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.130926    -0.423925     0.230614     1.700000     1.150000
   -0.402668    -0.125188    -0.661714     1.700000     1.150000
   -4.121721    -0.544101     0.651020     1.200000     1.150000
    0.324925    -0.126026    -1.477873     1.200000     1.150000
    1.104701     1.443859     0.614659     1.700000     1.150000
    0.402668     0.125188     0.661714     1.700000     1.150000
    0.570306     2.342088     0.897018     1.200000     1.150000
   -0.324925     0.126026     1.477873     1.200000     1.150000
    1.426074    -0.983719     0.936566     1.700000     1.150000
    1.827679    -0.830304     1.941515     1.200000     1.150000
    0.940037    -1.960429     0.934664     1.200000     1.150000
    2.402045     1.533898     0.142597     1.700000     1.150000
    2.860563     2.512891     0.054740     1.200000     1.150000
    3.130926     0.423925    -0.230614     1.700000     1.150000
    4.121721     0.544101    -0.651020     1.200000     1.150000
    2.584558    -0.957939    -0.060346     1.700000     1.150000
    2.258123    -1.348285    -1.033742     1.200000     1.150000
    3.376997    -1.634728     0.270453     1.200000     1.150000
   -1.104701    -1.443859    -0.614659     1.700000     1.150000
   -0.570306    -2.342088    -0.897018     1.200000     1.150000
   -2.402045    -1.533898    -0.142597     1.700000     1.150000
   -2.860563    -2.512891    -0.054740     1.200000     1.150000
   -1.426074     0.983719    -0.936566     1.700000     1.150000
   -0.940037     1.960429    -0.934664     1.200000     1.150000
   -1.827679     0.830304    -1.941515     1.200000     1.150000
   -2.584558     0.957939     0.060346     1.700000     1.150000
   -3.376997     1.634728    -0.270453     1.200000     1.150000
   -2.258122     1.348285     1.033742     1.200000     1.150000

--link1--
%chk=meso-TS_28_1_xp-110.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.130926  -0.423925   0.230614 
  6  -0.402668  -0.125188  -0.661714 
  1  -4.121721  -0.544101   0.651020 
  1   0.324925  -0.126026  -1.477873 
  6   1.104701   1.443859   0.614659 
  6   0.402668   0.125188   0.661714 
  1   0.570306   2.342088   0.897018 
  1  -0.324925   0.126026   1.477873 
  6   1.426074  -0.983719   0.936566 
  1   1.827679  -0.830304   1.941515 
  1   0.940037  -1.960429   0.934664 
  6   2.402045   1.533898   0.142597 
  1   2.860563   2.512891   0.054740 
  6   3.130926   0.423925  -0.230614 
  1   4.121721   0.544101  -0.651020 
  6   2.584558  -0.957939  -0.060346 
  1   2.258123  -1.348285  -1.033742 
  1   3.376997  -1.634728   0.270453 
  6  -1.104701  -1.443859  -0.614659 
  1  -0.570306  -2.342088  -0.897018 
  6  -2.402045  -1.533898  -0.142597 
  1  -2.860563  -2.512891  -0.054740 
  6  -1.426074   0.983719  -0.936566 
  1  -0.940037   1.960429  -0.934664 
  1  -1.827679   0.830304  -1.941515 
  6  -2.584558   0.957939   0.060346 
  1  -3.376997   1.634728  -0.270453 
  1  -2.258122   1.348285   1.033742 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.130926    -0.423925     0.230614     1.700000     1.100000
   -0.402668    -0.125188    -0.661714     1.700000     1.100000
   -4.121721    -0.544101     0.651020     1.200000     1.100000
    0.324925    -0.126026    -1.477873     1.200000     1.100000
    1.104701     1.443859     0.614659     1.700000     1.100000
    0.402668     0.125188     0.661714     1.700000     1.100000
    0.570306     2.342088     0.897018     1.200000     1.100000
   -0.324925     0.126026     1.477873     1.200000     1.100000
    1.426074    -0.983719     0.936566     1.700000     1.100000
    1.827679    -0.830304     1.941515     1.200000     1.100000
    0.940037    -1.960429     0.934664     1.200000     1.100000
    2.402045     1.533898     0.142597     1.700000     1.100000
    2.860563     2.512891     0.054740     1.200000     1.100000
    3.130926     0.423925    -0.230614     1.700000     1.100000
    4.121721     0.544101    -0.651020     1.200000     1.100000
    2.584558    -0.957939    -0.060346     1.700000     1.100000
    2.258123    -1.348285    -1.033742     1.200000     1.100000
    3.376997    -1.634728     0.270453     1.200000     1.100000
   -1.104701    -1.443859    -0.614659     1.700000     1.100000
   -0.570306    -2.342088    -0.897018     1.200000     1.100000
   -2.402045    -1.533898    -0.142597     1.700000     1.100000
   -2.860563    -2.512891    -0.054740     1.200000     1.100000
   -1.426074     0.983719    -0.936566     1.700000     1.100000
   -0.940037     1.960429    -0.934664     1.200000     1.100000
   -1.827679     0.830304    -1.941515     1.200000     1.100000
   -2.584558     0.957939     0.060346     1.700000     1.100000
   -3.376997     1.634728    -0.270453     1.200000     1.100000
   -2.258122     1.348285     1.033742     1.200000     1.100000

--link1--
%chk=meso-TS_28_1_xp-105.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.130926  -0.423925   0.230614 
  6  -0.402668  -0.125188  -0.661714 
  1  -4.121721  -0.544101   0.651020 
  1   0.324925  -0.126026  -1.477873 
  6   1.104701   1.443859   0.614659 
  6   0.402668   0.125188   0.661714 
  1   0.570306   2.342088   0.897018 
  1  -0.324925   0.126026   1.477873 
  6   1.426074  -0.983719   0.936566 
  1   1.827679  -0.830304   1.941515 
  1   0.940037  -1.960429   0.934664 
  6   2.402045   1.533898   0.142597 
  1   2.860563   2.512891   0.054740 
  6   3.130926   0.423925  -0.230614 
  1   4.121721   0.544101  -0.651020 
  6   2.584558  -0.957939  -0.060346 
  1   2.258123  -1.348285  -1.033742 
  1   3.376997  -1.634728   0.270453 
  6  -1.104701  -1.443859  -0.614659 
  1  -0.570306  -2.342088  -0.897018 
  6  -2.402045  -1.533898  -0.142597 
  1  -2.860563  -2.512891  -0.054740 
  6  -1.426074   0.983719  -0.936566 
  1  -0.940037   1.960429  -0.934664 
  1  -1.827679   0.830304  -1.941515 
  6  -2.584558   0.957939   0.060346 
  1  -3.376997   1.634728  -0.270453 
  1  -2.258122   1.348285   1.033742 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.130926    -0.423925     0.230614     1.700000     1.050000
   -0.402668    -0.125188    -0.661714     1.700000     1.050000
   -4.121721    -0.544101     0.651020     1.200000     1.050000
    0.324925    -0.126026    -1.477873     1.200000     1.050000
    1.104701     1.443859     0.614659     1.700000     1.050000
    0.402668     0.125188     0.661714     1.700000     1.050000
    0.570306     2.342088     0.897018     1.200000     1.050000
   -0.324925     0.126026     1.477873     1.200000     1.050000
    1.426074    -0.983719     0.936566     1.700000     1.050000
    1.827679    -0.830304     1.941515     1.200000     1.050000
    0.940037    -1.960429     0.934664     1.200000     1.050000
    2.402045     1.533898     0.142597     1.700000     1.050000
    2.860563     2.512891     0.054740     1.200000     1.050000
    3.130926     0.423925    -0.230614     1.700000     1.050000
    4.121721     0.544101    -0.651020     1.200000     1.050000
    2.584558    -0.957939    -0.060346     1.700000     1.050000
    2.258123    -1.348285    -1.033742     1.200000     1.050000
    3.376997    -1.634728     0.270453     1.200000     1.050000
   -1.104701    -1.443859    -0.614659     1.700000     1.050000
   -0.570306    -2.342088    -0.897018     1.200000     1.050000
   -2.402045    -1.533898    -0.142597     1.700000     1.050000
   -2.860563    -2.512891    -0.054740     1.200000     1.050000
   -1.426074     0.983719    -0.936566     1.700000     1.050000
   -0.940037     1.960429    -0.934664     1.200000     1.050000
   -1.827679     0.830304    -1.941515     1.200000     1.050000
   -2.584558     0.957939     0.060346     1.700000     1.050000
   -3.376997     1.634728    -0.270453     1.200000     1.050000
   -2.258122     1.348285     1.033742     1.200000     1.050000

--link1--
%chk=meso-TS_28_1_xp-100.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.130926  -0.423925   0.230614 
  6  -0.402668  -0.125188  -0.661714 
  1  -4.121721  -0.544101   0.651020 
  1   0.324925  -0.126026  -1.477873 
  6   1.104701   1.443859   0.614659 
  6   0.402668   0.125188   0.661714 
  1   0.570306   2.342088   0.897018 
  1  -0.324925   0.126026   1.477873 
  6   1.426074  -0.983719   0.936566 
  1   1.827679  -0.830304   1.941515 
  1   0.940037  -1.960429   0.934664 
  6   2.402045   1.533898   0.142597 
  1   2.860563   2.512891   0.054740 
  6   3.130926   0.423925  -0.230614 
  1   4.121721   0.544101  -0.651020 
  6   2.584558  -0.957939  -0.060346 
  1   2.258123  -1.348285  -1.033742 
  1   3.376997  -1.634728   0.270453 
  6  -1.104701  -1.443859  -0.614659 
  1  -0.570306  -2.342088  -0.897018 
  6  -2.402045  -1.533898  -0.142597 
  1  -2.860563  -2.512891  -0.054740 
  6  -1.426074   0.983719  -0.936566 
  1  -0.940037   1.960429  -0.934664 
  1  -1.827679   0.830304  -1.941515 
  6  -2.584558   0.957939   0.060346 
  1  -3.376997   1.634728  -0.270453 
  1  -2.258122   1.348285   1.033742 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.130926    -0.423925     0.230614     1.700000     1.000000
   -0.402668    -0.125188    -0.661714     1.700000     1.000000
   -4.121721    -0.544101     0.651020     1.200000     1.000000
    0.324925    -0.126026    -1.477873     1.200000     1.000000
    1.104701     1.443859     0.614659     1.700000     1.000000
    0.402668     0.125188     0.661714     1.700000     1.000000
    0.570306     2.342088     0.897018     1.200000     1.000000
   -0.324925     0.126026     1.477873     1.200000     1.000000
    1.426074    -0.983719     0.936566     1.700000     1.000000
    1.827679    -0.830304     1.941515     1.200000     1.000000
    0.940037    -1.960429     0.934664     1.200000     1.000000
    2.402045     1.533898     0.142597     1.700000     1.000000
    2.860563     2.512891     0.054740     1.200000     1.000000
    3.130926     0.423925    -0.230614     1.700000     1.000000
    4.121721     0.544101    -0.651020     1.200000     1.000000
    2.584558    -0.957939    -0.060346     1.700000     1.000000
    2.258123    -1.348285    -1.033742     1.200000     1.000000
    3.376997    -1.634728     0.270453     1.200000     1.000000
   -1.104701    -1.443859    -0.614659     1.700000     1.000000
   -0.570306    -2.342088    -0.897018     1.200000     1.000000
   -2.402045    -1.533898    -0.142597     1.700000     1.000000
   -2.860563    -2.512891    -0.054740     1.200000     1.000000
   -1.426074     0.983719    -0.936566     1.700000     1.000000
   -0.940037     1.960429    -0.934664     1.200000     1.000000
   -1.827679     0.830304    -1.941515     1.200000     1.000000
   -2.584558     0.957939     0.060346     1.700000     1.000000
   -3.376997     1.634728    -0.270453     1.200000     1.000000
   -2.258122     1.348285     1.033742     1.200000     1.000000

--link1--
%chk=meso-TS_28_1_xp-0975.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.130926  -0.423925   0.230614 
  6  -0.402668  -0.125188  -0.661714 
  1  -4.121721  -0.544101   0.651020 
  1   0.324925  -0.126026  -1.477873 
  6   1.104701   1.443859   0.614659 
  6   0.402668   0.125188   0.661714 
  1   0.570306   2.342088   0.897018 
  1  -0.324925   0.126026   1.477873 
  6   1.426074  -0.983719   0.936566 
  1   1.827679  -0.830304   1.941515 
  1   0.940037  -1.960429   0.934664 
  6   2.402045   1.533898   0.142597 
  1   2.860563   2.512891   0.054740 
  6   3.130926   0.423925  -0.230614 
  1   4.121721   0.544101  -0.651020 
  6   2.584558  -0.957939  -0.060346 
  1   2.258123  -1.348285  -1.033742 
  1   3.376997  -1.634728   0.270453 
  6  -1.104701  -1.443859  -0.614659 
  1  -0.570306  -2.342088  -0.897018 
  6  -2.402045  -1.533898  -0.142597 
  1  -2.860563  -2.512891  -0.054740 
  6  -1.426074   0.983719  -0.936566 
  1  -0.940037   1.960429  -0.934664 
  1  -1.827679   0.830304  -1.941515 
  6  -2.584558   0.957939   0.060346 
  1  -3.376997   1.634728  -0.270453 
  1  -2.258122   1.348285   1.033742 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.130926    -0.423925     0.230614     1.700000     0.975000
   -0.402668    -0.125188    -0.661714     1.700000     0.975000
   -4.121721    -0.544101     0.651020     1.200000     0.975000
    0.324925    -0.126026    -1.477873     1.200000     0.975000
    1.104701     1.443859     0.614659     1.700000     0.975000
    0.402668     0.125188     0.661714     1.700000     0.975000
    0.570306     2.342088     0.897018     1.200000     0.975000
   -0.324925     0.126026     1.477873     1.200000     0.975000
    1.426074    -0.983719     0.936566     1.700000     0.975000
    1.827679    -0.830304     1.941515     1.200000     0.975000
    0.940037    -1.960429     0.934664     1.200000     0.975000
    2.402045     1.533898     0.142597     1.700000     0.975000
    2.860563     2.512891     0.054740     1.200000     0.975000
    3.130926     0.423925    -0.230614     1.700000     0.975000
    4.121721     0.544101    -0.651020     1.200000     0.975000
    2.584558    -0.957939    -0.060346     1.700000     0.975000
    2.258123    -1.348285    -1.033742     1.200000     0.975000
    3.376997    -1.634728     0.270453     1.200000     0.975000
   -1.104701    -1.443859    -0.614659     1.700000     0.975000
   -0.570306    -2.342088    -0.897018     1.200000     0.975000
   -2.402045    -1.533898    -0.142597     1.700000     0.975000
   -2.860563    -2.512891    -0.054740     1.200000     0.975000
   -1.426074     0.983719    -0.936566     1.700000     0.975000
   -0.940037     1.960429    -0.934664     1.200000     0.975000
   -1.827679     0.830304    -1.941515     1.200000     0.975000
   -2.584558     0.957939     0.060346     1.700000     0.975000
   -3.376997     1.634728    -0.270453     1.200000     0.975000
   -2.258122     1.348285     1.033742     1.200000     0.975000

--link1--
%chk=meso-TS_28_1_xp-095.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.130926  -0.423925   0.230614 
  6  -0.402668  -0.125188  -0.661714 
  1  -4.121721  -0.544101   0.651020 
  1   0.324925  -0.126026  -1.477873 
  6   1.104701   1.443859   0.614659 
  6   0.402668   0.125188   0.661714 
  1   0.570306   2.342088   0.897018 
  1  -0.324925   0.126026   1.477873 
  6   1.426074  -0.983719   0.936566 
  1   1.827679  -0.830304   1.941515 
  1   0.940037  -1.960429   0.934664 
  6   2.402045   1.533898   0.142597 
  1   2.860563   2.512891   0.054740 
  6   3.130926   0.423925  -0.230614 
  1   4.121721   0.544101  -0.651020 
  6   2.584558  -0.957939  -0.060346 
  1   2.258123  -1.348285  -1.033742 
  1   3.376997  -1.634728   0.270453 
  6  -1.104701  -1.443859  -0.614659 
  1  -0.570306  -2.342088  -0.897018 
  6  -2.402045  -1.533898  -0.142597 
  1  -2.860563  -2.512891  -0.054740 
  6  -1.426074   0.983719  -0.936566 
  1  -0.940037   1.960429  -0.934664 
  1  -1.827679   0.830304  -1.941515 
  6  -2.584558   0.957939   0.060346 
  1  -3.376997   1.634728  -0.270453 
  1  -2.258122   1.348285   1.033742 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.130926    -0.423925     0.230614     1.700000     0.950000
   -0.402668    -0.125188    -0.661714     1.700000     0.950000
   -4.121721    -0.544101     0.651020     1.200000     0.950000
    0.324925    -0.126026    -1.477873     1.200000     0.950000
    1.104701     1.443859     0.614659     1.700000     0.950000
    0.402668     0.125188     0.661714     1.700000     0.950000
    0.570306     2.342088     0.897018     1.200000     0.950000
   -0.324925     0.126026     1.477873     1.200000     0.950000
    1.426074    -0.983719     0.936566     1.700000     0.950000
    1.827679    -0.830304     1.941515     1.200000     0.950000
    0.940037    -1.960429     0.934664     1.200000     0.950000
    2.402045     1.533898     0.142597     1.700000     0.950000
    2.860563     2.512891     0.054740     1.200000     0.950000
    3.130926     0.423925    -0.230614     1.700000     0.950000
    4.121721     0.544101    -0.651020     1.200000     0.950000
    2.584558    -0.957939    -0.060346     1.700000     0.950000
    2.258123    -1.348285    -1.033742     1.200000     0.950000
    3.376997    -1.634728     0.270453     1.200000     0.950000
   -1.104701    -1.443859    -0.614659     1.700000     0.950000
   -0.570306    -2.342088    -0.897018     1.200000     0.950000
   -2.402045    -1.533898    -0.142597     1.700000     0.950000
   -2.860563    -2.512891    -0.054740     1.200000     0.950000
   -1.426074     0.983719    -0.936566     1.700000     0.950000
   -0.940037     1.960429    -0.934664     1.200000     0.950000
   -1.827679     0.830304    -1.941515     1.200000     0.950000
   -2.584558     0.957939     0.060346     1.700000     0.950000
   -3.376997     1.634728    -0.270453     1.200000     0.950000
   -2.258122     1.348285     1.033742     1.200000     0.950000

