# why is the Cavity type SES in Vc calculation using vdw cavity
	this is for Ger calcultion, and ses is better for it.

# are volumes the same in Vc.log and Ger.log?
	yes

# why is the data in VGer.dat (used for EOS fitting) structure-independent

# user-specific spheres -- format differnt in Gcav.com from others
	noaddedsph for vdw cavity; the default is ses cavity

# can I run Gcav calculations before Ger calculations?

# cavitation energy, fixed Vc, mean pressure?

# pressure calculation formula, dG/dV

# are gaussian checkpoint files necessary

# Sept 22
  - multithreading support, add an option in input.jl
  - test system: 96 structures, b3lyp, 6-31g*, int=finegrid
    multi:    124 minutes
    parallel:  
    perl:     407 minutes
  - add guess=read for Ger jobs (excluding the first one that generates the .chk file)
  - add restart feature for Ger jobs. The program will scan through all the existing -Ger.log files and check whether the job is finished (by counting how many "SCF Done" fields are found), and restart the unfinished Ger jobs
