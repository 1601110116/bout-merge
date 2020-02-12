
Drift MHD(3f) and tranport coupling simulation model under BOUT++ framework


Source code
===========

* elm_pb_couple.cxx   -- 3f elm-pb code track the turbulence part
* trans_er_couple.cxx   -- transport code get the profile change


Python scripts
=========
* save-flux.py  -- collect data from file data-elm save flux to grid file
* save-pei.py  -- collect data from file  data-trans and save profile to grid file

Compiling
=========

To compile, run "make" and specify the location of BOUT++
> $ make BOUT_TOP=/path/to/BOUT/
and specify the Source code
> SOURCE = elm_pb_couple.cxx  or trans_er_couple.cxx

Code running
==========

 * To run, execute script
   "debug-elm0.sh" in this directory
  > $ sbatch debug-elm0.sh
 * the other job script can auto-submit
