#!/usr/bin/env csh
if ($?LD_LIBRARY_PATH>0) then
setenv LD_LIBRARY_PATH '/cm/shared/apps/openmpi/2.1.1/lib:/usr/lib64:/usr/lib64:/usr/lib64:/data2/dschulte/pytom-develop/pytom/pytomc/sh_alignment/SpharmonicKit27/:/data2/dschulte/pytom-develop/pytom/pytomc/sh_alignment/frm/swig/:/data2/dschulte/pytom-develop/pytom/external/lib/:/data2/dschulte/pytom-develop/pytom/pytomc/nufft/:/data2/dschulte/pytom-develop/pytom/pytomc/libs/libtomc/libs':$LD_LIBRARY_PATH
else
setenv LD_LIBRARY_PATH '/cm/shared/apps/openmpi/2.1.1/lib:/usr/lib64:/usr/lib64:/usr/lib64:/data2/dschulte/pytom-develop/pytom/pytomc/sh_alignment/SpharmonicKit27/:/data2/dschulte/pytom-develop/pytom/pytomc/sh_alignment/frm/swig/:/data2/dschulte/pytom-develop/pytom/external/lib/:/data2/dschulte/pytom-develop/pytom/pytomc/nufft/:/data2/dschulte/pytom-develop/pytom/pytomc/libs/libtomc/libs'
endif

if ($?PATH>0) then
setenv PATH '/cm/shared/apps/openmpi/2.1.1/bin:/cm/shared/apps/python/2.7/bin:/data2/dschulte/BachelorThesis/Scripts/:/cm/local/apps/cuda/libs/current/bin:/cm/shared/apps/cuda80/sdk/8.0.61/bin/x86_64/linux/release:/cm/shared/apps/cuda80/toolkit/8.0.61/bin:/cm/shared/apps/utilities/bin:/cm/shared/apps/slurm/17.02.2/sbin:/cm/shared/apps/slurm/17.02.2/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/sbin:/cm/local/apps/environment-modules/3.2.10/bin:/home/dschulte/.local/bin:/home/dschulte/bin':$PATH
else
setenv PATH '/cm/shared/apps/openmpi/2.1.1/bin:/cm/shared/apps/python/2.7/bin:/data2/dschulte/BachelorThesis/Scripts/:/cm/local/apps/cuda/libs/current/bin:/cm/shared/apps/cuda80/sdk/8.0.61/bin/x86_64/linux/release:/cm/shared/apps/cuda80/toolkit/8.0.61/bin:/cm/shared/apps/utilities/bin:/cm/shared/apps/slurm/17.02.2/sbin:/cm/shared/apps/slurm/17.02.2/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/sbin:/cm/local/apps/environment-modules/3.2.10/bin:/home/dschulte/.local/bin:/home/dschulte/bin'
endif

if ($?PYTHONPATH>0) then
setenv PYTHONPATH '/data2/dschulte/pytom-develop/pytom/pytomc:/data2/dschulte/pytom-develop/pytom/pytomc/sh_alignment/frm/swig/:/data2/dschulte/pytom-develop/pytom/external/lib/:/data2/dschulte/pytom-develop/pytom/external/lib/python2.7/site-packages/:/data2/dschulte/pytom-develop/pytom/pytomc/nufft/:/data2/dschulte/pytom-develop':'/data2/dschulte/pytom-develop/pytom/pytomc/swigModules':$PYTHONPATH
else
setenv PYTHONPATH '/data2/dschulte/pytom-develop/pytom/pytomc:/data2/dschulte/pytom-develop/pytom/pytomc/sh_alignment/frm/swig/:/data2/dschulte/pytom-develop/pytom/external/lib/:/data2/dschulte/pytom-develop/pytom/external/lib/python2.7/site-packages/:/data2/dschulte/pytom-develop/pytom/pytomc/nufft/:/data2/dschulte/pytom-develop':'/data2/dschulte/pytom-develop/pytom/pytomc/swigModules'
endif

