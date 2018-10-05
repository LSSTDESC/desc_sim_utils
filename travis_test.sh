#!/bin/bash
source scl_source enable devtoolset-6
source loadLSST.bash
setup lsst_sims
pip install nose
pip install coveralls
git clone https://github.com/lsst/obs_lsstCam.git
cd obs_lsstCam
git checkout u/jchiang/imsim_crosstalk
setup -r . -j
scons
cd ..
git clone https://github.com/LSSTDESC/imSim.git
cd imSim
setup -r . -j
scons
cd ..
eups declare desc_sim_utils -r ${TRAVIS_BUILD_DIR} -t current
setup desc_sim_utils
cd ${TRAVIS_BUILD_DIR}
scons
nosetests -s --with-coverage --cover-package=desc.sim_utils
