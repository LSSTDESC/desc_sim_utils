#!/bin/bash
source scl_source enable devtoolset-6
source loadLSST.bash
setup lsst_sims
pip install nose
pip install coveralls
git clone https://github.com/lsst/obs_lsst.git
cd obs_lsst
setup -r . -j
scons
cd ..
eups declare desc_sim_utils -r ${TRAVIS_BUILD_DIR} -t current
setup desc_sim_utils
cd ${TRAVIS_BUILD_DIR}
scons
nosetests -s --with-coverage --cover-package=desc.sim_utils
