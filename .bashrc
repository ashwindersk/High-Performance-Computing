# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
module add shared default-environment
module add tools/git-2.18.0
module add languages/python-2.7.6
module add languages/intel-compiler-16-u2
module add openmpi/intel/64/1.6.5  
module add openmpi/gcc/64/1.6.5
module add languages/gcc-7.1.0

module add tools/gnu_builds/tau-2.23.1-openmpi
module add tools/icc_builds/tau-2.23.1-openmpi
