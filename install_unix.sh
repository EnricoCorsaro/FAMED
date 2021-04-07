#!/bin/bash
# First check whether the installation of DIAMONDS is required. Give it priority with respect to
# the other software packages.
# Usage: ./install_unix.sh [--diamonds | -d] [--background | -b] [--peakbagging | -p] [--asymptotic | -a] [--parallel | -g]
flag1=0
flag2=0
flag3=0
flag4=0
flag5=0

if ! [ -x "$(command -v git)" ]; then
	echo "Error: git is not installed. Aborting..." >&2
	exit 1
fi

if ! [ -x "$(command -v cmake)" ]; then
	echo "Cmake is not installed. Trying to install it using apt-get..." >&2

	if ! [ -x "$(command -v apt-get)" ]; then
		echo "Error: apt-get is not installed. Aborting..." >&2
		exit 1
	else
		sudo apt-get install cmake
	fi
fi

while [ ! $# -eq 0 ]
do
	case "$1" in
		--diamonds | -d)
			flag1=1
			;;
		--background | -b)
			flag2=1
			;;
		--peakbagging | -p)
			flag3=1
			;;
		--asymptotic | -a)
			flag4=1
			;;
		--parallel | -g)
			flag5=1
			;;
		*) 
			echo "Flag $1 not recognized. Only flags -d, -b, -p, -a, -g -y are allowed. Aborting..."
			exit 1
			;;
	esac
	shift
done

if [ $flag1 -eq 1 ]; then
	echo "-----------------------------------"
	echo " Cloning and installing DIAMONDS..."
	echo "-----------------------------------"
	git clone https://github.com/EnricoCorsaro/DIAMONDS.git
	mkdir DIAMONDS/build
	cd DIAMONDS/build/
	cmake ..
	make -j 4
	echo "-----------------------------------"
	echo " Compiling and running test demo..."
	echo "-----------------------------------"
	cd ../../
	export LD_LIBRARY_PATH=${PWD}/DIAMONDS/build
	cd DIAMONDS/demos/
	g++ -o demoSingle2DGaussian demoSingle2DGaussian.cpp -L../build/ -I../include/ -ldiamonds -std=c++11
	./demoSingle2DGaussian
	cd ../../
	echo " "
	echo " "	
fi

if [ $flag2 -eq 1 ]; then
	echo "-------------------------------------------"
	echo " Cloning and installing Background code..."
	echo "-------------------------------------------"
	git clone https://github.com/EnricoCorsaro/Background.git
	mkdir Background/build
	cd Background/build/
	cmake ..
	make -j 4
	cd ../../
	echo " "
	echo " "
	echo "-------------------------------------------------------------------------"
	echo " Preparing localPath.txt file with local path of the Background folder..."
	echo "-------------------------------------------------------------------------"
	echo "${PWD}/Background/" > ${PWD}/Background/build/localPath.txt
	echo " "
	echo " "
	echo "--------------------------------------"
	echo " Setting up tutorial for Background..."
	echo "--------------------------------------"
	mkdir Background/data
	mkdir Background/results
	mv Background/tutorials/KIC012008916/KIC012008916.txt Background/data/
	mv Background/tutorials/KIC012008916 Background/results/
	mkdir Background/results/KIC012008916/00
	rm Background/results/KIC012008916/localPath.txt
	echo " "
	echo " "
fi

if [ $flag3 -eq 1 ]; then
	echo "-------------------------------------------"
	echo " Cloning and installing PeakBagging code..."
	echo "-------------------------------------------"
	git clone https://github.com/EnricoCorsaro/PeakBagging.git
	mkdir PeakBagging/build
	cd PeakBagging/build/
	cmake ..
	make -j 4
	cd ../../
	echo " "
	echo " "
	echo "----------------------------------------------------------------------------"
	echo " Changing localPath.txt content with local path of the PeakBagging folder..."
	echo "----------------------------------------------------------------------------"
	echo "${PWD}/PeakBagging/" > ${PWD}/PeakBagging/build/localPath.txt
	echo " "
	echo " "
	echo "---------------------------------------"
	echo " Setting up tutorial for PeakBagging..."
	echo "---------------------------------------"
	mkdir PeakBagging/data
	mkdir PeakBagging/results
	mv PeakBagging/tutorials/KIC012008916/KIC012008916.txt PeakBagging/data/
	mv PeakBagging/tutorials/KIC012008916 PeakBagging/results/
	mkdir PeakBagging/results/KIC012008916/pb/0
	mkdir PeakBagging/results/KIC012008916/pb/0A
	rm PeakBagging/results/KIC012008916/localPath.txt
	echo " "
	echo " "
fi

if [ $flag4 -eq 1 ]; then
	echo "-------------------------------------------"
	echo " Cloning and installing Asymptotic code... "
	echo "-------------------------------------------"
	git clone https://github.com/EnricoCorsaro/Asymptotic.git
	mkdir Asymptotic/build
	cd Asymptotic/build/
	cmake ..
	make -j 4
	cd ../../
	echo " "
	echo " "
	echo "----------------------------------------------------------------------------"
	echo " Changing localPath.txt content with local path of the Asymptotic folder..."
	echo "----------------------------------------------------------------------------"
	echo "${PWD}/PeakBagging/" > ${PWD}/Asymptotic/build/localPath.txt
fi

if [ $flag5 -eq 1 ]; then
	echo " "
	echo "-----------------------------------------------------------"
	echo " Downloading and installing GNU parallel (v. 2019-12-22)..."
	echo "-----------------------------------------------------------"
	if ! [ -x "$(command -v wget)" ]; then
		echo "Error: wget is not installed. Trying with curl..." >&2
		
		if ! [ -x "$(command -v curl)" ]; then
			echo "Error: curl is not installed. Aborting..." >&2
			exit 1
		else
			curl -O http://ftp.gnu.org/gnu/parallel/parallel-20191222.tar.bz2
		fi
	else 
		wget http://ftp.gnu.org/gnu/parallel/parallel-20191222.tar.bz2
	fi
	tar -xvjf parallel-20191222.tar.bz2
	cd parallel-20191222
	./configure --prefix=$HOME && make && make install
	cd ..
	rm -r parallel-20191222
	rm parallel-20191222.tar.bz2
fi

echo "----------------------------------------------"
echo " Cloning and installing the FAMED pipeline... "
echo "----------------------------------------------"
git clone https://github.com/EnricoCorsaro/FAMED.git
echo " "
echo " "
echo "----------------------------------------------------------------------"
echo " Changing configuring parameters of the pipeline with local folders..."
echo "----------------------------------------------------------------------"
sed -i.old "s^YOUR_LOCAL_ROOT_PATH_HERE^${PWD}^g" ${PWD}/FAMED/idl/famed_configuring_parameters.txt
sed -i.old "s^YOUR_LOCAL_ROOT_PATH_HERE^${PWD}^g" ${PWD}/FAMED/python/famed/famed_configuring_parameters.txt
sed -i.old "s^YOUR_LOCAL_ROOT_PATH_HERE^${PWD}^g" ${PWD}/FAMED/python/famed/famed_config.yml

echo "----------------------------------------------------------------------"
echo " Adding FAMED to your local Python library path..."
echo "----------------------------------------------------------------------"
MY_LOCAL_DIR=${PWD}
cd ~
echo "" >> .bashrc
echo "export PYTHONPATH="${MY_LOCAL_DIR}"/FAMED/python/" >> .bashrc