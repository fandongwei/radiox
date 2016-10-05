# @author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
# sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
# Institue: National Astronomical Observatories, Chinese Academy of Sciences
# and also from Chinese Virtual Observatory http://www.china-vo.org/
# please cite http://mnras.oxfordjournals.org/content/451/2/1299

Data files should be end with an empty line, especially the radioxresulttable. Sometimes, the program could not correctly identify the last data row if the row doesn't end with "\n".

1. Install Eigen
	wget http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
	tar xf 3.2.9.tar.bz2
	cd eigen-eigen-dc6cfdf9bcec
	mkdir build
	cd build
	cmake ../
	make install
	make check # check if eigen is in good condition

2. Compile radiox
	cd radiox
	#edit test/Makefile to change configurations
	#or edit test/tester.cpp to write new procedures
	make && make run

