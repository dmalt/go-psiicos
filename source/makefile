all:
	g++  -fopenmp ActSetViol.cpp -o ActSetViol  -L /opt/OpenBLAS/lib/ -I /opt/OpenBLAS/include/ -lopenblas
	mex BCD_fast.c -L /opt/OpenBLAS/lib -I/opt/OpenBLAS/include/ -lopenblas
	mex columnG_fast.c -L /opt/OpenBLAS/lib -I/opt/OpenBLAS/include/ -lopenblas
