# MPIによる並列計算

2次元ポアソン方程式をMPI並列で高速化する

1. シリアル版
diffusion.cpp

2. MPI並列版
diffusion_mpi.cpp

## How to build
1. シリアル版
```
icpc diffusion.cpp
```
or
``` 
g++ diffusion.cpp
```
2. MPI版
```
pmiicpc diffusion_mpi.cpp
```
## How to run
1. シリアル版
```
./a.out
```
2. MPI並列版
```
mpirun -np 並列数 ./a.out
```
現在はMAX8並列までしか対応していない
