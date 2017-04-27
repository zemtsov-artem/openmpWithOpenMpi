#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <ctime>
#include <vector>
#include <omp.h>
#include "mpi.h"
#define pathToTheFile "WP.txt"
#define NUMBER_OF_OPENMP_THR 4
#define NUMBER_OF_MPI_THR 4

long hashFunc(const std::string _inputString) {
	long sum = 0;
	for (int i = 0; i < _inputString.size(); ++i) {
		sum += (int)_inputString.at(i) * pow(101, _inputString.size() - i - 1);
	}
	return sum;
}

int RabinKarp(const std::string mainStr, const std::string sub, const int openMpThreadNum,const int openMpiThreadNum) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int diap = mainStr.size() / openMpiThreadNum;
	int rest = mainStr.size()%openMpiThreadNum;
	int start = diap * rank;
	(rank == openMpiThreadNum - 1) : diap = mainStr.size() - start + 1 : diap++;
	long desiredSubStrHash = hashFunc(sub);

		#pragma omp parallel num_threads(openMpThreadNum)  reduction(+:num) firstprivate(desiredSubStrHash)
		{
		#pragma omp for //nowait
			for (int diapIter = start; (diapIter < start + diap); diapIter++)
			{
				long _hs;
				if (diapIter == 0) {
					_hs = hashFunc(mainStr.substr(0, sub.size()));
				}
				if ((_hs == desiredSubStrHash) && (mainStr.substr(diapIter, sub.size()) == sub))
				{
					num++;
				}
				_hs = hashFunc(mainStr.substr(diapIter + 1, sub.size()));
			}
		}
	return num;
}

int main(int argc, char const *argv[]){
	std::string str,strTotal,ourSubString;
	std::ifstream in;
	int proc_num,rank,result = 0,subRes = 0;
	in.open(pathToTheFile);
	while ( in ) {
		getline(in,str);
 	  	strTotal += str;
	}
	in.close();

	printf("Enter your substring\n");
	std::cin >> ourSubString;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	omp_set_num_threads(NUMBER_OF_MPI_THR);

	clock_t begin_time = clock();

	result = RabinKarp(stdmainStr, stdsub, NUMBER_OF_OPENMP_THR,NUMBER_OF_MPI_THR);

	MPI_Reduce(&subRes, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	clock_t stop_time = clock();
	float parr_time = float(stop_time - begin_time) / CLOCKS_PER_SEC;
	MPI_Finalize();
	printf("Parallel time on %d threads \n %d processes: %d \n Result = %d \n",NUMBER_OF_OPENMP_THR,NUMBER_OF_MPI_THR,parr_time,Result)
	return 0;
}