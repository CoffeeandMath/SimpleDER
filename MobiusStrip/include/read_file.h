#ifndef READ_FILE_H
#define READ_FILE_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>
using namespace std;
#define MAXLINE 1024

class some_data {

public:
	some_data() {
	}
	;
	~some_data() {
	}
	;

	int Nn = 0;
	int Nstep = 0;
	Eigen::VectorXd BMod;
	double rho = 0.0;
	double T = 1.0;
	double LightIntensity = 0.0;
	double w0 = 0.0;
	double eta0 = 0.0;
	double viscosity;
	int NumThreads = 0;
	Eigen::Vector3d LightDirec = { 1.0, 0.0, 0.0 };
	Eigen::Vector3d fend = { 0.0, 0.0, 0.0 };
	int SolverTag = 0;
	char FolderName[MAXLINE];
	char FileName[MAXLINE];
	double width = 0.0;
};

class read_file {
public:
	read_file() {
	}
	;
	~read_file() {
	}
	;

	void readInputFile(char *filename, some_data &dat);

private:
	void getNextDataLine(FILE *const filePtr, char *nextLinePtr,
			int const maxSize, int *const endOfFileFlag);
};

#endif
