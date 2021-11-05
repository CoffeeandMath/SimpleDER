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
	//Eigen::Vector3d BMod = { 0.0, 0.0, 0.0 };
	double width = 0.0;
	double height = 0.0;
	double rho = 0.0;
	double length = 1.0;
	double youngsmod = 1.0;
	Eigen::Vector3d fend = { 0.0, 0.0, 0.0 };
	int SolverTag = 0;
	char FolderName[MAXLINE];
	char FileName[MAXLINE];
	double pertsc = 0.0;
	double phiend = 0.0;
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
