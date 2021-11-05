#ifndef READ_FILE_CC
#define READ_FILE_CC
#include "read_file.h"

#include <iostream>
#include <string.h>

void read_file::readInputFile(char *filename, some_data &dat) {
	FILE *fid;
	int endOfFileFlag;
	char nextLine[MAXLINE];

	int valuesWritten;
	bool fileReadErrorFlag = false;

	fid = std::fopen(filename, "r");
	if (fid == NULL) {
		std::cout << "Unable to open file \"" << filename << "\"" << std::endl;
		fileReadErrorFlag = true;
	} else {

		// read in an Nn
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%u", &(dat.Nn));

		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in Nstep
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%u", &(dat.Nstep));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in the Simulation time
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%lg", &(dat.T));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in the Light Intensity
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%lg", &(dat.LightIntensity));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// Read in the light propagation direction
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		double de1, de2, de3;
		valuesWritten = sscanf(nextLine, "%lg %lg %lg", &(de1), &(de2), &(de3));
		dat.LightDirec[0] = de1;
		dat.LightDirec[1] = de2;
		dat.LightDirec[2] = de3;
		if (valuesWritten != 3) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in w0
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%lg", &(dat.w0));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}
		// read in eta0
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%lg", &(dat.eta0));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in the viscosity
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%lg", &(dat.viscosity));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in the width
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%lg", &(dat.width));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// Read in the bending modulus
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		double BM1, BM2, BM3, BM4;
		dat.BMod.resize(4);
		valuesWritten = sscanf(nextLine, "%lg %lg %lg %lg", &(BM1), &(BM2),
				&(BM3), &(BM4));
		dat.BMod[0] = BM1;
		dat.BMod[1] = BM2;
		dat.BMod[2] = BM3;
		dat.BMod[3] = BM4;

		if (valuesWritten != 4) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in the gravitational force
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%lg", &(dat.rho));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// Read in the force on end of beam
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		double fe1, fe2, fe3;
		valuesWritten = sscanf(nextLine, "%lg %lg %lg", &(fe1), &(fe2), &(fe3));
		dat.fend[0] = fe1;
		dat.fend[1] = fe2;
		dat.fend[2] = fe3;

		if (valuesWritten != 3) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in the Number of threads. 0 means serial
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%u", &(dat.NumThreads));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in the solver tag
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%u", &(dat.SolverTag));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// Read in the Folder Name
		char tempchar[MAXLINE];
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%s", dat.FolderName);
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// Read in the File Name
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%s", dat.FileName);
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}
		/*
		 // Read in the char
		 getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		 valuesWritten = sscanf(nextLine, "%s", dat.someChar);
		 if (valuesWritten != 1)
		 {
		 fileReadErrorFlag = true;
		 goto fileClose;
		 }


		 // read in an int
		 getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		 valuesWritten = sscanf(nextLine, "%u", &(dat.someSize1));
		 if(valuesWritten != 1)
		 {
		 fileReadErrorFlag = true;
		 goto fileClose;
		 }

		 // read in another int int
		 getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		 valuesWritten = sscanf(nextLine, "%u", &(dat.someSize2));
		 if(valuesWritten != 1)
		 {
		 fileReadErrorFlag = true;
		 goto fileClose;
		 }

		 // Read in the doubles
		 getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		 valuesWritten = sscanf(nextLine, "%lg %lg", &(dat.someVal1), &dat.someVal2);
		 if(valuesWritten != 2)
		 {
		 fileReadErrorFlag = true;
		 goto fileClose;
		 }
		 */

		fileClose: {
			fclose(fid);
		}
	}

	if (fileReadErrorFlag) {
		// default parameter values
		std::cout << "Error reading input file, Exiting.\n" << std::endl;
		exit(1);
	} else
		std::cout << "Input file successfully read" << std::endl;

}

void read_file::getNextDataLine(FILE *const filePtr, char *nextLinePtr,
		int const maxSize, int *const endOfFileFlag) {
	*endOfFileFlag = 0;
	do {
		if (fgets(nextLinePtr, maxSize, filePtr) == NULL) {
			*endOfFileFlag = 1;
			break;
		}
		while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t')
				|| (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r')) {
			nextLinePtr = (nextLinePtr + 1);
		}
	} while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
}

#endif
