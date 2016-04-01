/*
 * search.cpp
 *
 *  Created on: Oct 14, 2010
 *      Author: xiaoliwang
 */

#include "Time.h"
#include "Gram.h"
#include "SeqDB.h"
#include "Query.h"

#include <vector>
#include <set>
#include <string>
#include <fstream>

using namespace std;

int gl_qmax = 5;
int gl_qmin = 2;
vector<string> dataset;
vector<string> dataset_outline;
unsigned gl_maxLen = 0;
unsigned gl_minLen = 65535;
bool output = false;

void usage()
{
	fprintf(stdout,"**********************************************\n");
	fprintf(stdout,"SSEARCH 1.0 Copyright 2010, NUS, Singapore\n");
	fprintf(stdout,"Build the index for the sequence search\n");
	fprintf(stdout,"----------------------------------------------\n\n");
	fprintf(stdout,"Usage: sindex [OPTION] ${INPUTDB}\n");
	fprintf(stdout,"-m --integer\n   upper gram length (default by 5).\n");
	fprintf(stdout,"-n --integer\n   lower gram length (default by 2).\n");
	fprintf(stdout,"-o\n   output the index levels separately or not.\n");
	fprintf(stdout,"\n----------------------------------------------\n");
	fprintf(stdout,"Author : Xiaoli Wang, SOC, NUS, Singapore.\n");
	fprintf(stdout,"GMail  : xiaolistue@gmail.com");
	fprintf(stdout,"\n**********************************************\n");
}

void parseOptions(int argc, const char* argv[])
{
	if (argc > 2){
		for (int i=1;argc-1>i;i++){
			if (strncmp(argv[i], "-m", 2) == 0) {
				if (1 != sscanf(argv[++i], "%d", &gl_qmax)) {
					fprintf(stdout, "Error occured while parsing qmax value.\n");
					exit(-1);
				}
				if (gl_qmax < 0) {
					fprintf(stdout, "Qmax value error: %d.\n", gl_qmax);
					exit(-1);
				}
			}

			if (strncmp(argv[i], "-n", 2) == 0) {
				if (1 != sscanf(argv[++i], "%d", &gl_qmin)) {
					fprintf(stdout, "Error occured while parsing qmin value.\n");
					exit(-1);
				}
				if (gl_qmin < 0) {
					fprintf(stdout, "Qmin value error: %d.\n", gl_qmin);
					exit(-1);
				}
			}

			if (gl_qmax < gl_qmin) {
				fprintf(stdout, "Error qmax and qmin vlues: qmax(%d) < qmin(%d). Swap them for queries!\n", gl_qmax, gl_qmin);
				int _tmp = gl_qmin;
				gl_qmin = gl_qmax;
				gl_qmax = _tmp;
			}

			if (gl_qmin >(gl_qmax + 1) / 2) {
				gl_qmin = (gl_qmax / 2) < 3 ? (gl_qmax / 2) : 3;
				fprintf(stdout, "The qmax and qmin values are reset to: qmax(%d) > qmin(%d) for algorithms!\n", gl_qmax, gl_qmin);
			}

			if (strncmp(argv[i],"-o",2) ==0){
				output = true;
			}
		}
	}
}

void read(const char *fileName, vector<string> &data)
{
	ifstream _fs;
	_fs.open(fileName, ios::in);
	if (_fs.fail()) {
		fprintf(stderr,"Error: Failed to open the data file - %s\n",fileName);
		fprintf(stderr,"Please check the file name and restart this program.\n");
		fprintf(stderr,"The program will be terminated!\n");
		exit(-1);
	}

	string line;
	while (!_fs.eof()) {
		getline(_fs, line);

		if (line.empty())
			continue;

		for (unsigned i = 0; line.size() > i; i++)
			if (line[i] >= 'A' && line[i] <= 'Z')
				line[i] += 32;

		if (line[line.size() - 1] == '\r')
			line.erase(line.size() - 1);

		if (gl_maxLen < line.size())
			gl_maxLen = line.size();

		if (gl_minLen > line.size())
			gl_minLen = line.size();

		data.push_back(line);
	}

	_fs.close();
}

int main(int argc, const char* argv[])
{
	if (argc < 2){
		usage();
		exit(-1);
	}

	// parameters
	parseOptions(argc, argv);

	// read data file
	read(argv[argc-1], dataset);
	if (dataset.size() == 0) {
		fprintf(stderr,"Error: Failed to open the dataset file - %s\n",argv[argc-1]);
		fprintf(stderr,"Please check the file name and restart this program.\n");
		fprintf(stderr,"The program will be terminated!\n");
		exit(-1);
	}

    /*
	// Set parameters
	if (gl_qmax > gl_minLen) {
		gl_qmax = gl_minLen;
		if (gl_qmin > (gl_qmax + 1) / 2) {
			gl_qmin = gl_qmax / 2;
			fprintf(stdout, "The qmax and qmin values are reset to: qmax(%d) > qmin(%d) for algorithms!\n", gl_qmax, gl_qmin);
		}
	}
	*/

	class TSINGHUA_CLIPSE_UTIL::TimeRecorder _time;

	CGram gramgenupper(gl_qmax, false);
	CGram gramgenlower(gl_qmin, false);
	CSeqDB<> db(gl_maxLen, 1, &dataset, &gramgenupper, &gramgenlower);
	db.buildindex();

	_time.check();
	fprintf(stdout, "Total construction time: %-10.6f(sec)\n", _time.diffTime(0, 1));

	string index = argv[argc-1];
	index += ".ix";
	db.save(index.c_str(), output);

	//print the index: upper and lower
	FILE* fr;
	fr = fopen("indexfile", "w");
	if (fr == NULL)
        fprintf(stdout, "index file open failed\n");
    else
        db.print(fr);
    fclose(fr);
}
