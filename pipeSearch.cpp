/*
 * search.cpp
 *
 *  Created on: Oct 14, 2010
 *      Author: xiaoliwang
 */
#include "Time.h"
#include "Gram.h"
#include "CountFilter.h"
#include "SeqDB.h"
#include "Query.h"

#include<iostream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <pthread.h>
//add the unistd.h for sleep in 173
#include<unistd.h>

using namespace std;

int gl_qmax = 5;
int gl_qmin = 2;
int gl_topk = 1;
int gl_topk_f = 1;
bool output = false;
vector<string> dataset;
vector<string> queryset;
unsigned gl_maxLen = 0;

//double queryPrep_time = 0.0;
double query_time = 0.0;

void usage()
{
	fprintf(stdout,"**********************************************\n");
	fprintf(stdout,"SSEARCH 1.0 Copyright 2010, NUS, Singapore\n");
	fprintf(stdout,"Sequence search algorithms based on CA\n");
	fprintf(stdout,"----------------------------------------------\n\n");
	fprintf(stdout,"Usage: psearch [OPTION] ${INPUTDB} ${QUERY}\n");
	fprintf(stdout,"-k --integer\n   k value (default by 1).\n");
	fprintf(stdout,"-f --integer\n   k value for the f-queue (default by 1).\n");
	fprintf(stdout,"-o\n   ouput the result set of sequence ids.\n");
	fprintf(stdout,"\n----------------------------------------------\n");
	fprintf(stdout,"Author : Xiaoli Wang, SOC, NUS, Singapore.\n");
	fprintf(stdout,"GMail  : xiaolistue@gmail.com");
	fprintf(stdout,"\n**********************************************\n");
}

//compare the strings between psearch and dataset, and get some possible arguments
void parseOptions(int argc, const char* argv[])
{
	if (argc > 3){
		for (int i=1;argc-2>i;i++){
            //compare first size chars of string1 and string2, if true return 0
			if (strncmp(argv[i],"-k",2) == 0){
				if (1 != sscanf(argv[++i],"%d",&gl_topk)){
					fprintf(stdout,"Error occured while parsing the value of k.\n");
					exit(-1);
				}
				if (gl_topk<=0){
					fprintf(stdout,"The top k value error : %d.\n",gl_topk);
					exit(-1);
				}
			}

			if (strncmp(argv[i],"-f",2) ==0){
				if (1 != sscanf(argv[++i],"%d",&gl_topk_f)){
					fprintf(stdout,"Error occured while parsing the value of k for f-queue.\n");
					exit(-1);
				}
				if (gl_topk_f<=0){
					fprintf(stdout,"The queue k value error : %d.\n",gl_topk_f);
					exit(-1);
				}
			}

			if (strncmp(argv[i],"-o",2) ==0){
				output = true;
			}
		}
	}
}

//string : file(strings) -> data
void read(const char *fileName, vector<string> &data)
{
	ifstream _fs;
	_fs.open(fileName, ios::in);
	//failed to read the file
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

		data.push_back(line);
	}

	_fs.close();
}

static CSeqDB<> *pThis;
volatile int threadnum = 0;
pthread_mutex_t running_mutex = PTHREAD_MUTEX_INITIALIZER;
#define PTHEAD_NUM 2

// The fuction of the thread
void *thr_fn_knn(void *arg)
{
    //pthread_self() : use it for the thread id
	int iret = pthread_detach(pthread_self());
	if(iret != 0)
	{
		fprintf(stderr, "Can't detach at my thread!\n");
		pthread_mutex_lock(&running_mutex);
		threadnum--;
		pthread_mutex_unlock(&running_mutex);
		pthread_exit(NULL);   // Exit the thread
		return ((void *)0);
	}

	// Accumulate frequency
	long ged = (long)arg;
	pThis->accumulateFrequency(ged);

	pthread_mutex_lock(&running_mutex);
	threadnum--;
	pthread_mutex_unlock(&running_mutex);
	pthread_exit(NULL);
	return ((void *)0);
}

int pipe_knn(CQuery* query) {
	pThis->theQuery = query;
	pthread_t threads[PTHEAD_NUM];
	int error = 0;
	for (int i = 0; i < PTHEAD_NUM; i++)
	{
		pthread_mutex_lock(&running_mutex);
		threadnum++;
		pthread_mutex_unlock(&running_mutex);
		error = pthread_create(&threads[i], NULL, thr_fn_knn, (void *)i);
		if (0 != error)
		{
			fprintf(stderr, "Cannot create thread:  %d!\n", error);
			return 0;
		}
	}

	// Check if all threads are closed.
	while (true)
	{
		if (threadnum <= 0)
		{
			// Try post precessing here with more approximate bounds
			pThis->knn_postprocess();
			// Then end of whole processing here
			return 1;
		}
		sleep(0.0001);
	}
	return 0;
}


int main(int argc, const char* argv[])
{
    //./psearch dataset query
    //the length of the array argv is at least 3
	if (argc < 3){
		usage();
		exit(-1);
	}

	// parameters
	parseOptions(argc, argv);

	// read data file
	read(argv[argc-2], dataset);
	if (dataset.size() == 0) {
		fprintf(stderr,"Error: Failed to open the dataset file - %s\n",argv[argc-2]);
		fprintf(stderr,"Please check the file name and restart this program.\n");
		fprintf(stderr,"The program will be terminated!\n");
		exit(-1);
	}

	// read query file
	read(argv[argc-1], queryset);
	if (queryset.size() == 0) {
		fprintf(stderr,"Error: Failed to open the query file - %s\n",argv[argc-1]);
		fprintf(stderr,"Please check the file name and restart this program.\n");
		fprintf(stderr,"The program will be terminated!\n");
		exit(-1);
	}

    //in case of we build the index and use it in the same program
    CGram gramgenupper(gl_qmax, false);
    CGram gramgenlower(gl_qmin, false);
    CCountFilter thefilter(&gramgenupper, gl_maxLen, PTHEAD_NUM);
    CSeqDB<> db_load(gl_maxLen, PTHEAD_NUM, &dataset, &gramgenupper, &gramgenlower, &thefilter);
    db_load.buildindex();

    //in case of we need to load the index from index file
    /*
    string index = argv[argc-2];
	index += ".ix";
    //build the upper and lower index, load the dataset and the .ix file
	CSeqDB<> db_load(gl_maxLen, PTHEAD_NUM, &dataset);
	db_load.load(index.c_str());
	CGram gramgenupper(db_load.gramGenUpper_len, false);
	CGram gramgenlower(db_load.gramGenLower_len, false);
	CCountFilter thefilter(&gramgenupper, gl_maxLen, db_load.gram_maxed);
	db_load.initParas((unsigned)gl_topk_f, &gramgenupper, &gramgenlower, &thefilter);
    */

	pThis = &db_load;
	long processednum = 0;

	for (unsigned i = 0; i < queryset.size(); i++) {//get each query
		if (queryset[i].size() < db_load.gramGenUpper_len) {
			fprintf(stdout, "We support queries with length larger than %d\n", db_load.gramGenUpper_len);
			fprintf(stdout, "Too small query size for [%s]\n", queryset[i].c_str());
			continue;
		}
        //cout the processing query index
        cout << "processing the " << i << "th" << endl;
		db_load.reset();//reset the db_load

		class TSINGHUA_CLIPSE_UTIL::TimeRecorder _time;
		CQuery query(queryset[i].c_str(), db_load.gramGenUpper, (unsigned)gl_topk);

		db_load.init_threshold(query);
		db_load.filter->setQueryLb(query);
		if (query.threshold != 0){
            pipe_knn(&query);
		}
		_time.check();
		query_time += _time.diffTime(0, 1);

		if (output) {
			fprintf(stdout, "Finish the query for [%s]:\n", queryset[i].c_str());
			while(!db_load.m_queue.empty()) {
				const queue_entry &entry = db_load.m_queue.top();
				fprintf(stdout, "%d: [%s]\t%d\n", entry.m_sid, dataset[entry.m_sid].c_str(), entry.m_dist);
				db_load.m_queue.pop();
			}
		}
		processednum += db_load.processed;
	}

	fprintf(stdout, "The average access number and query time: %d\t%-10.6f (sec)\n", (int)((float)processednum/queryset.size()), (float)(query_time)/(float)queryset.size());
}
