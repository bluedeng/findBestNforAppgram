//#2015-11, dtz
//Using the result of simple_index_knn_query and the idea of double_index_knn_query
//building the lower index using n1
//using the lower bound of upper index of 2 * n1 for n2
//using the lower bound of upper index of 2 * best_n2 for n3, while n1 still the lower index
//increasing the nk.... until nk >= min(query.length)
//nk <= query_minLength

#include "Time.h"
#include "Gram.h"
#include "CountFilter.h"
#include "SeqDB.h"
#include "Query.h"

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <fstream>
#include <pthread.h>
#include <unistd.h>

using namespace std;

int gl_qmax;
int gl_qmin;
int gl_topk = 5;
vector<string> dataset;
vector<string> queryset;
unsigned gl_maxLen = 0;
unsigned query_minLen = 1000;
bool output = false;

double indextime = 0.0;
//double queryPreptime = 0.0;
double querytime = 0.0;

map<unsigned, unsigned> candidate;
vector<double> candidata_average;
map<unsigned, unsigned> max_ed;
vector<double> max_ed_average;
vector<double> index_time;
//vector<double> queryPrep_time;
//vector<double> queryPrep_time_average;
vector<double> query_time;
vector<double> query_time_average;
vector<double> time_all;
vector<double> time_all_average;
//recording the n-value for each level
vector<unsigned> nk;

//usage instruction
void usage(){
	cout << "************************************************" << endl;
	cout << "First, build the index for the sequence search  " << endl;
	cout << "Second, do the search process and output results" << endl;
	cout << "------------------------------------------------" << endl;
	cout << "Usage: multi_index [OPTION] ${INPUTDB} ${QUERY} " << endl;
	cout << "-k val     set the value for top-k heap         " << endl;
	cout << "-o         output the result set of sequence ids" << endl;
	cout << "------------------------------------------------" << endl;
	cout << "************************************************" << endl;
}

//get the input parameters: output
void parseOptions(int argc, char* argv[]){
	if (argc > 3)
	    for (int i = 1; i < argc - 2; i++){
            if (strncmp(argv[i], "-k", 2) == 0){
                gl_topk = atoi(argv[i + 1]);
                i++;
            }
            if (strncmp(argv[i], "-o", 2) == 0)
                output = true;
	    }
}

//read the input dataset file
void read(string filename, vector<string>& data){
	ifstream _fs;
	_fs.open(filename.c_str(), ios::in);
	if (_fs.fail()){
		cerr << "Error: Failed to open the data file: " << filename << endl;
		cerr << "Please check the file name and restart this program" << endl;
		cerr << "The program will be terminated!" << endl;
		exit(-1);
	}

	string line;
	while (!_fs.eof()){
		getline(_fs, line);

		if (line.empty())
			continue;

		for (unsigned i = 0; i < line.length(); i++)
			if (line[i] >= 'A' && line[i] <= 'Z')
				line[i] += 32;

		gl_maxLen = gl_maxLen > line.length() ? gl_maxLen : line.length();

		data.push_back(line);
	}

	_fs.close();
}

void getQuery_minlength(vector<string>& data){
	for (unsigned i = 0; i < data.size(); i++)
		query_minLen = query_minLen < data[i].length() ? query_minLen : data[i].length();
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
		    class TSINGHUA_CLIPSE_UTIL::TimeRecorder _time;
			// Try post precessing here with more approximate bounds
			pThis->knn_postprocess();
			pThis->old_version_knn_postprocess();
			// Then end of whole processing here
            _time.check();
            querytime += _time.diffTime(0, 1);
			return 1;
		}
		sleep(0.0001);
	}
	return 0;
}

int main(int argc, char* argv[]){
	if (argc < 3){
		usage();
		exit(-1);
	}

	parseOptions(argc, argv);

	//read the data file
	read(argv[argc - 2], dataset);
	if (dataset.size() == 0){
		cerr << "Error: Failed to open the data file: " << argv[argc - 2] << endl;
		cerr << "Please check the file name and restart this program" << endl;
		cerr << "The program will be terminated!" << endl;
		exit(-1);
	}

	//read the query file
	read(argv[argc - 1], queryset);
	if (queryset.size() == 0){
		cerr << "Error: Failed to open the data file: " << argv[argc - 1] << endl;
		cerr << "Please check the file name and restart this program" << endl;
		cerr << "The program will be terminated!" << endl;
		exit(-1);
	}

	//get the min length of query
	getQuery_minlength(queryset);

	string stemp = argv[argc - 1];
	ofstream fout((stemp + "_multi_index_analysis").c_str());

	cout << "Please input the n of the lower index" << endl;
	cin >> gl_qmin;

	fout << "The lower index n is " << gl_qmin << endl;
	fout << "The shortest query_string's length is " << query_minLen << endl;
	fout << endl;
	fout << "Upper_gram_length, Index_time, Average_max_ed, Average_candidate_number, Average_query_time/*, Average_queryPerp_time*/, Average_all_time" << endl;
    fout << endl;
	//for multi_index test, choosing best nk for each level
	unsigned gl_temp = gl_qmin;
	while (gl_temp <= query_minLen){
		nk.push_back(gl_temp);
		fout << "The best n-value for " << nk.size() << "th level index is " << nk.back() << endl;
		fout << "Starting the " << nk.size() + 1 << "th level index test (bottom up!)" << endl;

		unsigned countnumber = 0;
		double min_time = 1000.0;

		max_ed_average.clear();
		candidata_average.clear();
		query_time_average.clear();
		time_all_average.clear();

		//outer loop: for index building
		for (unsigned i = 2 * gl_temp; i <= query_minLen; i++){
			if (countnumber >= 3){
				fout.close();
				break;
			}

			indextime = 0.0;
			gl_qmax = i;

			class TSINGHUA_CLIPSE_UTIL::TimeRecorder time1;

			CGram gramgenupper(gl_qmax, false);
			CGram gramgenlower(gl_qmin, false);
			CCountFilter _fltTable(&gramgenupper, gl_maxLen, PTHEAD_NUM);
			CSeqDB<> db_load(gl_maxLen, PTHEAD_NUM, &dataset, &gramgenupper, &gramgenlower, &_fltTable);
			db_load.buildindex();

			time1.check();
			indextime += time1.diffTime(0, 1);
			index_time.push_back(indextime);

			max_ed.clear();
			candidate.clear();
			query_time.clear();
			//queryPrep_time.clear();
			time_all.clear();

			pThis = &db_load;

			//inner loop: for querying
			for (unsigned j = 0; j < queryset.size(); j++){
				//queryPreptime = 0.0;
				querytime = 0.0;
				cout << "Processing the " << j << "th" << endl;
				pThis->reset();

				class TSINGHUA_CLIPSE_UTIL::TimeRecorder time2;
				CQuery query(queryset[j].c_str(), &gramgenupper, (unsigned)gl_topk);
				pThis->init_threshold(query);
				pThis->filter->setQueryLb(query);
				//time2.check();
				//queryPreptime += time2.diffTime(0, 1);

				if (pipe_knn(&query)){
                    time2.check();
					time_all.push_back(time2.diffTime(0, 1));

					max_ed[pThis->m_queue.top().m_dist] += 1;
					candidate[pThis->processed] += 1;
					query_time.push_back(querytime);
					//queryPrep_time.push_back(queryPreptime);
					//can't compute the queryPrep_time for each query, so comput the time_all and query_time instead

					if (output){
						cout << "Finish the query for: " << queryset[j] << endl;
						while (!pThis->m_queue.empty()){
							const queue_entry &entry = pThis->m_queue.top();
							cout << entry.m_sid << " " << dataset[entry.m_sid] << " " << entry.m_dist << endl;
							pThis->m_queue.pop();
						}
					}
				}
			}

			//analysis and static
			map<unsigned, unsigned>::iterator iter;
			unsigned temp = 0;
			unsigned len = queryset.size();

			for (iter = max_ed.begin(); iter != max_ed.end(); iter++)
				temp += iter->first * iter->second;
			max_ed_average.push_back((double)temp / len);

			temp = 0;
			for (iter = candidate.begin(); iter != candidate.end(); iter++)
				temp += iter->first * iter->second;
			candidata_average.push_back((double)temp / len);

			double t1 = 0.0;
			//double t2 = 0.0;
			double t3 = 0.0;
			for (temp = 0; temp < len; temp++){
				t1 += query_time[temp];
				//t2 += queryPrep_time[temp];
				t3 += time_all[temp];
			}
			query_time_average.push_back(t1 / len);
			//queryPrep_time_average.push_back(t2 / len);
			time_all_average.push_back(t3 / len);

			if (min_time < time_all_average.back())
				countnumber++;
			else {
				min_time = time_all_average.back();
				gl_temp = i;
				countnumber = 0;
			}

			//Gram_length, Index_time, Average_max_ed, Average_candidate_number,
			//Average_query_time, /*Average_queryPerp_time,*/ Average_all_time
			fout << i << "  " << index_time.back() << "  " << max_ed_average.back() << "  " << candidata_average.back() << "  "
				<< query_time_average.back() /*<< "  " << queryPrep_time_average.back()*/ << "  " << time_all_average.back() << endl;

			cout << "Processed number and average time: " << candidata_average.back() << "    " << time_all_average.back() << endl;
		}
	}
	return 0;
}
