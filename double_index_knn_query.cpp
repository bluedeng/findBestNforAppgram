//2015-11, dtz
//Using the result of simple_index_knn_query, build the lower index using n1
//compute the best n2 for the upper index
//for double_index_knn_query
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
#include <math.h>
#include <sstream>
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
bool output = false;

double indextime = 0.0;
double queryPreptime = 0.0;
double querytime = 0.0;

map<unsigned, unsigned> candidate;
vector<double> candidate_average;
map<unsigned, unsigned> result;
vector<double> result_average;
map<unsigned, unsigned> max_ed;
vector<double> max_ed_average;
vector<double> index_time;
//vector<double> queryPrep_time;
//vector<double> queryPrep_time_average;
vector<double> query_time;
vector<double> query_time_average;
vector<double> time_all;
vector<double> time_all_average;

//usage instruction
void usage(){
	cout << "************************************************" << endl;
	cout << "First, build the index for the sequence search  " << endl;
	cout << "Second, do the search process and output results" << endl;
	cout << "------------------------------------------------" << endl;
	cout << "Usage: double_index [OPTION] ${INPUTDB} ${QUERY}" << endl;
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
void read(string filename, vector<string>& data, bool b){
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

//string compare
bool isCompare(string s1, string s2){
	for (int i = 0; i < 3; i++)
		if (s1[i] != s2[i])
			return false;

	return true;
}

int get_ed(string s){
	//the directory shoud be change for difference place
	ifstream fin((s.substr(0, 34) + "max_ed_file").c_str());

	string line;
	while (!fin.eof()){
		getline(fin, line);
		if (isCompare(line, s.substr(34, 3))){
			getline(fin, line);
			stringstream ss;
			ss << line;
			int ed;
			ss >> ed;
			fin.close();

			return ed;
		}
		getline(fin, line);
	}
	fin.close();
	return 0;
}

double calculateAverage(map<unsigned, unsigned>& m, unsigned l){
    unsigned temp = 0;
    map<unsigned, unsigned>::iterator iter;
    for (iter = m.begin(); iter != m.end(); iter++)
        temp += iter->first * iter->second;
    return (double)temp / l;
}

int main(int argc, char* argv[]){
	if (argc < 3){
		usage();
		exit(-1);
	}

	parseOptions(argc, argv);

	//read the data file
	read(argv[argc - 2], dataset, true);
	if (dataset.size() == 0){
		cerr << "Error: Failed to open the data file: " << argv[argc - 2] << endl;
		cerr << "Please check the file name and restart this program" << endl;
		cerr << "The program will be terminated!" << endl;
		exit(-1);
	}

	//read the query file
	read(argv[argc - 1], queryset, false);
	if (queryset.size() == 0){
		cerr << "Error: Failed to open the data file: " << argv[argc - 1] << endl;
		cerr << "Please check the file name and restart this program" << endl;
		cerr << "The program will be terminated!" << endl;
		exit(-1);
	}

	unsigned countnumber = 0;
	double min_time = 1000.0;

	string stemp = argv[argc - 1];
	ofstream fout((stemp + "_double_index_analysis").c_str());
	//ofstream fdetail((stemp + "_double_TS").c_str());

	cout << "Please input the n of the lower index" << endl;
	cin >> gl_qmin;

	fout << "The lower index n is " << gl_qmin << endl;
	fout << endl;
	fout << "Upper_gram_length, Index_time, Average_max_ed, Average_candidate_number, ";
	fout << "Average_result_number, Average_query_time/*, Average_queryPerp_time*/, Average_all_time" << endl;

	//get the max_ed from file, be set as a threshold in the code
	//int ed = log(get_ed(stemp)) / log(2) + 1;
    //int ed = get_ed(stemp) / 2;
	int ed = 0.618 * get_ed(stemp) + 1;

	//this ed should be change along with the formula above
	if (ed == 0){
		cout << "We can not read the ed from the max_ed_file, please check it" << endl;
		exit(-1);
	}

	//outer loop: for index building
	for (unsigned i = 2 * gl_qmin; ; i++){
		if (countnumber >= 3){
			fout.close();
			//fdetail.close();
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
		result.clear();
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

			//pThis->init_threshold(query);
			query.threshold = ed;
			pThis->filter->setQueryLb(query);
			//time2.check();
			//queryPreptime += time2.diffTime(0, 1);

			if (pipe_knn(&query)){
                time2.check();
				time_all.push_back(time2.diffTime(0, 1));

				max_ed[pThis->m_queue.top().m_dist] += 1;
				candidate[pThis->processed] += 1;
				result[pThis->m_queue.size()] += 1;
				query_time.push_back(querytime);
				//queryPrep_time.push_back(queryPreptime);
				//can't compute the queryPrep_time for each query, so comput the time_all and query_time instead

				if (output){
					cout << "Finish the query for: " << queryset[j] << endl;
					//fdetail << "Query: " << queryset[j] << endl;
					while (!pThis->m_queue.empty()){
						const queue_entry &entry = pThis->m_queue.top();
						cout << entry.m_sid << " " << dataset[entry.m_sid] << " " << entry.m_dist << endl;
						//fdetail << dataset[entry.m_sid] << " " << entry.m_dist << endl;
						//fdetail << entry.m_dist << " ";
						pThis->m_queue.pop();
					}
				}
				//fdetail << endl;
			}
		}

		//analysis and static
		unsigned temp = 0;
		unsigned len = queryset.size();

		max_ed_average.push_back(calculateAverage(max_ed, len));
		candidate_average.push_back(calculateAverage(candidate, len));
        result_average.push_back(calculateAverage(result, len));

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
			countnumber = 0;
		}

		//Gram_length, Index_time, Average_max_ed, Average_candidate_number,
		//Average_query_time, /*Average_queryPerp_time,*/ Average_all_time
		fout    << i << "  " << index_time.back() << "  " << max_ed_average.back() << "  "
                << candidate_average.back() << "  " << result_average.back() << "  "
                << query_time_average.back() /*<< "  " << queryPrep_time_average.back()*/ << "  " << time_all_average.back() << endl;

		cout << "Processed number and average time: " << candidate_average.back() << "    " << time_all_average.back() << endl;
	}
	return 0;
}
