//2015-11, dtz
//testing file, remove the multi-threads
//for fixing some bugs in double_index_knn_query files

//Using the result of simple_index_knn_query, build the lower index using n1
//compute the best n2 for the upper index

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

using namespace std;

int gl_qmax;
int gl_qmin;
int gl_topk = 5;
vector<string> dataset;
vector<string> queryset;
unsigned dataset_maxLen = 0;
unsigned query_maxLen = 0;
bool output = false;

double indextime = 0.0;
double queryPreptime = 0.0;
double querytime = 0.0;

map<unsigned, unsigned> candidate;
vector<double> candidata_average;
map<unsigned, unsigned> max_ed;
vector<double> max_ed_average;
vector<double> index_time;
vector<double> queryPrep_time;
vector<double> queryPrep_time_average;
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
	unsigned maxLen = 0;
	while (!_fs.eof()){
		getline(_fs, line);

		if (line.empty())
			continue;

		for (unsigned i = 0; i < line.length(); i++)
			if (line[i] >= 'A' && line[i] <= 'Z')
				line[i] += 32;

		maxLen = maxLen > line.length() ? maxLen : line.length();

		data.push_back(line);
	}

	if (b)
		dataset_maxLen = maxLen;
	else
		query_maxLen = maxLen;

	_fs.close();
}

static CSeqDB<> *pThis;

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
	ofstream fout((stemp + "_double_test").c_str());

	cout << "Please input the n of the lower index" << endl;
	cin >> gl_qmin;

	fout << "The lower index n is " << gl_qmin << endl;
	fout << endl;
	fout << "Upper_gram_length, Index_time, Average_max_ed, Average_candidate_number, Average_query_time, Average_queryPerp_time, Average_all_time" << endl;

	//outer loop: for index building
	for (unsigned i = 2 * gl_qmin;; i++){
		if (countnumber >= 5){
			fout.close();
			break;
		}

		indextime = 0.0;
		gl_qmax = i;

		class TSINGHUA_CLIPSE_UTIL::TimeRecorder time1;

		CGram gramgenupper(gl_qmax, false);
		CGram gramgenlower(gl_qmin, false);
		CCountFilter _fltTable(&gramgenupper);
		CSeqDB<> db_load(&dataset, &gramgenupper, &gramgenlower, &_fltTable, dataset_maxLen, query_maxLen);
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
			queryPreptime = 0.0;
			querytime = 0.0;
			cout << "Processing the " << j << "th" << endl;
			pThis->reset();

			class TSINGHUA_CLIPSE_UTIL::TimeRecorder time2;
			CQuery query(queryset[j].c_str(), &gramgenupper, (unsigned)gl_topk);
			_fltTable.setQueryLb(query);
			pThis->theQuery = &query;
			pThis->accumulateFrequency(1);
			time2.check();
			queryPreptime += time2.diffTime(0, 1);

			class TSINGHUA_CLIPSE_UTIL::TimeRecorder time3;
			pThis->knn_pipesearch();
			time3.check();
			querytime += time3.diffTime(0, 1);

			max_ed[pThis->m_queue.top().m_dist] += 1;
			candidate[pThis->processed] += 1;
			query_time.push_back(querytime);
			queryPrep_time.push_back(queryPreptime);
			time_all.push_back(querytime + queryPreptime);

			if (output){
				cout << "Finish the query for: " << queryset[j] << endl;
				while (!pThis->m_queue.empty()){
					const queue_entry &entry = pThis->m_queue.top();
					cout << entry.m_sid << " " << dataset[entry.m_sid] << " " << entry.m_dist << endl;
					pThis->m_queue.pop();
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
		double t2 = 0.0;
		double t3 = 0.0;
		for (temp = 0; temp < len; temp++){
			t1 += query_time[temp];
			t2 += queryPrep_time[temp];
			t3 += time_all[temp];
		}
		query_time_average.push_back(t1 / len);
		queryPrep_time_average.push_back(t2 / len);
		time_all_average.push_back(t3 / len);

		if (min_time < time_all_average.back())
			countnumber++;
		else {
			min_time = time_all_average.back();
			countnumber = 0;
		}

		//Gram_length, Index_time, Average_max_ed, Average_candidate_number,
		//Average_query_time, Average_queryPerp_time, Average_all_time
		fout << i << "  " << index_time.back() << "  " << max_ed_average.back() << "  " << candidata_average.back() << "  "
			<< query_time_average.back() << "  " << queryPrep_time_average.back() << "  " << time_all_average.back() << endl;

		cout << "Processed number and average time: " << candidata_average.back() << "    " << time_all_average.back() << endl;
	}
	return 0;
}
