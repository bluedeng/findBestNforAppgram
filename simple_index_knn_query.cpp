//2015-11, dtz
//for simple_index_knn_query
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

using namespace std;

int gl_qmax = 1;
int gl_topk = 5;
vector<string> dataset;
vector<string> queryset;
unsigned gl_maxLen = 0;
bool output = false;

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
	cout << "Usage: simple_index [OPTION] ${INPUTDB} ${QUERY}" << endl;
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
    ofstream fout((stemp + "_simple_index_analysis").c_str());
    //fout << "Gram_length, Index_time, Average_max_ed, Average_candidate_number, Average_query_time, Average_queryPerp_time, Average_all_time, Max_ed in results" << endl;
    fout << "Gram_length, Index_time, Average_candidate_number, Average_queryPerp_time, Least_time_recorded" << endl;

	//get the max_ed from file, be set as a threshold in the code
	int ed = log(get_ed(stemp)) / log(2) + 1;
    //int ed = get_ed(stemp) / 2;
	//int ed = 0.618 * get_ed(stemp) + 1;

	//this ed should be change along with the formula above
	if (ed == 0){
		cout << "We can not read the ed from the max_ed_file, please check it" << endl;
		exit(-1);
	}

    for (unsigned i = 1; ; i++){
        if (i >= 10 && countnumber >= 3){
            fout.close();
            break;
        }

        //set the gram length
        gl_qmax = i;

        double indextime = 0.0;
        double queryPreptime = 0.0;
        double querytime = 0.0;

        class TSINGHUA_CLIPSE_UTIL::TimeRecorder time1;

        CGram gramgenupper(gl_qmax, false);
        CCountFilter _fltTable(&gramgenupper, gl_maxLen, 1);
        CSeqDB<> db_load(gl_maxLen, 1, &dataset, &gramgenupper, NULL, &_fltTable);
        db_load.buildSimpleindex();

        time1.check();
        indextime += time1.diffTime(0, 1);
        index_time.push_back(indextime);

        max_ed.clear();
        candidate.clear();
        query_time.clear();
        queryPrep_time.clear();
        time_all.clear();

        pThis = &db_load;

        //store the max_ed in the query results
        int max_ed_in_results = 0;

        for (unsigned j = 0; j < queryset.size(); j++){

            queryPreptime = 0.0;
            querytime = 0.0;
            cout << "processing the " << j << "th" << endl;
            pThis->reset();

            //queryPrep: compose the query to grams and get the candidates
            class TSINGHUA_CLIPSE_UTIL::TimeRecorder time2;
            CQuery query(queryset[j].c_str(), &gramgenupper, (unsigned)gl_topk);

            query.threshold = ed;
            _fltTable.setQueryLb(query);

            pThis->theQuery = &query;
            //return time of querying for query j
			querytime += pThis->stringQuery(*(pThis->theQuery), ed);
            time2.check();
            queryPreptime += time2.diffTime(0, 1);
            //minus the query time
            queryPreptime -= querytime;

            /* //we only compute the number of the candidate and the filter time only
            max_ed_in_results = max_ed_in_results > pThis->m_queue.top().m_dist ? max_ed_in_results : pThis->m_queue.top().m_dist;
            max_ed[pThis->m_queue.top().m_dist] += 1;
            */
            candidate[pThis->processed] += 1;
            query_time.push_back(querytime);
            queryPrep_time.push_back(queryPreptime);
            time_all.push_back(querytime + queryPreptime);

            /*
			//just add it for verify
			if (pThis->processed - pThis->sizeofcandis >= 5 || pThis->sizeofcandis - pThis->processed >= 5){
				cout << "There may be something wrong with the processed number" << endl;
				fout << "There may be something wrong with the processed number" << endl;
			}
			*/

            if (output){
                cout << "Finish the query for: " << queryset[j] << endl;
                /*
                while (!pThis->m_queue.empty()){
                    const queue_entry &entry = pThis->m_queue.top();
                    cout << entry.m_sid << " " << dataset[entry.m_sid] << " " << entry.m_dist << endl;
                    //fdetail << dataset[entry.m_sid] << " " << entry.m_dist << endl;
                    pThis->m_queue.pop();
                }
                */
            }
        }

        //analysis and static
        map<unsigned, unsigned>::iterator iter;
        unsigned temp = 0;
        unsigned len = queryset.size();

        /*
        for (iter = max_ed.begin(); iter != max_ed.end(); iter++)
            temp += iter->first * iter->second;
        max_ed_average.push_back((double)temp / len);
        */

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

        /*
        if (min_time < time_all_average.back())
            countnumber++;
        else {
            min_time = time_all_average.back();
            countnumber = 0;
        }
        */
        if (min_time < queryPrep_time_average.back())
            countnumber++;
        else {
            min_time = queryPrep_time_average.back();
            countnumber = 0;
        }

        /*
        //Gram_length, Index_time, Average_max_ed, Average_candidate_number, Average_query_time, Average_queryPerp_time, Average_all_time
        fout << i << "  " << index_time.back() << "  " << max_ed_average.back() << "  " << candidata_average.back() << "  "
             << query_time_average.back() << "  " << queryPrep_time_average.back() << "  " << time_all_average.back() << "  "
             << max_ed_in_results << endl;
        */
        //Gram_length, Index_time, Average_candidate_number, Average_queryPrep_time, Least_time_recorded
        fout << i << "  " << index_time.back() << "  " << candidata_average.back()
        << "  " << queryPrep_time_average.back() << "  " << min_time << endl;

        cout << "Processed number and average time: " << candidata_average.back() << "    " << time_all_average.back() << endl;
    }
	return 0;
}
