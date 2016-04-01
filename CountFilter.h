/*
 * CountFilter.h
 *
 *  Created on: Oct 22, 2010
 *      Author: xiaoliwang
 */

#ifndef _COUNTFILTER_H_
#define _COUNTFILTER_H_

#include "Gram.h"
#include "Query.h"
#include <vector>
#include <string>

class CCountFilter
{
public:
	int** tabUp;
	int** tabUpQuery;
	unsigned maxlen;
	unsigned maxed;
	unsigned maxged;

	CCountFilter(CGram* gram, const unsigned& data_maxlen, const unsigned& gram_maxed);
	~CCountFilter(void){
		for (unsigned i = 1; i < maxlen; i++) {
			if (this->tabUp[i] != NULL) {
				delete [] this->tabUp[i];
				this->tabUp[i] = NULL;
			}
			if (this->tabUpQuery[i] != NULL) {
				delete [] this->tabUpQuery[i];
				this->tabUpQuery[i] = NULL;
			}
		}
		delete [] this->tabUp;
		this->tabUp = NULL;
		delete [] this->tabUpQuery;
		this->tabUpQuery = NULL;
	}
	void setQueryLb(CQuery& query);
	void print(FILE *pf);
};

#endif
