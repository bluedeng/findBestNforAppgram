/*
 * CountFilter.cpp
 *
 *  Created on: Oct 22, 2010
 *      Author: xiaoliwang
 */

#include "CountFilter.h"

using namespace std;

CCountFilter::CCountFilter(CGram* gram, const unsigned& data_maxlen, const unsigned& gram_maxed)
{
	// lower bound
	unsigned j, k;
	int gramLen = (int)gram->getGramLength();
	this->maxlen = data_maxlen+1;
	this->maxged = gram_maxed;
	this->tabUp = new int* [this->maxlen];
	this->tabUpQuery = new int* [this->maxlen];
	this->tabUp[0] = NULL;
	this->tabUpQuery[0] = NULL;
	for (j = 1; j < this->maxlen; j++) {
		this->tabUp[j] = new int[gram_maxed];
		this->tabUpQuery[j] = new int[gram_maxed];
		for (k = 0; k < gram_maxed; k++) {
			int affectGramNum = (gramLen - 2*k) > 1 ? (gramLen - 2*k) : 1;
			this->tabUp[j][k] = affectGramNum + (gramLen-k)*(j-1);
		}
	}
}

void CCountFilter::setQueryLb(CQuery &query)
{
	unsigned i, j;
	this->maxed = query.threshold+1;
	for (i = 1; i < this->maxed; i++) {
		for (j = 0; j < this->maxged; j++) {
			this->tabUpQuery[i][j] = ((int)(query.gramsize) - this->tabUp[i][j]) > 0 ? ((int)(query.gramsize) - this->tabUp[i][j]) : 0;
		}
	}
}

void CCountFilter::print(FILE *pf)
{
	unsigned i, j;
	fprintf(pf, "\n****************Mismatch filter table**********************\n");
	fprintf(pf, "The upper table: \n");
	for (i = 1; i < this->maxed; i++) {
		fprintf(pf, "(%d)\t[", i);
		for (j = 0; j < this->maxged; j++)
			fprintf(pf, " %d", this->tabUpQuery[i][j]);
		fprintf(pf, "]\n");
	}
	fprintf(pf, "\n");
}
