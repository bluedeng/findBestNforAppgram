/*
 * Gram.cpp
 *
 *  Created on: Sep 14, 2010
 *      Author: xiaoliwang
 */

#include "Gram.h"

using namespace std;
using namespace tr1;

#if _WIN32
const hash_win32 CGram::hashString = hash_win32();
#else
const hash<string> CGram::hashString = hash<string>();
#endif

//return n stands for n-grams
unsigned CGram::getGramLength() const
{
	return this->n;
}

void CGram::decompose(
  const string &s,
  vector<string> &grams,
  vector<unsigned> &res,
  uchar st,
  uchar en)
  const
{
 if(prePost) {
    const string sPad = string(n - 1, st) + s + string(n - 1, en);
    for (unsigned i = 0; i < s.length() + n - 1; i++) {
      string substring = sPad.substr(i, n);
	  grams.push_back(substring);
      res.push_back(hashString(substring));
    }
  }
  else {
    if(s.length() < n) {
	  grams.push_back(s);
      res.push_back(hashString(s));
      return;
    }

    for (unsigned i = 0; i < s.length() - n + 1; i++) {
      string substring = s.substr(i, n);
	  grams.push_back(substring);
      res.push_back(hashString(substring));
    }
  }
}

void CGram::decompose(
  const string &s,
  vector<unsigned> &res,
  uchar st,
  uchar en)
  const
{
 if(prePost) {
    const string sPad = string(n - 1, st) + s + string(n - 1, en);
    for (unsigned i = 0; i < s.length() + n - 1; i++) {
      res.push_back(hashString(sPad.substr(i, n)));
    }
  }
  else {
    if(s.length() < n) {
      res.push_back(hashString(s));
      return;
    }

    for (unsigned i = 0; i < s.length() - n + 1; i++) {
      res.push_back(hashString(s.substr(i, n)));
    }
  }
}

void CGram::decompose(
  const string &s,
  set<unsigned> &res,
  uchar st,
  uchar en)
  const
{
 if(prePost) {
    const string sPad = string(n - 1, st) + s + string(n - 1, en);
    for (unsigned i = 0; i < s.length() + n - 1; i++) {
      res.insert(hashString(sPad.substr(i, n)));
    }
  }
  else {
    if(s.length() < n) {
      res.insert(hashString(s));
      return;
    }

    for (unsigned i = 0; i < s.length() - n + 1; i++) {
      res.insert(hashString(s.substr(i, n)));
    }
  }
}
