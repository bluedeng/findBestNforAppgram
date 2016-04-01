/*
 * GramList.h
 *
 *  Created on: Sep 14, 2010
 *      Author: xiaoliwang
 */
#ifndef _GREMLIST_H_
#define _GREMLIST_H_

#include <vector>

using namespace std;

template <typename InvList = vector<unsigned> >
class CGramList {
 protected:
  vector<unsigned> invertedList;
  
 public:
  CGramList() {};
  InvList* getArray() { return &invertedList; }
  ~CGramList() {};
  void free() { delete this; }
  void clear() { }
};

template<class InvList = vector<unsigned> >
class QueryGramList {
public:
  unsigned gramCode;
  CGramList<InvList>* gl;  

  QueryGramList(unsigned gramCode, CGramList<InvList>* gl)
    : gramCode(gramCode), gl(gl) {}
};

#endif
