/*
 * Gram.h
 *
 *  Created on: Sep 14, 2010
 *      Author: xiaoliwang
 */

#ifndef _GREM_H_
#define _GREM_H_

#include <vector>
#include <set>
#include <string>

#ifdef _WIN32
#include <functional>
#else
#include <tr1/functional>
#endif

typedef unsigned char uchar;

typedef enum
{
  GT_FIXED,
  GT_WORDS
} GramType;

#ifdef _WIN32
class hash_win32
{
public:
  std::size_t operator()(const std::string& s) const
  {
    const char* first = s.data();
    std::size_t length = s.length();
    std::size_t result = 0;
    for (; length > 0; --length)
      result = (result * 131) + *first++;
    return result;
  }
};
#endif

class CGram
{
private:
	unsigned n;

protected:
#ifdef _WIN32
  static const hash_win32 hashString;
#else
  static const std::tr1::hash<std::string> hashString;
#endif

	static const uchar PREFIXCHAR = 156; // pound
	static const uchar SUFFIXCHAR = 190; // yen
	GramType gramType;

public:
	// Using decomposition with prefix and suffix or not
	bool prePost;

	CGram(unsigned length, bool usePrePost = true):n(length),prePost(usePrePost){}
	~CGram(void){}

	unsigned getGramLength() const;

    //substr->grams, hashcode->res
	void decompose(
    const std::string& s,
	std::vector<std::string>& grams,
    std::vector<unsigned>& res,
    uchar st = PREFIXCHAR,
    uchar en = SUFFIXCHAR)
    const;

	void decompose(
    const std::string& s,
    std::vector<unsigned>& res,
    uchar st = PREFIXCHAR,
    uchar en = SUFFIXCHAR)
    const;

	void decompose(
    const std::string& s,
    std::set<unsigned>& res,
    uchar st = PREFIXCHAR,
    uchar en = SUFFIXCHAR)
    const;
};

#endif
