#ifndef __PGUI_HELPER_H__
#define __PGUI_HELPER_H__

#include <string>
#include <vector>

bool fileExists(const std::string & f);

std::vector<std::string> char_split( const std::string & s , const char c , bool empty = true );
std::vector<std::string> quoted_char_split( const std::string & s , const char c , bool empty );
std::string int2str(int);

#endif
