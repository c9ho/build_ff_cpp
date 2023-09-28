#ifndef TOOL_H
#define TOOL_H

#include <iostream>
#include <vector>

std::vector<std::string> split(std::string input, std::string delimiter = " ");
double innerProduct(const double vec_1[3], const double vec_2[3]);
void outerProduct(const double vec_1[3], const double vec_2[3], double result[3]);

#endif /* tool.h */