#include <iostream>
#include <vector>
using namespace std;

vector<string> split(string input, string delimiter = " ") {
    vector<string> splitted;
    size_t pos = 0;
    string token;
    while (pos != std::string::npos) {
        pos = input.find(delimiter);
        if (pos != 0) {
            token = input.substr(0, pos);
            if (token != "") {
                splitted.push_back(token);
            }
        }
        input.erase(0, pos + delimiter.length());
    }
    return splitted;
}

double innerProduct(const double vec_1[3], const double vec_2[3]) {
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2]);
}

void outerProduct(const double vec_1[3], const double vec_2[3], double result[3]) {
    result[0] = vec_1[1] * vec_2[2] - vec_1[2] * vec_2[1];
    result[1] = vec_1[2] * vec_2[0] - vec_1[0] * vec_2[2];
    result[2] = vec_1[0] * vec_2[1] - vec_1[1] * vec_2[0];
}