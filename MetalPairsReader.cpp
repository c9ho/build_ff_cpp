#include "MetalPairsReader.h"

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

MetalPairsReader::MetalPairsReader() {
    fstream file("metal_pairs");
    string line;
    istringstream iss;
    string atom1, atom2;

    cout << "** Reading metal_pairs **" << endl;
    if (!file.is_open()) {
        throw runtime_error("** Error: cannot open metal_pairs **");
    }

    getline(file, line);
    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        iss.clear();
        iss.str(line);
        iss >> atom1 >> atom2;
        metalPairs.push_back({atom1, atom2});
    }

    file.close();
}

ostream& operator<<(ostream& stream, const array<std::string, 2>& array) {
    stream << "(" << array[0]
           << ", " << array[1]
           << ")";

    return stream;
}