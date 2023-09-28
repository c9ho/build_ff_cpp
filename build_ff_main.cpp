#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "Processor.h"

using namespace std;

vector<string> findFileNames(int argumentCount, char** argumentVector) {
    vector<string> fileNames;

    for (auto i = 1; i < argumentCount; i++) {
        fileNames.push_back(argumentVector[i]);
    }
    return fileNames;
}

int main(int argc, char** argv) {
    const chrono::time_point now = chrono::system_clock::now();
    const time_t timeNow = chrono::system_clock::to_time_t(now);

    Processor processor(findFileNames(argc, argv));
    processor.processMassPart();
    processor.processAtomPart();
    processor.processBondPart();
    processor.processAnglePart();
    processor.processDihedralPart();
    processor.processImproperPart();
    processor.processVDWPart();
    processor.processBoxInfo();

    ofstream headerFile("HEADER");

    headerFile << "Created on " << ctime(&timeNow) << "\n";
    headerFile << processor.atomLine << "\n";
    headerFile << processor.bondCountLine << "\n";
    headerFile << processor.angleCountLine << "\n";
    headerFile << processor.dihedralCountLine << "\n";
    headerFile << processor.improperCountLine << "\n";
    headerFile << "\n";
    headerFile << processor.massLine << "\n";
    headerFile << processor.bondIndexLine << "\n";
    headerFile << processor.angleIndexLine << "\n";
    headerFile << processor.dihedralIndexLine << "\n";
    headerFile << processor.improperIndexLine << "\n";
    headerFile << processor.aLine << "\n";
    headerFile << processor.bLine << "\n";
    headerFile << processor.cLine << "\n";
    headerFile << processor.tiltLine << "\n"
               << "\n";

    headerFile.close();

    system("cat HEADER MASS ATOM BOND ANGLE DIHEDRAL IMPROPER > data.system");
    system("cat bond_coeff angle_coeff dihedral_coeff improper_coeff > coefficients");
    system("rm HEADER MASS ATOM BOND ANGLE DIHEDRAL IMPROPER bond_coeff angle_coeff dihedral_coeff improper_coeff");

    return 0;
}
