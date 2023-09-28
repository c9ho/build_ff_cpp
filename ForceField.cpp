#include "ForceField.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include "ForceFieldAngle.h"
#include "ForceFieldAtom.h"
#include "ForceFieldBond.h"
#include "ForceFieldDihedral.h"
#include "ForceFieldImproper.h"
#include "tool.h"

using namespace std;

ForceField::ForceField() {
    ifstream forceFieldFile = checkFileAndOpen();

    vector<string> allForceFieldLines;
    string line;
    istringstream iss;
    string title;
    int number;

    while (getline(forceFieldFile, line)) {
        allForceFieldLines.push_back(line);
    }

    int pointer = -1;
    for (auto& line : allForceFieldLines) {
        pointer++;
        if (line.starts_with("ATOMS")) {
            number = getNumber(line);
            prepareAtoms(allForceFieldLines, pointer, number);
        }
        if (line.starts_with("BONDS")) {
            number = getNumber(line);
            prepareBonds(allForceFieldLines, pointer, number);
        }
        if (line.starts_with("ANGLES")) {
            number = getNumber(line);
            prepareAngles(allForceFieldLines, pointer, number);
        }
        if (line.starts_with("DIHEDRALS")) {
            number = getNumber(line);
            prepareDihedrals(allForceFieldLines, pointer, number);
        }
        if (line.starts_with("IMPROPERS")) {
            number = getNumber(line);
            prepareImpropers(allForceFieldLines, pointer, number);
        }
        if (line.starts_with("VDW")) {
            number = getNumber(line);
            prepareVDWs(allForceFieldLines, pointer, number);
        }
    }
    forceFieldFile.close();
}

void ForceField::prepareVDWs(const vector<string>& allForceFieldLines, const int pointer, const int numberOfVDWs) {
    string line;
    istringstream iss;

    for (int i = 1; i <= numberOfVDWs; i++) {
        iss.clear();
        iss.str(allForceFieldLines[pointer + i]);
        ForceFieldVDW forceFieldVDW;

        iss >> forceFieldVDW;
        vdws.push_back(forceFieldVDW);
    }
}

void ForceField::prepareImpropers(const vector<string>& allForceFieldLines, const int pointer, const int numberOfImpropers) {
    string line;
    istringstream iss;

    for (int i = 1; i <= numberOfImpropers; i++) {
        iss.clear();
        iss.str(allForceFieldLines[pointer + i]);
        ForceFieldImproper forceFieldImproper;

        iss >> forceFieldImproper;
        impropers.push_back(forceFieldImproper);
    }
}

void ForceField::prepareDihedrals(const std::vector<std::string>& allForceFieldLines, const int pointer, const int numberOfDihedrals) {
    string line;
    istringstream iss;

    for (int i = 1; i <= numberOfDihedrals; i++) {
        iss.clear();
        iss.str(allForceFieldLines[pointer + i]);
        ForceFieldDihedral forceFieldDihedral;

        iss >> forceFieldDihedral;
        dihedrals.push_back(forceFieldDihedral);
    }
}

void ForceField::prepareAngles(const vector<string>& allForceFieldLines, const int pointer, const int numberOfAngles) {
    vector<string> splitted;

    for (int i = 1; i <= numberOfAngles; i++) {
        splitted = split(allForceFieldLines[pointer + i], " ");
        string* atom1 = &splitted[0];
        string* atomC = &splitted[1];
        string* atom2 = &splitted[2];
        string* functional = &splitted[3];
        double forceConstant = stod(splitted[4]);
        double eqAngle = stod(splitted[5]);

        ForceFieldAngle forceFieldAngle(*atom1, *atomC, *atom2, *functional, forceConstant, eqAngle);
        if (splitted.size() >= 9) {
            string restrictSign = splitted[8];
            transform(restrictSign.begin(), restrictSign.end(), restrictSign.begin(), ::toupper);
            if (restrictSign == "R") {
                forceFieldAngle.isRestricted = true;
                forceFieldAngle.restrictLowerAngle = stod(splitted[9]);
                forceFieldAngle.restrictUpperAngle = stod(splitted[10]);
            } else {
                throw runtime_error("Invalid restriction setting.");
            }
        }
        angles.push_back(forceFieldAngle);
    }
}

void ForceField::prepareBonds(const vector<string>& allForceFieldLines, const int pointer, const int numberOfBonds) {
    string line;
    istringstream iss;

    for (int i = 1; i <= numberOfBonds; i++) {
        iss.clear();
        iss.str(allForceFieldLines[pointer + i]);
        ForceFieldBond forceFieldBond;

        iss >> forceFieldBond;
        bonds.push_back(forceFieldBond);
    }
}

void ForceField::prepareAtoms(const vector<string>& allForceFieldLines, const int pointer, const int numberOfAtoms) {
    string line;
    istringstream iss;

    for (int i = 1; i <= numberOfAtoms; i++) {
        iss.clear();
        iss.str(allForceFieldLines[pointer + i]);
        ForceFieldAtom forceFieldAtom;

        iss >> forceFieldAtom;
        atoms.push_back(forceFieldAtom);
    }
}

int ForceField::getNumber(const string& line) {
    istringstream iss;
    string title;
    int number;

    iss.clear();
    iss.str(line);
    iss >> title >> number;

    return number;
}

ifstream ForceField::checkFileAndOpen() {
    ifstream forceFieldFile;

    cout << "** Reading force_field.dat **" << endl;
    forceFieldFile.open("force_field.dat");
    if (!forceFieldFile.is_open()) {
        throw runtime_error("** Error: cannot open force_field.dat **");
    }

    return forceFieldFile;
}
