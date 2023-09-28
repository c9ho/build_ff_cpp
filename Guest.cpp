#include "Guest.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

// GuestMolecules::GuestMolecules(const vector<string>& fileNames) {
//     string line;

//     for (auto fileName : fileNames) {
//         auto guestMolecule = GuestMolecule(fileName);
//         for (auto existingMolecule : this->guestMolecules) {
//             if (existingMolecule.moleculeName == guestMolecule.moleculeName) {
//                 line = "** Error: Duplicated type of guest molecule for file " + fileName + " and " + existingMolecule.moleculeName;
//                 throw runtime_error(line);
//             }
//         }
//         this->guestMolecules.push_back(guestMolecule);
//     }
// }

GuestMolecule::GuestMolecule(const string& fileName) {
    analyzeMolecule(fileName);
}

void GuestMolecule::analyzeMolecule(const string& fileName) {
    ifstream file(fileName);
    vector<string> coordinateLines;
    istringstream iss;
    string line;
    string name;

    if (!file.is_open()) {
        line = "** Error: cannot open the file " + fileName + " **";
        throw runtime_error(line);
    }

    getline(file, line);
    iss.clear();
    iss.str(line);
    iss >> name;
    transform(name.begin(), name.end(), name.begin(), ::tolower);

    while (getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        coordinateLines.push_back(line);
    }

    findAndSetMolecule(name, coordinateLines);

    file.close();
}

void GuestMolecule::findAndSetMolecule(const string& name, const vector<string>& allCoordinateLines) {
    int numberOfAtoms = allCoordinateLines.size();
    string atom;
    istringstream iss;
    double x, y, z;

    if (name == "h2o") {
        this->moleculeName = "h2o";
        this->atomsPerMolecule = 3;
        for (int i = 0; i < numberOfAtoms; i++) {
            iss.clear();
            iss.str(allCoordinateLines[i]);
            iss >> atom;

            if (i % 3 == 0 && atom == "O") {
                this->atoms.push_back("OW");
            } else if ((i % 3 == 1 || i % 3 == 2) && atom == "H") {
                this->atoms.push_back("HW");
            } else {
                throw runtime_error("** Error: incorrect order of H2O. It should be O H H **");
            }
        }

    } else if (name == "co2") {
        this->moleculeName = "co2";
        this->atomsPerMolecule = 3;
        for (int i = 0; i < numberOfAtoms; i++) {
            iss.clear();
            iss.str(allCoordinateLines[i]);
            iss >> atom;

            if (i % 3 == 0 && atom == "C") {
                this->atoms.push_back("CC");
            } else if ((i % 3 == 1 || i % 3 == 2) && atom == "O") {
                this->atoms.push_back("OC");
            } else {
                throw runtime_error("** Error: incorrect order of CO2. It should be C O O **");
            }
        }
    } else {
        throw runtime_error("** Error: unknown guest molecule. **");
    }

    if (numberOfAtoms % (this->atomsPerMolecule) == 0) {
        this->numberOfMolecule = numberOfAtoms / (this->atomsPerMolecule);
    } else {
        throw runtime_error("** Error: invalid number of molecules. **");
    }

    for (auto line : allCoordinateLines) {
        iss.clear();
        iss.str(line);
        iss >> atom >> x >> y >> z;

        array coordinate = {x, y, z};
        this->coordinates.push_back(coordinate);
    }
}

ostream& operator<<(ostream& stream, const array<double, 3>& array) {
    stream << "(" << array[0]
           << ", " << array[1]
           << ", " << array[2]
           << ")";

    return stream;
}