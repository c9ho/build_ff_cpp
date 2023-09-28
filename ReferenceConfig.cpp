#include "ReferenceConfig.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "tool.h"

using namespace std;

ReferenceAtom::ReferenceAtom(string& atomName, double xyzCoordinate[3], double abcCoordinate[3], int index) {
    this->atomName = atomName;
    for (int i = 0; i < 3; i++) {
        this->xyzCoordinate[i] = xyzCoordinate[i];
        this->abcCoordinate[i] = abcCoordinate[i];
        // cout << abcCoordinate[i] << endl;
        // cout << xyzCoordinate[i] << endl;
    }
    this->atomIndex = index;
}

void ReferenceConfig::getZeroToOne(double numbers[3]) {
    for (int i = 0; i < 3; i++)
        numbers[i] = numbers[i] - floor(numbers[i]);
}

double ReferenceConfig::calcDistance(const double abcCoor1[3], const double abcCoor2[3]) const {
    auto [vecA, vecB, vecC] = this->box;
    double diff[3] = {0.0};
    double diffX, diffY, diffZ = {0.0};

    for (int i = 0; i < 3; i++)
        diff[i] = abcCoor2[i] - abcCoor1[i];

    // cout << diff[0] << " " << diff[1] << " " << diff[2];
    getFirstBZone(diff);

    diffX = diff[0] * vecA[0] + diff[1] * vecB[0] + diff[2] * vecC[0];
    diffY = diff[0] * vecA[1] + diff[1] * vecB[1] + diff[2] * vecC[1];
    diffZ = diff[0] * vecA[2] + diff[1] * vecB[2] + diff[2] * vecC[2];

    return sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);
}

double ReferenceConfig::calcAngle(const double abcCoor1[3], const double abcCoor2[3], const double abcCoor3[3]) const {
    auto [vecA, vecB, vecC] = this->box;
    double vec_21[3];
    double vec_23[3];
    double vec_21_xyz[3];
    double vec_23_xyz[3];
    double vec_21_norm;
    double vec_23_norm;

    for (int i = 0; i < 3; i++) {
        vec_21[i] = abcCoor2[i] - abcCoor1[i];
        vec_23[i] = abcCoor2[i] - abcCoor3[i];
    }
    getFirstBZone(vec_21);
    getFirstBZone(vec_23);
    for (int i = 0; i < 3; i++) {
        vec_21_xyz[i] = vec_21[0] * vecA[i] + vec_21[1] * vecB[i] + vec_21[2] * vecC[i];
        vec_23_xyz[i] = vec_23[0] * vecA[i] + vec_23[1] * vecB[i] + vec_23[2] * vecC[i];
        vec_21_norm += vec_21_xyz[i] * vec_21_xyz[i];
        vec_23_norm += vec_23_xyz[i] * vec_23_xyz[i];
    }
    vec_21_norm = sqrt(vec_21_norm);
    vec_23_norm = sqrt(vec_23_norm);

    double cosAngle = innerProduct(vec_21_xyz, vec_23_xyz) / vec_21_norm / vec_23_norm;

    return acos(cosAngle) * 180 / PI;
}

ReferenceConfig::ReferenceConfig() {
    fstream referenceFile("reference.xyz");

    cout << "** Reading reference.xyz **" << endl;
    if (referenceFile.is_open()) {
        string line;
        istringstream iss;
        int lineCounter = 0;

        while (getline(referenceFile, line)) {
            if (line.empty()) {
                break;
            }

            lineCounter++;
            iss.clear();
            iss.str(line);

            if (lineCounter == 1) {
                prepareNumberOfAtoms(iss);
            } else if (lineCounter == 2) {
                prepareBox(iss);
                prepareReciprocalVectors();
            } else {
                int atomIndex = lineCounter - 2;
                prepareAtoms(iss, atomIndex);
            }
        }
        referenceFile.close();

    } else {
        throw runtime_error("** Error: cannot open reference.xyz **");
    }
}

void ReferenceConfig::prepareReciprocalVectors() {
    auto [vecA, vecB, vecC] = this->box;
    double bCrossC[3];

    outerProduct(vecB, vecC, bCrossC);
    double volume = innerProduct(vecA, bCrossC);
    // cout << volume << endl;

    outerProduct(vecB, vecC, this->reciprocalVectors[0]);
    outerProduct(vecC, vecA, this->reciprocalVectors[1]);
    outerProduct(vecA, vecB, this->reciprocalVectors[2]);
    for (int i = 0; i < 3; i++) {
        reciprocalVectors[0][i] /= volume;
        reciprocalVectors[1][i] /= volume;
        reciprocalVectors[2][i] /= volume;
    }
}

void ReferenceConfig::prepareNumberOfAtoms(istringstream& iss) {
    iss >> this->numberOfAtoms;
}

void ReferenceConfig::prepareBox(istringstream& iss) {
    iss >> this->box[0][0] >> this->box[0][1] >> this->box[0][2] >>
        this->box[1][0] >> this->box[1][1] >> this->box[1][2] >>
        this->box[2][0] >> this->box[2][1] >> this->box[2][2];
}

void ReferenceConfig::prepareAtoms(istringstream& iss, const int atomIndex) {
    string atomName;
    double x, y, z;

    iss >> atomName >> x >> y >> z;
    double xyzCoors[3] = {x, y, z};

    double a = innerProduct(xyzCoors, this->reciprocalVectors[0]);
    double b = innerProduct(xyzCoors, this->reciprocalVectors[1]);
    double c = innerProduct(xyzCoors, this->reciprocalVectors[2]);
    double abcCoors[3] = {a, b, c};
    getZeroToOne(abcCoors);

    ReferenceAtom referenceAtom(atomName, xyzCoors, abcCoors, atomIndex);
    referenceAtoms.push_back(referenceAtom);
}

void ReferenceConfig::getFirstBZone(double numbers[3]) {
    for (int i = 0; i < 3; i++) {
        if (numbers[i] > 0.5)
            numbers[i] -= 1.0;
        else if (numbers[i] < -0.5)
            numbers[i] += 1.0;
    }
}