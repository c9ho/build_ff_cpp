#include "ForceFieldAngle.h"

#include <istream>
#include <ostream>
#include <string>

using namespace std;

ForceFieldAngle::ForceFieldAngle(const string& atom1, const string& atomC, const string& atom2, const string& functional, const double forceConstant, const double eqAngle) {
    this->atom1 = atom1;
    this->atomC = atomC;
    this->atom2 = atom2;
    this->functional = functional;
    this->forceConstant = forceConstant;
    this->eqAngle = eqAngle;
};

ostream& operator<<(ostream& stream, const ForceFieldAngle& forceFieldAngle) {
    stream << "(" << forceFieldAngle.atom1
           << ", " << forceFieldAngle.atomC
           << ", " << forceFieldAngle.atom2
           << ", " << forceFieldAngle.functional
           << ", " << forceFieldAngle.forceConstant
           << ", " << forceFieldAngle.eqAngle
           << ")";

    return stream;
}

istream& operator>>(istream& stream, ForceFieldAngle& forceFieldAngle) {
    stream >>
        forceFieldAngle.atom1 >>
        forceFieldAngle.atomC >>
        forceFieldAngle.atom2 >>
        forceFieldAngle.functional >>
        forceFieldAngle.forceConstant >>
        forceFieldAngle.eqAngle;

    return stream;
}