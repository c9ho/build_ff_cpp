#include "ForceFieldBond.h"

#include <istream>
#include <ostream>
#include <string>

using namespace std;

ForceFieldBond::ForceFieldBond(const string& atom1, const string& atom2, const string& functional, const double forceConstant, const double eqDistance) {
    this->atom1 = atom1;
    this->atom2 = atom2;
    this->functional = functional;
    this->forceConstant = forceConstant;
    this->eqDistance = eqDistance;
};

ostream& operator<<(ostream& stream, const ForceFieldBond& forceFieldBond) {
    stream << "(" << forceFieldBond.atom1
           << ", " << forceFieldBond.atom2
           << ", " << forceFieldBond.functional
           << ", " << forceFieldBond.forceConstant
           << ", " << forceFieldBond.eqDistance
           << ")";
    return stream;
}

istream& operator>>(istream& stream, ForceFieldBond& forceFieldBond) {
    stream >>
        forceFieldBond.atom1 >>
        forceFieldBond.atom2 >>
        forceFieldBond.functional >>
        forceFieldBond.forceConstant >>
        forceFieldBond.eqDistance;

    return stream;
}