#include "ForceFieldDihedral.h"

#include <istream>
#include <ostream>
#include <string>

using namespace std;

ForceFieldDihedral::ForceFieldDihedral(
    const string& atom1, const string& atom2, const string& atom3, const string& atom4,
    const string& functional, const double forceConstant, const double eqAngle, const short multiplicity, const short oneFourSetter) {
    this->atom1 = atom1;
    this->atom2 = atom2;
    this->atom3 = atom3;
    this->atom4 = atom4;
    this->functional = functional;
    this->forceConstant = forceConstant;
    this->eqAngle = eqAngle;
    this->multiplicity = multiplicity;
    this->oneFourSetter = oneFourSetter;
};

ostream& operator<<(ostream& stream, const ForceFieldDihedral& forceFieldDihedral) {
    stream << "(" << forceFieldDihedral.atom1
           << ", " << forceFieldDihedral.atom2
           << ", " << forceFieldDihedral.atom3
           << ", " << forceFieldDihedral.atom4
           << ", " << forceFieldDihedral.functional
           << ", " << forceFieldDihedral.forceConstant
           << ", " << forceFieldDihedral.eqAngle
           << ", " << forceFieldDihedral.multiplicity
           << ", " << forceFieldDihedral.oneFourSetter
           << ")";

    return stream;
}

istream& operator>>(istream& stream, ForceFieldDihedral& forceFieldDihedral) {
    stream >>
        forceFieldDihedral.atom1 >>
        forceFieldDihedral.atom2 >>
        forceFieldDihedral.atom3 >>
        forceFieldDihedral.atom4 >>
        forceFieldDihedral.functional >>
        forceFieldDihedral.forceConstant >>
        forceFieldDihedral.eqAngle >>
        forceFieldDihedral.multiplicity >>
        forceFieldDihedral.oneFourSetter;

    return stream;
}