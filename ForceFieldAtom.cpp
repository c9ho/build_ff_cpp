#include "ForceFieldAtom.h"

#include <istream>
#include <ostream>
#include <string>

using namespace std;

ForceFieldAtom::ForceFieldAtom(const string& specificAtomType, const string& generalAtomType, const double charge, const double mass) {
    this->specificAtomType = specificAtomType;
    this->generalAtomType = generalAtomType;
    this->charge = charge;
    this->mass = mass;
}

ostream& operator<<(ostream& stream, const ForceFieldAtom& forceFieldAtom) {
    stream << "(" << forceFieldAtom.specificAtomType
           << ", " << forceFieldAtom.generalAtomType
           << ", " << forceFieldAtom.charge
           << ", " << forceFieldAtom.mass
           << ")";
    return stream;
}
istream& operator>>(istream& stream, ForceFieldAtom& forceFieldAtom) {
    stream >>
        forceFieldAtom.specificAtomType >>
        forceFieldAtom.generalAtomType >>
        forceFieldAtom.charge >>
        forceFieldAtom.mass;

    return stream;
}

bool ForceFieldAtom::operator==(const ForceFieldAtom& other) const {
    return ((specificAtomType == other.specificAtomType) && (generalAtomType == other.generalAtomType) && (charge == other.charge) && (mass == other.mass));
}

// size_t ForceFieldAtomHash::operator()(const ForceFieldAtom& forceFieldAtom) const {
//     string specific = forceFieldAtom.specificAtomType;
//     string general = forceFieldAtom.generalAtomType;
//     double charge = forceFieldAtom.charge;
//     double mass = forceFieldAtom.mass;

//     string line = specific + general + to_string(charge) + to_string(mass);

//     return hash<string>()(line);
// }