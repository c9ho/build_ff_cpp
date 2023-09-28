#include "ForceFieldVDW.h"

#include <istream>
#include <ostream>
#include <string>

using namespace std;

ForceFieldVDW::ForceFieldVDW(const string& atom, const string& type, const double rOver2, const double epsilon) {
    this->atom = atom;
    this->type = type;
    this->rOver2 = rOver2;
    this->epsilon = epsilon;
}

ostream& operator<<(ostream& stream, const ForceFieldVDW& forceFieldVDW) {
    stream << "(" << forceFieldVDW.atom
           << ", " << forceFieldVDW.type
           << ", " << forceFieldVDW.rOver2
           << ", " << forceFieldVDW.epsilon
           << ")";

    return stream;
}

istream& operator>>(istream& stream, ForceFieldVDW& forceFieldVDW) {
    stream >>
        forceFieldVDW.atom >>
        forceFieldVDW.type >>
        forceFieldVDW.rOver2 >>
        forceFieldVDW.epsilon;

    return stream;
}