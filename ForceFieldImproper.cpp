#include "ForceFieldImproper.h"

#include <istream>
#include <ostream>
#include <string>

using namespace std;

ForceFieldImproper::ForceFieldImproper(
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
}

std::ostream& operator<<(std::ostream& stream, const ForceFieldImproper& forceFieleImproper) {
    stream << "(" << forceFieleImproper.atom1
           << ", " << forceFieleImproper.atom2
           << ", " << forceFieleImproper.atom3
           << ", " << forceFieleImproper.atom4
           << ", " << forceFieleImproper.functional
           << ", " << forceFieleImproper.forceConstant
           << ", " << forceFieleImproper.eqAngle
           << ", " << forceFieleImproper.multiplicity
           << ", " << forceFieleImproper.oneFourSetter
           << ")";

    return stream;
}
std::istream& operator>>(std::istream& stream, ForceFieldImproper& forceFieleImproper) {
    stream >>
        forceFieleImproper.atom1 >>
        forceFieleImproper.atom2 >>
        forceFieleImproper.atom3 >>
        forceFieleImproper.atom4 >>
        forceFieleImproper.functional >>
        forceFieleImproper.forceConstant >>
        forceFieleImproper.eqAngle >>
        forceFieleImproper.multiplicity >>
        forceFieleImproper.oneFourSetter;

    return stream;
}