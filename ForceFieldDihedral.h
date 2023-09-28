#ifndef FORCEFIELDDIHEDRAL_H
#define FORCEFIELDDIHEDRAL_H

#include <istream>
#include <ostream>
#include <string>

class ForceFieldDihedral {
   public:
    std::string atom1;
    std::string atom2;
    std::string atom3;
    std::string atom4;
    std::string functional;
    double forceConstant;
    double eqAngle;
    short multiplicity;
    short oneFourSetter;

    ForceFieldDihedral() = default;
    ForceFieldDihedral(const std::string& atom1, const std::string& atom2, const std::string& atom3, const std::string& atom4,
                       const std::string& functional, const double forceConstant, const double eqAngle, const short multiplicity,
                       const short oneFourSetter);
};

std::ostream& operator<<(std::ostream& stream, const ForceFieldDihedral& forceFieleDihedral);
std::istream& operator>>(std::istream& stream, ForceFieldDihedral& forceFieleDihedral);

#endif