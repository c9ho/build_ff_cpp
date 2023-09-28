#ifndef FORCEFIELDBOND_H
#define FORCEFIELDBOND_H

#include <istream>
#include <ostream>
#include <string>

class ForceFieldBond {
   public:
    std::string atom1;
    std::string atom2;
    std::string functional;
    double forceConstant;
    double eqDistance;

    ForceFieldBond() = default;
    ForceFieldBond(const std::string& atom1, const std::string& atom2, const std::string& functional, const double forceConstant, const double eqDistance);
};

std::ostream& operator<<(std::ostream& stream, const ForceFieldBond& forceFieldBond);
std::istream& operator>>(std::istream& stream, ForceFieldBond& forceFieldBond);

#endif