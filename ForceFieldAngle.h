#ifndef FORCEFIELDANGLE_H
#define FORCEFIELDANGLE_H

#include <string>

class ForceFieldAngle {
   public:
    std::string atom1;
    std::string atomC;
    std::string atom2;
    std::string functional;
    double forceConstant;
    double eqAngle;
    bool isRestricted = false;
    double restrictUpperAngle = 0;
    double restrictLowerAngle = 0;

    ForceFieldAngle() = default;
    ForceFieldAngle(const std::string& atom1, const std::string& atomC, const std::string& atom2, const std::string& functional, const double forceConstant, const double eqAngle);
};

std::ostream& operator<<(std::ostream& stream, const ForceFieldAngle& forceFieleAngle);
std::istream& operator>>(std::istream& stream, ForceFieldAngle& forceFieldAngle);

#endif