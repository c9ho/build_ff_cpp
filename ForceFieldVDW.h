#ifndef FORCEFIELDVDW_H
#define FORCEFIELDVDW_H

#include <istream>
#include <ostream>
#include <string>

class ForceFieldVDW {
   public:
    std::string atom;
    std::string type;
    double rOver2;
    double epsilon;

    ForceFieldVDW() = default;
    ForceFieldVDW(const std::string& atom, const std::string& type, const double rOver2, const double epsilon);
};

std::ostream& operator<<(std::ostream& stream, const ForceFieldVDW& forceFieldVDW);
std::istream& operator>>(std::istream& stream, ForceFieldVDW& forceFieldVDW);

#endif