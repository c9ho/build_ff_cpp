#ifndef FORCEFIELDATOM_H
#define FORCEFIELDATOM_H

#include <istream>
#include <ostream>
#include <string>

class ForceFieldAtom {
   public:
    std::string specificAtomType;
    std::string generalAtomType;
    double charge;
    double mass;

    ForceFieldAtom() = default;
    ForceFieldAtom(const std::string& specificAtomType, const std::string& generalAtomType, const double charge, const double mass);
    bool operator==(const ForceFieldAtom& other) const;
};

// class ForceFieldAtomHash {
//    public:
//     size_t operator()(const ForceFieldAtom& forceFieldAtom) const;
// };

std::ostream& operator<<(std::ostream& stream, const ForceFieldAtom& forceFieldAtom);
std::istream& operator>>(std::istream& stream, ForceFieldAtom& forceFieldAtom);

#endif