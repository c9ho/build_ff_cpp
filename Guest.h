#ifndef GUEST_H
#define GUEST_H

#include <array>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

class GuestMolecule {
   public:
    std::string moleculeName;
    int atomsPerMolecule;
    std::vector<std::string> atoms;
    std::vector<std::array<double, 3>> coordinates;
    int numberOfMolecule;

    GuestMolecule(const std::string& fileName);

   private:
    void analyzeMolecule(const std::string& fileName);
    void findAndSetMolecule(const std::string& name, const std::vector<std::string>& coordinateLines);
};

// class GuestMolecules {
//    public:
//     std::vector<GuestMolecule> guestMolecules;

//     GuestMolecules(const std::vector<std::string>& fileNames);
// };

std::ostream& operator<<(std::ostream& stream, const std::array<double, 3>& array);

#endif