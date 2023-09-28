#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include <fstream>
#include <string>
#include <vector>

#include "ForceFieldAngle.h"
#include "ForceFieldAtom.h"
#include "ForceFieldBond.h"
#include "ForceFieldDihedral.h"
#include "ForceFieldImproper.h"
#include "ForceFieldVDW.h"

class ForceField {
   public:
    std::vector<ForceFieldAtom> atoms;
    std::vector<ForceFieldBond> bonds;
    std::vector<ForceFieldAngle> angles;
    std::vector<ForceFieldDihedral> dihedrals;
    std::vector<ForceFieldImproper> impropers;
    std::vector<ForceFieldVDW> vdws;

    ForceField();

   private:
    std::ifstream checkFileAndOpen();
    int getNumber(const std::string& line);
    void prepareAtoms(const std::vector<std::string>& allForceFieldLines, const int pointer, const int numberOfAtoms);
    void prepareBonds(const std::vector<std::string>& allForceFieldLines, const int pointer, const int numberOfBonds);
    void prepareAngles(const std::vector<std::string>& allForceFieldLines, const int pointer, const int numberOfAngles);
    void prepareDihedrals(const std::vector<std::string>& allForceFieldLines, const int pointer, const int numberOfDihedrals);
    void prepareImpropers(const std::vector<std::string>& allForceFieldLines, const int pointer, const int numberOfImpropers);
    void prepareVDWs(const std::vector<std::string>& allForceFieldLines, const int pointer, const int numberOfVDWs);
};

#endif