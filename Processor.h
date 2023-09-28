#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <unordered_map>

#include "ForceField.h"
#include "ForceFieldAtom.h"
#include "Guest.h"
#include "ReferenceConfig.h"

class Processor {
   public:
    ForceField forceField;
    ReferenceConfig referenceConfig;
    std::vector<GuestMolecule> guestMolecules;
    std::vector<std::array<std::string, 2>> metalPairs;

    std::string atomLine;
    std::string massLine;
    std::string bondCountLine;
    std::string bondIndexLine;
    std::string angleCountLine;
    std::string angleIndexLine;
    std::string dihedralCountLine;
    std::string dihedralIndexLine;
    std::string improperCountLine;
    std::string improperIndexLine;
    std::string aLine;
    std::string bLine;
    std::string cLine;
    std::string tiltLine;

    Processor(const std::vector<std::string>& fileNames);
    void processMassPart();
    void processAtomPart();
    void processBondPart();
    void processAnglePart();
    void processDihedralPart();
    void processImproperPart();
    void processVDWPart();
    void processBoxInfo();

    std::vector<ForceFieldAtom> processedForceFieldAtoms;
    std::unordered_map<std::string, int> generalToIndex;
    std::unordered_map<int, std::string> indexToGeneral;
    std::unordered_map<std::string, std::string> specToGeneral;
    std::unordered_map<std::string, double> specToCharge;
    std::unordered_map<std::string, double> generalToMass;

   private:
    void findGuestMolecules(const std::vector<std::string>& fileNames);
    void findProcessedForceFieldAtoms();
    void findSpecificCorrespondance();
    void prepareMetalPairs();
    void changeCutoffIfNeeded(const std::string& atom1, const std::string& atom2, double& cutoff);
};

#endif