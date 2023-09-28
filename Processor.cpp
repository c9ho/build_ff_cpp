#include "Processor.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_set>

#include "Guest.h"
#include "MetalPairsReader.h"

using namespace std;

Processor::Processor(const vector<string>& fileNames) {
    findGuestMolecules(fileNames);
    prepareMetalPairs();
    findProcessedForceFieldAtoms();
    findSpecificCorrespondance();
}

void Processor::processBoxInfo() {
    auto [vecA, vecB, vecC] = referenceConfig.box;
    auto [xhi, a2, a3] = vecA;
    auto [xy, yhi, b3] = vecB;
    auto [xz, yz, zhi] = vecC;

    stringstream sstreamA, sstreamB, sstreamC, sstreamTilt;
    sstreamA << setw(19) << "0.000000" << setw(12) << fixed << setprecision(7) << xhi << " xlo  xhi";
    sstreamB << setw(19) << "0.000000" << setw(12) << fixed << setprecision(7) << yhi << " ylo  yhi";
    sstreamC << setw(19) << "0.000000" << setw(12) << fixed << setprecision(7) << zhi << " zlo  zhi";
    sstreamTilt << setw(19) << fixed << setprecision(7) << xy << setw(12) << xz << setw(12) << yz << " xy  xz  yz";
    aLine = sstreamA.str();
    bLine = sstreamB.str();
    cLine = sstreamC.str();
    tiltLine = sstreamTilt.str();

    if (a2 != 0.0 || a3 != 0.0 || b3 != 0.0)
        cout << "Alert! The lattice vectors do not follow the definitions in LAMMPS. You need to manually fix it." << endl;
}

void Processor::processVDWPart() {
    const double twoToTheOneOverSix = 1.122462048;

    const double OWROver2 = 1.77287268;
    const double OWEpsilon = 0.18519;
    const double HWROver2 = 0.0;
    const double HWEpsilon = 0.0;
    const double CCROver2 = 1.571446867;
    const double CCEpsilon = 0.05363150574;
    const double OCROver2 = 1.711754623;
    const double OCEpsilon = 0.1569218131;

    ofstream pairFile("PAIR_COEFFICIENTS");

    cout << "Starting vdw part." << endl;
    pairFile << "# MOF pair coefficients\n";
    for (auto& item1 : forceField.vdws) {
        int atom1Index = generalToIndex[item1.atom];
        for (auto& item2 : forceField.vdws) {
            int atom2Index = generalToIndex[item2.atom];
            if (atom1Index > atom2Index)
                continue;
            double pairSigma = (item1.rOver2 + item2.rOver2) / twoToTheOneOverSix;
            double pairEpsilon = sqrt(item1.epsilon * item2.epsilon);
            if (item1.rOver2 == 0.0 || item2.rOver2 == 0.0)
                pairSigma = 0.0;
            pairFile << right
                     << "pair_coeff"
                     << setw(5)
                     << atom1Index
                     << setw(10)
                     << atom2Index
                     << setw(10)
                     << "lj/cut"
                     << setw(15)
                     << fixed
                     << setprecision(6)
                     << pairEpsilon
                     << setw(15)
                     << pairSigma
                     << "\n";
        }
    }
    for (auto& guestMolecule : guestMolecules) {
        if (guestMolecule.moleculeName == "h2o")
            pairFile << "\n"
                     << "# H2O MOF Mixed\n";
        else if (guestMolecule.moleculeName == "co2")
            pairFile << "\n"
                     << "# CO2 MOF Mixed\n";
        for (auto& item : forceField.vdws) {
            int atomIndex = generalToIndex[item.atom];
            if (guestMolecule.moleculeName == "h2o") {
                int OWIndex = generalToIndex["OW"];
                int HWIndex = generalToIndex["HW"];

                double pairOWsigma = (item.rOver2 == 0.0) ? 0.0 : (item.rOver2 + OWROver2) / twoToTheOneOverSix;
                double pairHWsigma = (item.rOver2 == 0.0 || HWROver2 == 0.0) ? 0.0 : (item.rOver2 + HWROver2) / twoToTheOneOverSix;
                double pairOWEpsilon = sqrt(item.epsilon * OWEpsilon);
                double pairHWEpsilon = sqrt(item.epsilon * HWEpsilon);

                pairFile << right
                         << "pair_coeff"
                         << setw(5)
                         << atomIndex
                         << setw(10)
                         << OWIndex
                         << setw(10)
                         << "lj/cut"
                         << setw(15)
                         << fixed
                         << setprecision(6)
                         << pairOWEpsilon
                         << setw(15)
                         << pairOWsigma
                         << "\n";
                pairFile << right
                         << "pair_coeff"
                         << setw(5)
                         << atomIndex
                         << setw(10)
                         << HWIndex
                         << setw(10)
                         << "lj/cut"
                         << setw(15)
                         << fixed
                         << setprecision(6)
                         << pairHWEpsilon
                         << setw(15)
                         << pairHWsigma
                         << "\n";
            } else if (guestMolecule.moleculeName == "co2") {
                int CCIndex = generalToIndex["CC"];
                int OCIndex = generalToIndex["OC"];

                double pairCCsigma = (item.rOver2 == 0.0) ? 0.0 : (item.rOver2 + CCROver2) / twoToTheOneOverSix;
                double pairOCsigma = (item.rOver2 == 0.0) ? 0.0 : (item.rOver2 + OCROver2) / twoToTheOneOverSix;
                double pairCCEpsilon = sqrt(item.epsilon * CCEpsilon);
                double pairOCEpsilon = sqrt(item.epsilon * OCEpsilon);

                pairFile << right
                         << "pair_coeff"
                         << setw(5)
                         << atomIndex
                         << setw(10)
                         << CCIndex
                         << setw(10)
                         << "lj/cut"
                         << setw(15)
                         << fixed
                         << setprecision(6)
                         << pairCCEpsilon
                         << setw(15)
                         << pairCCsigma
                         << "\n";
                pairFile << right
                         << "pair_coeff"
                         << setw(5)
                         << atomIndex
                         << setw(10)
                         << OCIndex
                         << setw(10)
                         << "lj/cut"
                         << setw(15)
                         << fixed
                         << setprecision(6)
                         << pairOCEpsilon
                         << setw(15)
                         << pairOCsigma
                         << "\n";
            }
        }
    }

    pairFile.close();
}

void Processor::processImproperPart() {
    ofstream improperFile("IMPROPER");
    ofstream improperCofficientFile("improper_coeff");

    int improperCount = 0;
    int improperIndex = 0;

    cout << "Starting improper part." << endl;
    improperFile << "Impropers\n"
                 << "\n";

    for (auto& item : forceField.impropers) {
        vector<ReferenceAtom> atom1;
        vector<ReferenceAtom> atom2;
        vector<ReferenceAtom> atom3;
        vector<ReferenceAtom> atom4;
        double cutoff12 = 1.7;
        double cutoff13 = 1.7;
        double cutoff14 = 1.7;
        changeCutoffIfNeeded(item.atom1, item.atom2, cutoff12);
        changeCutoffIfNeeded(item.atom2, item.atom3, cutoff13);
        changeCutoffIfNeeded(item.atom3, item.atom4, cutoff14);
        improperIndex += 1;
        int improperCountBefore = improperCount;

        for (auto& atom : referenceConfig.referenceAtoms) {
            if (specToGeneral[atom.atomName] == item.atom1)
                atom1.push_back(atom);
            if (specToGeneral[atom.atomName] == item.atom2)
                atom2.push_back(atom);
            if (specToGeneral[atom.atomName] == item.atom3)
                atom3.push_back(atom);
            if (specToGeneral[atom.atomName] == item.atom4)
                atom4.push_back(atom);
        }
        for (auto& item1 : atom1) {
            for (auto& item2 : atom2) {
                if (item1.atomIndex == item2.atomIndex)
                    continue;
                double distance12 = referenceConfig.calcDistance(item1.abcCoordinate, item2.abcCoordinate);
                if (distance12 < cutoff12) {
                    for (auto& item3 : atom3) {
                        if (item3.atomIndex == item1.atomIndex || item3.atomIndex == item2.atomIndex)
                            continue;
                        if (item.atom2 == item.atom3 && item2.atomIndex >= item3.atomIndex)
                            continue;
                        double distance13 = referenceConfig.calcDistance(item1.abcCoordinate, item3.abcCoordinate);
                        if (distance13 < cutoff13) {
                            for (auto& item4 : atom4) {
                                if (item4.atomIndex == item1.atomIndex || item4.atomIndex == item2.atomIndex || item4.atomIndex == item3.atomIndex)
                                    continue;
                                if (item.atom3 == item.atom4 && item3.atomIndex >= item4.atomIndex)
                                    continue;
                                if (item.atom2 == item.atom4 && item2.atomIndex >= item4.atomIndex)
                                    continue;
                                double distance14 = referenceConfig.calcDistance(item1.abcCoordinate, item4.abcCoordinate);
                                if (distance14 < cutoff14) {
                                    improperCount += 1;
                                    improperFile << right
                                                 << improperCount
                                                 << setw(7)
                                                 << improperIndex
                                                 << setw(7)
                                                 << item2.atomIndex
                                                 << setw(7)
                                                 << item1.atomIndex
                                                 << setw(7)
                                                 << item3.atomIndex
                                                 << setw(7)
                                                 << item4.atomIndex
                                                 << "\n";
                                }
                            }
                        }
                    }
                }
            }
        }
        if (improperCount != improperCountBefore) {
            int dValue = -1;
            if (static_cast<int>(item.eqAngle) == 0)
                dValue = 1;

            improperCofficientFile << right
                                   << "improper_coeff"
                                   << setw(10)
                                   << improperIndex
                                   << setw(15)
                                   << fixed
                                   << setprecision(5)
                                   << item.forceConstant
                                   << setw(15)
                                   << dValue
                                   << setw(10)
                                   << item.multiplicity
                                   << " # "
                                   << item.atom2
                                   << "-"
                                   << item.atom1
                                   << "-"
                                   << item.atom3
                                   << "-"
                                   << item.atom4
                                   << "\n";
        } else {
            improperIndex -= 1;
        }
    }

    stringstream numberStream1, numberStream2;
    numberStream1 << setw(12) << improperCount;
    numberStream2 << setw(12) << improperIndex;
    improperCountLine = numberStream1.str() + " impropers";
    improperIndexLine = numberStream2.str() + " improper types";

    improperFile << "\n";

    improperFile.close();
    improperCofficientFile.close();
}

void Processor::processDihedralPart() {
    ofstream dihedralFile("DIHEDRAL");
    ofstream dihedralCoefficientFile("dihedral_coeff");

    int dihedralCount = 0;
    int dihedralIndex = 0;

    cout << "Starting dihedrals part." << endl;
    dihedralFile << "Dihedrals\n"
                 << "\n";

    for (auto& item : forceField.dihedrals) {
        vector<ReferenceAtom> atom1;
        vector<ReferenceAtom> atom2;
        vector<ReferenceAtom> atom3;
        vector<ReferenceAtom> atom4;
        double cutoff12 = 1.7;
        double cutoff23 = 1.7;
        double cutoff34 = 1.7;
        changeCutoffIfNeeded(item.atom1, item.atom2, cutoff12);
        changeCutoffIfNeeded(item.atom2, item.atom3, cutoff23);
        changeCutoffIfNeeded(item.atom3, item.atom4, cutoff34);
        dihedralIndex += 1;
        int dihedralCountBefore = dihedralCount;

        for (auto& atom : referenceConfig.referenceAtoms) {
            if (specToGeneral[atom.atomName] == item.atom1)
                atom1.push_back(atom);
            if (specToGeneral[atom.atomName] == item.atom2)
                atom2.push_back(atom);
            if (specToGeneral[atom.atomName] == item.atom3)
                atom3.push_back(atom);
            if (specToGeneral[atom.atomName] == item.atom4)
                atom4.push_back(atom);
        }

        for (auto& item1 : atom1) {
            for (auto& item2 : atom2) {
                if (item1.atomIndex == item2.atomIndex)
                    continue;
                double distance12 = referenceConfig.calcDistance(item1.abcCoordinate, item2.abcCoordinate);
                if (distance12 < cutoff12) {
                    for (auto& item3 : atom3) {
                        if (item3.atomIndex == item1.atomIndex || item3.atomIndex == item2.atomIndex)
                            continue;
                        double distance23 = referenceConfig.calcDistance(item2.abcCoordinate, item3.abcCoordinate);
                        if (distance23 < cutoff23) {
                            for (auto& item4 : atom4) {
                                if (item4.atomIndex == item1.atomIndex || item4.atomIndex == item2.atomIndex || item4.atomIndex == item3.atomIndex)
                                    continue;
                                if (item.atom1 == item.atom4 && item.atom2 == item.atom3 && item1.atomIndex >= item4.atomIndex)
                                    continue;
                                double distance34 = referenceConfig.calcDistance(item3.abcCoordinate, item4.abcCoordinate);
                                if (distance34 < cutoff34) {
                                    dihedralCount += 1;
                                    dihedralFile << right
                                                 << dihedralCount
                                                 << setw(7)
                                                 << dihedralIndex
                                                 << setw(7)
                                                 << item1.atomIndex
                                                 << setw(7)
                                                 << item2.atomIndex
                                                 << setw(7)
                                                 << item3.atomIndex
                                                 << setw(7)
                                                 << item4.atomIndex
                                                 << "\n";
                                }
                            }
                        }
                    }
                }
            }
        }
        if (dihedralCount != dihedralCountBefore) {
            int dValue = -1;
            if (static_cast<int>(item.eqAngle) == 0)
                dValue = 1;

            dihedralCoefficientFile << right
                                    << "dihedral_coeff"
                                    << setw(10)
                                    << dihedralIndex
                                    << setw(15)
                                    << fixed
                                    << setprecision(5)
                                    << item.forceConstant
                                    << setw(15)
                                    << dValue
                                    << setw(10)
                                    << item.multiplicity
                                    << " # "
                                    << item.atom1
                                    << "-"
                                    << item.atom2
                                    << "-"
                                    << item.atom3
                                    << "-"
                                    << item.atom4
                                    << "\n";
        } else {
            dihedralIndex -= 1;
        }
    }

    stringstream numberStream1, numberStream2;
    numberStream1 << setw(12) << dihedralCount;
    numberStream2 << setw(12) << dihedralIndex;
    dihedralCountLine = numberStream1.str() + " dihedrals";
    dihedralIndexLine = numberStream2.str() + " dihedral types";

    dihedralFile << "\n";

    dihedralFile.close();
    dihedralCoefficientFile.close();
}

void Processor::processAnglePart() {
    ofstream angleFile("ANGLE");
    ofstream angleCoefficientFile("angle_coeff");

    int angleCount = 0;
    int angleIndex = 0;

    int h2oAngles = 0;
    int co2Angles = 0;

    cout << "Starting angles part." << endl;
    angleFile << "Angles\n"
              << "\n";

    for (auto& item : forceField.angles) {
        vector<ReferenceAtom> atom1;
        vector<ReferenceAtom> atomC;
        vector<ReferenceAtom> atom2;
        double cutoff1c = 1.7;
        double cutoff2c = 1.7;
        changeCutoffIfNeeded(item.atom1, item.atomC, cutoff1c);
        changeCutoffIfNeeded(item.atom2, item.atomC, cutoff2c);
        angleIndex += 1;
        int angleCountBefore = angleCount;

        for (auto& atom : referenceConfig.referenceAtoms) {
            if (specToGeneral[atom.atomName] == item.atom1)
                atom1.push_back(atom);
            if (specToGeneral[atom.atomName] == item.atomC)
                atomC.push_back(atom);
            if (specToGeneral[atom.atomName] == item.atom2)
                atom2.push_back(atom);
        }
        for (auto& item1 : atom1) {
            for (auto& itemC : atomC) {
                double distance1c = referenceConfig.calcDistance(item1.abcCoordinate, itemC.abcCoordinate);
                if (distance1c >= cutoff1c)
                    continue;
                for (auto& item2 : atom2) {
                    if (item1.atomIndex == itemC.atomIndex || item2.atomIndex == itemC.atomIndex || item1.atomIndex == item2.atomIndex)
                        continue;
                    else if (item.atom1 == item.atom2 && item1.atomIndex >= item2.atomIndex)
                        continue;
                    double distance2c = referenceConfig.calcDistance(item2.abcCoordinate, itemC.abcCoordinate);
                    if (distance2c < cutoff2c) {
                        if (item.isRestricted) {
                            double angle = referenceConfig.calcAngle(item1.abcCoordinate, itemC.abcCoordinate, item2.abcCoordinate);
                            if (angle < item.restrictLowerAngle || angle > item.restrictUpperAngle)
                                continue;
                        }
                        angleCount += 1;
                        angleFile << right
                                  << angleCount
                                  << setw(7)
                                  << angleIndex
                                  << setw(7)
                                  << item1.atomIndex
                                  << setw(7)
                                  << itemC.atomIndex
                                  << setw(7)
                                  << item2.atomIndex
                                  << "\n";
                    }
                }
            }
        }
        if (angleCount != angleCountBefore) {
            angleCoefficientFile << right
                                 << "angle_coeff"
                                 << setw(5)
                                 << angleIndex
                                 << setw(16)
                                 << fixed
                                 << setprecision(5)
                                 << item.forceConstant / 2
                                 << setw(16)
                                 << item.eqAngle
                                 << " # "
                                 << item.atom1
                                 << "-"
                                 << item.atomC
                                 << "-"
                                 << item.atom2
                                 << "\n";
        } else {
            angleIndex -= 1;
        }
    }

    int numberOfAngleMOF = angleCount;
    int MOFIndex = angleIndex;
    int totalAtoms = referenceConfig.numberOfAtoms;

    for (auto& molecule : guestMolecules) {
        if (molecule.moleculeName == "h2o") {
            angleIndex += 1;
            int h2oInitial = angleCount;
            for (int i = 0; i < molecule.numberOfMolecule; i++) {
                totalAtoms += 1;
                int OWNumber = totalAtoms;
                totalAtoms += 1;
                int HW1Number = totalAtoms;
                totalAtoms += 1;
                int HW2Number = totalAtoms;
                angleCount += 1;
                angleFile << right
                          << angleCount
                          << setw(7)
                          << angleIndex
                          << setw(7)
                          << HW1Number
                          << setw(7)
                          << OWNumber
                          << setw(7)
                          << HW2Number
                          << "\n";
            }
            int h2oFinal = angleCount;
            h2oAngles = h2oFinal - h2oInitial;
        } else if (molecule.moleculeName == "co2") {
            angleIndex += 1;
            int co2Initial = angleCount;
            for (int i = 0; i < molecule.numberOfMolecule; i++) {
                totalAtoms += 1;
                int CCNumber = totalAtoms;
                totalAtoms += 1;
                int OC1Number = totalAtoms;
                totalAtoms += 1;
                int OC2Number = totalAtoms;
                angleCount += 1;
                angleFile << right
                          << angleCount
                          << setw(7)
                          << angleIndex
                          << setw(7)
                          << OC1Number
                          << setw(7)
                          << CCNumber
                          << setw(7)
                          << OC2Number
                          << "\n";
            }
            int co2Final = angleCount;
            co2Angles = co2Final - co2Initial;
        }
    }

    stringstream numberStream1, numberStream2;
    numberStream1 << setw(12) << angleCount;
    numberStream2 << setw(12) << angleIndex;
    angleCountLine = numberStream1.str() + " angles " + " # MOF " + to_string(numberOfAngleMOF);
    angleIndexLine = numberStream2.str() + " angle types " + " # MOF " + to_string(MOFIndex);

    for (auto& molecule : guestMolecules) {
        if (molecule.moleculeName == "h2o") {
            angleCountLine += " water " + to_string(h2oAngles);
            angleIndexLine += " water 1";
        } else if (molecule.moleculeName == "co2") {
            angleCountLine += " CO2 " + to_string(co2Angles);
            angleIndexLine += " CO2 1";
        }
    }
    angleFile << "\n";

    angleFile.close();
    angleCoefficientFile.close();
}

void Processor::processBondPart() {
    ofstream bondFile("BOND");
    ofstream bondCoefficientFile("bond_coeff");

    int bondCount = 0;
    int bondIndex = 0;

    int h2oBonds = 0;
    int co2Bonds = 0;

    cout << "Starting bonds part." << endl;
    bondFile << "Bonds\n"
             << "\n";

    for (auto& item : forceField.bonds) {
        vector<ReferenceAtom> atom1;
        vector<ReferenceAtom> atom2;
        double cutoff = 1.7;
        changeCutoffIfNeeded(item.atom1, item.atom2, cutoff);
        bondIndex += 1;
        int bondCountBefore = bondCount;

        for (auto& atom : referenceConfig.referenceAtoms) {
            if (specToGeneral[atom.atomName] == item.atom1) {
                atom1.push_back(atom);
            }
            if (specToGeneral[atom.atomName] == item.atom2) {
                atom2.push_back(atom);
            }
        }

        for (auto& item1 : atom1) {
            for (auto& item2 : atom2) {
                string atom1Name = specToGeneral[item1.atomName];
                string atom2Name = specToGeneral[item2.atomName];
                if (atom1Name == atom2Name && item1.atomIndex >= item2.atomIndex)
                    continue;
                double distance = referenceConfig.calcDistance(item1.abcCoordinate, item2.abcCoordinate);
                if (distance < cutoff) {
                    bondCount += 1;
                    bondFile << right
                             << bondCount
                             << setw(7)
                             << bondIndex
                             << setw(7)
                             << item1.atomIndex
                             << setw(7)
                             << item2.atomIndex
                             << "\n";
                }
            }
        }
        if (bondCountBefore != bondCount) {
            bondCoefficientFile << right
                                << "bond_coeff"
                                << setw(5)
                                << bondIndex
                                << setw(16)
                                << fixed
                                << setprecision(5)
                                << item.forceConstant / 2
                                << setw(16)
                                << item.eqDistance
                                << setw(5)
                                << "# "
                                << item.atom1
                                << setw(5)
                                << item.atom2
                                << "\n";
        } else {
            bondIndex -= 1;
        }
    }

    int numberOfBondMOF = bondCount;
    int MOFIndex = bondIndex;
    int totalAtoms = referenceConfig.numberOfAtoms;

    for (auto& molecule : guestMolecules) {
        if (molecule.moleculeName == "h2o") {
            int h2oInitial = bondCount;
            bondIndex += 1;
            for (int i = 0; i < molecule.numberOfMolecule; i++) {
                totalAtoms += 1;
                int OWNumber = totalAtoms;
                totalAtoms += 1;
                int HW1Number = totalAtoms;
                totalAtoms += 1;
                int HW2Number = totalAtoms;
                bondCount += 1;
                bondFile << right
                         << bondCount
                         << setw(5)
                         << bondIndex
                         << setw(5)
                         << OWNumber
                         << setw(5)
                         << HW1Number
                         << "\n";
                bondCount += 1;
                bondFile << right
                         << bondCount
                         << setw(5)
                         << bondIndex
                         << setw(5)
                         << OWNumber
                         << setw(5)
                         << HW2Number
                         << "\n";
            }
            int h2oFinal = bondCount;
            h2oBonds = h2oFinal - h2oInitial;
        } else if (molecule.moleculeName == "co2") {
            int co2Initial = bondCount;
            bondIndex += 1;
            for (int i = 0; i < molecule.numberOfMolecule; i++) {
                totalAtoms += 1;
                int CCNumber = totalAtoms;
                totalAtoms += 1;
                int OC1Number = totalAtoms;
                totalAtoms += 1;
                int OC2Number = totalAtoms;
                bondCount += 1;
                bondFile << right
                         << bondCount
                         << setw(5)
                         << bondIndex
                         << setw(5)
                         << CCNumber
                         << setw(5)
                         << OC1Number
                         << "\n";
                bondCount += 1;
                bondFile << right
                         << bondCount
                         << setw(5)
                         << bondIndex
                         << setw(5)
                         << CCNumber
                         << setw(5)
                         << OC2Number
                         << "\n";
            }
            int co2Final = bondCount;
            co2Bonds = co2Final - co2Initial;
        }
    }
    stringstream numberString1, numberString2;
    numberString1 << setw(12) << bondCount;
    numberString2 << setw(12) << bondIndex;
    bondCountLine = numberString1.str() + " bonds " + " # MOF " + to_string(numberOfBondMOF);
    bondIndexLine = numberString2.str() + " bond types " + " # MOF " + to_string(MOFIndex);
    for (auto& molecule : guestMolecules) {
        if (molecule.moleculeName == "h2o") {
            bondCountLine += " water " + to_string(h2oBonds);
            bondIndexLine += " water 1";
        } else if (molecule.moleculeName == "co2") {
            bondCountLine += " co2 " + to_string(co2Bonds);
            bondIndexLine += " CO2 1";
        }
    }

    bondFile << "\n";

    bondFile.close();
    bondCoefficientFile.close();
}

void Processor::processAtomPart() {
    ofstream outputFile("ATOM");
    double totalCharge = 0.0;
    int belongingIndex = 1;
    int indexNumber = 0;

    cout << "Starting Atom part" << endl;
    outputFile << "Atoms\n"
               << "\n";
    for (auto& item : referenceConfig.referenceAtoms) {
        indexNumber++;
        outputFile << right
                   << setw(5)
                   << item.atomIndex
                   << setw(5)
                   << belongingIndex
                   << setw(5)
                   << generalToIndex[specToGeneral[item.atomName]]
                   << setw(15)
                   << fixed
                   << setprecision(7)
                   << specToCharge[item.atomName]
                   << setw(15)
                   << item.xyzCoordinate[0]
                   << setw(15)
                   << item.xyzCoordinate[1]
                   << setw(15)
                   << item.xyzCoordinate[2]
                   << "\n";
        totalCharge += specToCharge[item.atomName];
    }
    for (auto& guestMolecule : guestMolecules) {
        double charge = 0.0;
        int atomCount = 0;
        for (int i = 0; i < guestMolecule.numberOfMolecule * guestMolecule.atomsPerMolecule; i++) {
            indexNumber++;
            if (i % guestMolecule.atomsPerMolecule == 0) {
                belongingIndex++;
            }
            auto atomNumber = generalToIndex[guestMolecule.atoms[i]];
            auto x = guestMolecule.coordinates[i][0];
            auto y = guestMolecule.coordinates[i][1];
            auto z = guestMolecule.coordinates[i][2];

            outputFile << right
                       << setw(5)
                       << indexNumber
                       << setw(5)
                       << belongingIndex
                       << setw(5)
                       << atomNumber
                       << setw(15)
                       << fixed
                       << setprecision(7)
                       << charge
                       << setw(15)
                       << x
                       << setw(15)
                       << y
                       << setw(15)
                       << z
                       << "\n";
        }
    }
    cout << "Total charge is " << fixed << scientific << setprecision(5) << totalCharge << ". Be sure that makes sense." << endl;

    stringstream numberStream;
    numberStream << setw(12) << indexNumber;
    atomLine = numberStream.str() + " atoms" + " # MOF " + to_string(referenceConfig.numberOfAtoms);
    for (auto& guestMolecule : guestMolecules) {
        if (guestMolecule.moleculeName == "h2o")
            atomLine += " water " + to_string(guestMolecule.numberOfMolecule * guestMolecule.atomsPerMolecule);
        else if (guestMolecule.moleculeName == "co2")
            atomLine += " CO2 " + to_string(guestMolecule.numberOfMolecule * guestMolecule.atomsPerMolecule);
    }

    outputFile << "\n";
    outputFile.close();
}

void Processor::processMassPart() {
    ofstream outputFile("MASS");
    int MOFCount = 0;

    cout << "Starting Mass part" << endl;
    outputFile << "Masses\n"
               << "\n";

    for (int i = 1; i <= indexToGeneral.size(); i++) {
        auto atom = indexToGeneral[i];
        auto mass = generalToMass[atom];

        outputFile << right << setw(5) << i << fixed << setprecision(5) << setw(16) << mass << " # ";

        if (atom == "OW") {
            outputFile << "O of water\n";
        } else if (atom == "HW") {
            outputFile << "H of water\n";
        } else if (atom == "CC") {
            outputFile << "C of CO2\n";
        } else if (atom == "OC") {
            outputFile << "O of CO2\n";
        } else {
            outputFile << atom << "\n";
            MOFCount++;
        }
    }

    stringstream numberStream;
    numberStream << setw(12) << indexToGeneral.size();
    massLine = numberStream.str() + " atom types # MOF " + to_string(MOFCount);
    for (auto& guestMolecule : guestMolecules) {
        if (guestMolecule.moleculeName == "h2o")
            massLine += " water 2";
        else if (guestMolecule.moleculeName == "co2")
            massLine += " CO2 2";
    }

    outputFile << "\n";
    outputFile.close();
}

void Processor::prepareMetalPairs() {
    MetalPairsReader metalPairsReader;

    this->metalPairs = metalPairsReader.metalPairs;
}

void Processor::changeCutoffIfNeeded(const string& atom1, const string& atom2, double& cutoff) {
    for (auto& item : metalPairs) {
        if ((item[0] == atom1 && item[1] == atom2) || (item[0] == atom2 && item[1] == atom1)) {
            cutoff = 2.8;
        }
    }
}

void Processor::findGuestMolecules(const vector<string>& fileNames) {
    string line;

    for (auto& fileName : fileNames) {
        auto guestMolecule = GuestMolecule(fileName);
        for (auto& existingMolecule : this->guestMolecules) {
            if (existingMolecule.moleculeName == guestMolecule.moleculeName) {
                line = "** Error: Duplicated type of guest molecule for file " + fileName + " and " + existingMolecule.moleculeName;
                throw runtime_error(line);
            }
        }
        this->guestMolecules.push_back(guestMolecule);
    }
}

void Processor::findProcessedForceFieldAtoms() {
    unordered_set<string> referenceAtomSet;

    for (auto& item : referenceConfig.referenceAtoms) {
        referenceAtomSet.insert(item.atomName);
    }

    for (auto& item : forceField.atoms) {
        if (referenceAtomSet.contains(item.specificAtomType)) {
            processedForceFieldAtoms.push_back(item);
        }
    }

    for (auto& guestMolecule : guestMolecules) {
        if (guestMolecule.moleculeName == "h2o") {
            processedForceFieldAtoms.push_back(ForceFieldAtom("OW", "OW", 0.0, 15.99490));
            processedForceFieldAtoms.push_back(ForceFieldAtom("HW", "HW", 0.0, 1.007900));
        } else if (guestMolecule.moleculeName == "co2") {
            processedForceFieldAtoms.push_back(ForceFieldAtom("CC", "CC", 0.0, 12.01100));
            processedForceFieldAtoms.push_back(ForceFieldAtom("OC", "OC", 0.0, 15.99490));
        }
    }
}

void Processor::findSpecificCorrespondance() {
    vector<string> existingGeneralAtoms;
    int count = 0;

    for (auto& item : processedForceFieldAtoms) {
        specToGeneral[item.specificAtomType] = item.generalAtomType;
        specToCharge[item.specificAtomType] = item.charge;

        bool isPresent = find(existingGeneralAtoms.begin(), existingGeneralAtoms.end(), item.generalAtomType) != existingGeneralAtoms.end();
        if (!isPresent) {
            existingGeneralAtoms.push_back(item.generalAtomType);
            count++;
            generalToIndex[item.generalAtomType] = count;
            indexToGeneral[count] = item.generalAtomType;
            generalToMass[item.generalAtomType] = item.mass;
        }
    }
}