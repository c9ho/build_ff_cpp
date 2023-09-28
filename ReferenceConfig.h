#ifndef REFERENCECONFIG_H
#define REFERENCECONFIG_H

#include <sstream>
#include <string>
#include <vector>

class ReferenceAtom {
   public:
    std::string atomName;
    double xyzCoordinate[3];
    double abcCoordinate[3];
    int atomIndex;

    ReferenceAtom(std::string& atomName, double xyzCoordinate[3], double abcCorrdinate[3], int index);
};

class ReferenceConfig {
   public:
    int numberOfAtoms;
    double box[3][3];
    std::vector<ReferenceAtom> referenceAtoms;

    ReferenceConfig();

    double calcDistance(const double abcCoor1[3], const double abcCoor2[3]) const;
    double calcAngle(const double abcCoor1[3], const double abcCoor2[3], const double abcCoor3[3]) const;
    static void getZeroToOne(double numbers[3]);
    static void getFirstBZone(double numbers[3]);

   private:
    const double PI = 3.141592653589793238L;
    double reciprocalVectors[3][3];
    void prepareNumberOfAtoms(std::istringstream& iss);
    void prepareBox(std::istringstream& iss);
    void prepareAtoms(std::istringstream& iss, const int atomIndex);
    void prepareReciprocalVectors();
};

#endif