#ifndef METALPAIRSREADER_H
#define METALPAIRSREADER_H

#include <array>
#include <ostream>
#include <vector>

class MetalPairsReader {
   public:
    std::vector<std::array<std::string, 2>> metalPairs;
    MetalPairsReader();
};

std::ostream& operator<<(std::ostream& stream, const std::array<std::string, 2>& array);

#endif