#ifndef SSE_V2_GEOMETRY_H
#define SSE_V2_GEOMETRY_H
#include <vector>

class Geometry{
private:
    std::vector<std::vector<int>> sites;
    int num_sites;
    int num_bonds;

public:
    Geometry(const int &N, const std::vector<std::vector<int>> &bonds_to_sites);
};

#endif //SSE_V2_GEOMETRY_H
