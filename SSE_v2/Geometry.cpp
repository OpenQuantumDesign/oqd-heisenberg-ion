#include "Geometry.h"

Geometry::Geometry(const int &N, const std::vector<std::vector<int>> &bonds_to_sites) {
    num_sites = N;
    num_bonds = (N * (N-1))/2;
    sites = bonds_to_sites;
}
