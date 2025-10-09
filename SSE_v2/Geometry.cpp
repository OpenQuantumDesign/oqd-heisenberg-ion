#include<vector>

class Geometry{
private:
    std::vector<std::vector<int>> sites;
    int num_sites;
    int num_bonds;

public:
    Geometry(const int &N, const std::vector<std::vector<int>> &bonds_to_sites){
        num_sites = N;
        num_bonds = (N * (N-1))/2;
        sites = bonds_to_sites;
    };
};
