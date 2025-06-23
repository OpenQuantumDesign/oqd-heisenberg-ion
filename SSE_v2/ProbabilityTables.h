#ifndef SSE_V2_PROBABILITYTABLES_H
#define SSE_V2_PROBABILITYTABLES_H
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include "SimulationParameters.h"
#include "VertexTypes.h"

class ProbabilityTables {
private:
    std::string diagonal_probabilities_file;
    std::string vertex_weights_file;
    std::string loop_update_probabilities_file;
    std::string max_norm_probabilities_file;
    std::string geometry_file;

    std::string readTabularFile(const std::string &file_path, std::vector<std::vector<double>> &data_container);
    std::string readTabularFile(const std::string &file_path, std::vector<std::vector<int>> &data_container);
    std::string readVectorFile(const std::string &file_path, std::vector<double> &data_container);
    void extractHeaderEntry(const std::string &header_string, const std::string &identifier, double &member_to_update);
    void normalizeLoopProbs(const VertexTypes &vertex_types, const SimulationParameters &sim_params);

public:
    std::vector<std::vector<int>> lattice_sites;

    std::vector<std::vector<double>> diagonal_probabilities;
    std::vector<double> max_norm_probabilities;
    double max_diagonal_norm;

    std::vector<std::vector<double>> vertex_weights;

    std::vector<std::vector<double>> loop_update_probabilities;

    double spectrum_offset;

    explicit ProbabilityTables(const SimulationParameters &sim_params, const VertexTypes &vertex_types);
};


#endif //SSE_V2_PROBABILITYTABLES_H
