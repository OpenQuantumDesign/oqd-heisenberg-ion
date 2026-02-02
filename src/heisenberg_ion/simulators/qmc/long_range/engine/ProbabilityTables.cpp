#include "ProbabilityTables.h"

std::string ProbabilityTables::readTabularFile(const std::string &file_path,
                                               std::vector<std::vector<double>> &data_container){

    std::ifstream file_stream(file_path);
    if (!file_stream.good())
        throw std::runtime_error("Could not open file with path: " + file_path);

    std::string line, entry;
    std::vector<double> row_data;

    int lineCount = 0;
    std::string header;
    while (std::getline(file_stream, line)) {

        lineCount += 1;
        if (lineCount == 1) {
            header = line;
            header.erase(std::remove_if(header.begin(), header.end(),
                [](unsigned char c){ return (c == '#' || c == ' '); }), header.end());
            std::cout << header << std::endl;
            continue;
        }

        std::stringstream row_stream(line);

        row_data.clear();
        while (std::getline(row_stream, entry, ',')) {
            double table_val = std::stod(entry);
            if (isnan(table_val)) {
                throw std::runtime_error("Invalid value encountered\n");
            }
            else if (table_val < 0.0){
                if (table_val > -1e-15){
                    throw std::runtime_error("Small negative value encountered\n");
                }
                else {
                    throw std::runtime_error("Large negative value encountered\n");
                }
            }
            row_data.push_back(table_val);
        }

        data_container.push_back(row_data);
    }
    file_stream.close();

    return header;
}

std::string ProbabilityTables::readTabularFile(const std::string &file_path,
                                               std::vector<std::vector<int>> &data_container,
                                               const bool &negatives_allowed) {

    std::ifstream file_stream(file_path);
    if (!file_stream.good())
        throw std::runtime_error("Could not open file with path: " + file_path);

    std::string line, entry;
    std::vector<int> row_data;

    int lineCount = 0;
    std::string header;
    while (std::getline(file_stream, line)) {

        lineCount += 1;
        if (lineCount == 1) {
            header = line;
            header.erase(std::remove_if(header.begin(), header.end(),
                [](unsigned char c){ return (c == '#' || c == ' '); }), header.end());
            std::cout << header << std::endl;
            continue;
        }

        std::stringstream rowStream(line);

        row_data.clear();
        while (std::getline(rowStream, entry, ',')) {
            int table_val = std::stoi(entry);
            if (isnan(table_val)) {
                throw std::runtime_error("Invalid value encountered\n");
            }
            else if (table_val < 0 and !negatives_allowed) {
                throw std::runtime_error("Negative value encountered\n");
            }
            row_data.push_back(table_val);
        }

        data_container.push_back(row_data);
    }
    file_stream.close();

    return header;
}

std::string ProbabilityTables::readVectorFile(const std::string &file_path, std::vector<double> &data_container) {

    std::ifstream file_stream(file_path);
    if (!file_stream.good())
        throw std::runtime_error("Could not open file with path: " + file_path);

    std::string line, entry;
    double row_data;

    int lineCount = 0;
    std::string header;
    while (std::getline(file_stream, line)) {

        lineCount += 1;
        if (lineCount == 1) {
            header = line;
            header.erase(std::remove_if(header.begin(), header.end(),
                [](unsigned char c){ return (c == '#' || c == ' '); }), header.end());
            std::cout << header << std::endl;
            continue;
        }

        std::stringstream rowStream(line);

        row_data = std::stod(line);

        if (isnan(row_data)) {
            throw std::runtime_error("Invalid value encountered\n");
        }
        else if (row_data < 0.0){
            if (row_data > -1e-15){
                throw std::runtime_error("Small negative value encountered\n");
            }
            else {
                throw std::runtime_error("Large negative value encountered\n");
            }
        }
        data_container.push_back(row_data);
    }
    file_stream.close();

    return header;
}

void ProbabilityTables::extractHeaderEntry(const std::string &header_string, const std::string &identifier,
                                           double &member_to_update) {

    std::string key_val_pair;
    std::stringstream ss_1(header_string);
    bool norm_found = false;
    member_to_update = 0.0;

    while(std::getline(ss_1, key_val_pair, ','))
    {
        if (norm_found) {
            break;
        }

        std::stringstream ss_2(key_val_pair);
        std::string kvp_segment;

        while(std::getline(ss_2, kvp_segment, '=')) {

            if (norm_found) {
                member_to_update = std::stod(kvp_segment);
            }
            else if (kvp_segment == identifier){
                norm_found = true;
            }
        }
    }

    if (!norm_found) {
        throw std::runtime_error("Value not found in probability tables for " + identifier + "\n");
    }

}

ProbabilityTables::ProbabilityTables(const SimulationParameters &sim_params, const VertexTypes &vertex_types) {

    if (sim_params.hamiltonian_type == 2){

        diagonal_probabilities_file = sim_params.simulation_subfolder + "/probability_densities/diag_probs.csv";
        vertex_weights_file = sim_params.simulation_subfolder + "/probability_densities/vertex_weights.csv";
        loop_update_probabilities_file = sim_params.simulation_subfolder + "/probability_densities/" +
                sim_params.loop_type + "_off_diag_table.csv";

        std::string header_diag_prob_file = readTabularFile(diagonal_probabilities_file,
                                                            diagonal_probabilities);
        std::string header_loop_update_file = readTabularFile(loop_update_probabilities_file,
                                                              loop_update_probabilities);
        std::string  header_vertex_weights_file = readTabularFile(vertex_weights_file,
                                                                  vertex_weights);

        normalizeLoopProbs(vertex_types, sim_params);
    }

    max_norm_probabilities_file = sim_params.simulation_subfolder + "/probability_densities/max_over_states.csv";
    geometry_file = sim_params.simulation_subfolder + "/probability_densities/geometry.csv";

    std::string header_max_norm_file = readVectorFile(max_norm_probabilities_file,
                                                      max_norm_probabilities);

    extractHeaderEntry(header_max_norm_file, "norm", max_diagonal_norm);
    extractHeaderEntry(header_max_norm_file, "spectrum_offset", spectrum_offset);

    readTabularFile(geometry_file, lattice_sites, true);
}

void ProbabilityTables::normalizeLoopProbs(const VertexTypes &vertex_types, const SimulationParameters &sim_params) {

    int a_1 = vertex_types.num_legs_per_vertex*(vertex_types.num_legs_per_vertex-1);
    int a_2 = (vertex_types.num_legs_per_vertex-1);

    int a_3 = vertex_types.num_legs_per_vertex * vertex_types.num_legs_per_vertex;
    int a_4 = vertex_types.num_legs_per_vertex;

    for (int b=0; b<sim_params.num_bonds; b++){
        for (int i=0; i<vertex_types.num_vertices; i++){
            for (int j=0; j<vertex_types.num_legs_per_vertex; j++){
                double prob_sum = 0.0;
                int c_1 = i*a_1 + j*a_2;
                int c_2 = i*a_3 + j*a_4;
                for (int k=0; k<vertex_types.num_legs_per_vertex-1; k++){
                    int l_x = vertex_types.allowed_exit_legs.at(c_1 + k);
                    double prob = loop_update_probabilities.at( c_2 + l_x).at(b);
                    prob_sum += prob;
                }
                if (prob_sum != 0){
                    for (int k=0; k<vertex_types.num_legs_per_vertex-1; k++){
                        int l_x = vertex_types.allowed_exit_legs.at(c_1 + k);
                        loop_update_probabilities.at(c_2 + l_x).at(b) =
                                loop_update_probabilities.at(c_2 + l_x).at(b)/prob_sum;
                    }
                }
            }
        }
    }
}