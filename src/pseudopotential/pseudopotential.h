#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Eigen/dense>

std::vector<std::string> splitStringToWords(const std::string& str) {
    std::istringstream iss(str);
    std::vector<std::string> words;
    std::string word;

    while (iss >> word) {
        words.push_back(word);
    }

    return words;
}

class GoedeckerPseudopotential {

public:
    std::string element_name;
    int zeff;
    int zion;
    double rloc;
    double alpha_pp;
    int nloc;
    Eigen::VectorXd c;
    std::vector<double> rl;
    int lmax;
    std::vector<Eigen::MatrixXd> h;
    int nsep

    GoedeckerPseudopotential(std::string filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }
        std::string line;
        std::vector<std::string> words;
        // first line is irrelevant
        std::getline(file, line);
        // read zatom and zeff
        std::getline(file, line);
        words = splitStringToWords(line);
        zion = words[0];
        zeff = words[1];
        // read rloc, nloc, and the c coefficients
        std::getline(file, line);
        words = splitStringToWords(line);
        rloc = words[0];
        nloc = words[1];
        c = Eigen::VectorXd::Zero(nloc);
        for (int i = 0; i < nloc; i++) {
            c[i] = words[2 + i];
        }
        // read nsep
        std::getline(file, line);
        words = splitStringToWords(line);
        nsep = words[0];
        // check if nsep is between 1 and 3
        if (nsep < 1 || nsep > 3) {
            std::cerr << "Error: nsep must be between 1 and 3" << std::endl;
            return;
        }
        
        // read s projector if file not over yet:
        if (!std::getline(file, line)) return;
        words = splitStringToWords(line);
        rl.push_back(words[0]);
        h.push_back(Eigen::MatrixXd::Zero(nsep, nsep));
        for (int i = 0; i < nsep; i++) {
            h[0](0, i) = words[1 + i];
            h[0](i, 0) = h[0](0, i);
        }
        for (int i = 1; i < nsep; i++) {
            std::getline(file, line);
            words = splitStringToWords(line);
            for (int j = i; j < nsep; j++){
                h[0](i, j) = words[j];
                h[0](j, i) = h[0](i, j);
            }
        }
        // read p projector if file not over yet:
        if (!std::getline(file, line)) return;
        words = splitStringToWords(line);
        rl.push_back(words[0]);
        h.push_back(Eigen::MatrixXd::Zero(nsep, nsep));
        for (int i = 0; i < nsep; i++) {
            h[0](0, i) = words[1 + i];
            h[0](i, 0) = h[0](0, i);
        }
        for (int i = 1; i < nsep; i++) {
            std::getline(file, line);
            words = splitStringToWords(line);
            for (int j = i; j < nsep; j++){
                h[0](i, j) = words[j];
                h[0](j, i) = h[0](i, j);
            }
        }
        // read d projector if file not over yet:
        if (!std::getline(file, line)) return;
        words = splitStringToWords(line);
        rl.push_back(words[0]);
        h.push_back(Eigen::MatrixXd::Zero(nsep, nsep));
        for (int i = 0; i < nsep; i++) {
            h[0](0, i) = words[1 + i];
            h[0](i, 0) = h[0](0, i);
        }
        for (int i = 1; i < nsep; i++) {
            std::getline(file, line);
            words = splitStringToWords(line);
            for (int j = i; j < nsep; j++){
                h[0](i, j) = words[j];
                h[0](j, i) = h[0](i, j);
            }
        }
    }
};

