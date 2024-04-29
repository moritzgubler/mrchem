#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Eigen/Dense>

std::vector<std::string> splitStringToWords(const std::string& str) {
    std::istringstream iss(str);
    std::vector<std::string> words;
    std::string word;

    while (iss >> word) {
        words.push_back(word);
    }

    return words;
}

class PseudopotentialData {

public:
    int zeff;
    int zion;
    double rloc;
    double alpha_pp;
    int nloc;
    Eigen::VectorXd c;
    std::vector<double> rl;
    int lmax;
    int nProjectors;
    std::vector<Eigen::MatrixXd> h;
    int nsep;

    PseudopotentialData(std::string filename) {
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
        zion = std::stoi(words[0]);
        zeff = std::stoi(words[1]);

        // third position is lmax
        std::getline(file, line);
        words = splitStringToWords(line);
        lmax = std::stoi(words[2]);
        nProjectors = lmax + 1;

        // read rloc, nloc, and the c coefficients
        std::getline(file, line);
        words = splitStringToWords(line);
        rloc = std::stod(words[0]);
        nloc = std::stoi(words[1]);
        alpha_pp = 1.0 / (std::sqrt(2.0) * rloc);

        c = Eigen::VectorXd::Zero(nloc);
        for (int i = 0; i < nloc; i++) {
            c[i] = std::stod(words[2 + i]);
        }
        // read nsep
        std::getline(file, line);
        words = splitStringToWords(line);
        nsep = std::stoi(words[0]);
        // check if nsep is between 1 and 3
        if (nsep < 0 || nsep > 3) {
            std::cerr << "Error: nsep must be between 1 and 3" << std::endl;
            return;
        }

        for (int l = 0; l < nProjectors; l++) {
            if (!std::getline(file, line)) {
                std::cerr << "Error: Could not read projector" << std::endl;
                return;
            }
            words = splitStringToWords(line);
            rl.push_back(std::stod(words[0]));
            h.push_back(Eigen::MatrixXd::Zero(nsep, nsep));
            for (int i = 0; i < nsep; i++) {
                h[l](0, i) = std::stod(words[1 + i]);
                h[l](i, 0) = h[l](0, i);
            }
            for (int i = 1; i < nsep; i++) {
                std::getline(file, line);
                words = splitStringToWords(line);
                for (int j = i; j < nsep; j++){
                    h[l](i, j) = std::stod(words[j-i]);
                    h[l](j, i) = h[l](i, j);
                }
            }
        }
    }

    void print() {
        std::cout << "Zeff: " << zeff << std::endl;
        std::cout << "Zion: " << zion << std::endl;
        std::cout << "rloc: " << rloc << std::endl;
        std::cout << "alpha_pp: " << alpha_pp << std::endl;
        std::cout << "nloc: " << nloc << std::endl;
        std::cout << "c: " << c.transpose() << std::endl;
        std::cout << "nsep: " << nsep << std::endl;
        std::cout << "lmax " << lmax << std::endl;

        for (int l = 0; l < nProjectors; l++)
        {
            std::cout << "l: " << l << std::endl;
            std::cout << "rl: " << rl[l] << std::endl;
            std::cout << "h: " << h[l] << std::endl;
        }
    }

};

