#ifndef GoedeckerPseudopotential_h
#define GoedeckerPseudopotential_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

struct GoedeckerPseudopotential {
    std::string element_name;
    int zatom;
    int zion;
    double rloc;
    double alpha_pp;
    std::vector<double> nloc;
    std::vector<double> c;
    std::vector<double> rs;
    std::vector<double> ns;
    std::vector<double> hs;
    std::vector<double> rp;
    std::vector<double> np;
    std::vector<double> hp;
    std::vector<double> kp;
    std::vector<double> rd;
    std::vector<double> nd;
    std::vector<double> hd;
    std::vector<double> kd;
};

bool readGoedeckerPseudopotential(const std::string& filename, GoedeckerPseudopotential& pseudopotential) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::string line;
    // Read element name
    if (!std::getline(file, line)) {
        std::cerr << "Error: Unable to read element name from " << filename << std::endl;
        return false;
    }
    pseudopotential.element_name = line;

    // Read atomic number, charge, and pseudopotential type
    if (!std::getline(file, line)) {
        std::cerr << "Error: Unable to read atomic number, charge, and pseudopotential type from " << filename << std::endl;
        return false;
    }
    std::istringstream iss1(line);
    if (!(iss1 >> pseudopotential.zatom >> pseudopotential.zion)) {
        std::cerr << "Error: Invalid format in " << filename << std::endl;
        return false;
    }

    // Skip the line with pspcod, pspxc, lmax, lloc, mmax, r2well
    if (!std::getline(file, line)) {
        std::cerr << "Error: Unable to read pspcod, pspxc, lmax, lloc, mmax, r2well from " << filename << std::endl;
        return false;
    }

    // Read rloc, nloc, c
    double rloc_val, nloc_val, c_val;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> rloc_val >> nloc_val >> c_val)) {
            break; // End of the section
        }
        pseudopotential.rloc = rloc_val;
        pseudopotential.alpha_pp = 1.0 / (std::sqrt(2.0) * rloc_val);
        pseudopotential.nloc.push_back(nloc_val);
        pseudopotential.c.push_back(c_val);
    }

    // Read rs, ns, hs
    double rs_val, ns_val, hs_val;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> rs_val >> ns_val >> hs_val)) {
            break; // End of the section
        }
        pseudopotential.rs.push_back(rs_val);
        pseudopotential.ns.push_back(ns_val);
        pseudopotential.hs.push_back(hs_val);
    }

    // Read rp, np, hp, kp
    double rp_val, np_val, hp_val, kp_val;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> rp_val >> np_val >> hp_val >> kp_val)) {
            break; // End of the section
        }
        pseudopotential.rp.push_back(rp_val);
        pseudopotential.np.push_back(np_val);
        pseudopotential.hp.push_back(hp_val);
        pseudopotential.kp.push_back(kp_val);
    }

    // Read rd, nd, hd, kd
    double rd_val, nd_val, hd_val, kd_val;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> rd_val >> nd_val >> hd_val >> kd_val)) {
            break; // End of the section
        }
        pseudopotential.rd.push_back(rd_val);
        pseudopotential.nd.push_back(nd_val);
        pseudopotential.hd.push_back(hd_val);
        pseudopotential.kd.push_back(kd_val);
    }

    return true;
}

void buildCorePotential(const GoedeckerPseudopotential& pseudopotential) {
    std::cout << "Element Name: " << pseudopotential.element_name << std::endl;
    std::cout << "Atomic Number: " << pseudopotential.zatom << std::endl;
    std::cout << "Charge: " << pseudopotential.zion << std::endl;

    std::cout << "Core Potential:\n";
    std::cout << "rloc: " << pseudopotential.rloc << ", nloc: " << pseudopotential.nloc[0] << ", c: " << pseudopotential.c[0] << std::endl;
}

void performCalculation(const GoedeckerPseudopotential& pseudopotential) {
    std::cout << "Performing calculation using valence orbitals:\n";
    // You can implement your calculation using the valence orbitals here
    // For demonstration purposes, let's just print out the valence orbitals
    std::cout << "rs:\n";
    for (size_t i = 0; i < pseudopotential.rs.size(); ++i) {
        std::cout << "rs: " << pseudopotential.rs[i] << ", ns: " << pseudopotential.ns[i] << ", hs: " << pseudopotential.hs[i] << std::endl;
    }

    std::cout << "rp:\n";
    for (size_t i = 0; i < pseudopotential.rp.size(); ++i) {
        std::

cout << "rp: " << pseudopotential.rp[i] << ", np: " << pseudopotential.np[i] << ", hp: " << pseudopotential.hp[i] << ", kp: " << pseudopotential.kp[i] << std::endl;
    }

    std::cout << "rd:\n";
    for (size_t i = 0; i < pseudopotential.rd.size(); ++i) {
        std::cout << "rd: " << pseudopotential.rd[i] << ", nd: " << pseudopotential.nd[i] << ", hd: " << pseudopotential.hd[i] << ", kd: " << pseudopotential.kd[i] << std::endl;
    }
}

// int main() {
//     std::string filename = "goedecker_pp.dat"; // Change this to the path of your pseudopotential file
//     GoedeckerPseudopotential pseudopotential;
//     if (readGoedeckerPseudopotential(filename, pseudopotential)) {
//         buildCorePotential(pseudopotential);
//         performCalculation(pseudopotential);
//     }
//     return 0;
// }

#endif
