/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#pragma once

#include "NuclearFunction.h"

#include "utils/math_utils.h"
#include "pseudopotential/HGH.hpp"

namespace mrchem {

class PPNucleus : public NuclearFunction {
public:

    PPNucleus(std::vector<GoedeckerPseudopotential> pps) : NuclearFunction(){
        this->pps = pps;
    }
    

    // zero order, just take constant
    double evalf(const mrcpp::Coord<3> &r) const override {
        double result = 0.0;
        double temp_exp;
        double temp_square;
        // std::cerr << this->c1 << " " << this->c2 << " " << this->c3 << " " << this->c4 << std::endl;
        for (int i = 0; i < this->nuclei.size(); i++) {
            const auto &R = this->nuclei[i].getCoord();
            auto R1 = math_utils::calc_distance(R, r);

            temp_exp = std::exp(-this->pps[i].alpha_pp * this->pps[i].alpha_pp * R1 * R1);

            result -= this->pps[i].zatom / R1 * std::erf(this->pps[i].alpha_pp * R1);


            if (this->pps[i].c.size() > 0)
            {
                result += this->pps[i].c[0] * temp_exp;
            }

            if (this->pps[i].c.size() > 1)
            {
                temp_square = 2.0 * R1 * R1 * this->pps[i].alpha_pp * this->pps[i].alpha_pp;
                result += this->pps[i].c[1] * temp_exp * temp_square;
            }
            if (this->pps[i].c.size() > 2)
            {
                temp_square = temp_square * temp_square;
                result += this->pps[i].c[2] * temp_exp * temp_square;
            }
            if (this->pps[i].c.size() > 3)
            {
                temp_square = temp_square * 2.0 * R1 * R1 * this->pps[i].alpha_pp * this->pps[i].alpha_pp;
                result += this->pps[i].c[3] * temp_exp * temp_square;
            }
            
        }
        return result;
    }

    std::string getParamName1() const { return "Precision"; }
    std::string getParamName2() const { return "Smoothing"; }
    double calcParam1(double prec, const Nucleus &nuc) const { return prec; }
    double calcParam2(double prec, const Nucleus &nuc) const {
        double tmp = 0.0;
        return prec;
    }

    protected:
    std::vector<GoedeckerPseudopotential> pps;

};

} // namespace mrchem
