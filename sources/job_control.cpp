#include "job_control.h"

Job_control::Job_control(const Input_data& in) {
	std::string key;

	key = boost::to_lower_copy(in.keys["SOLVE_BOUND"]);
	if (key == "y")
		compute_bound = true;
	else if (key == "n")
		compute_bound = false;
	else
		throw std::runtime_error("Invalid argument for SOLVE_BOUND.");

	key = boost::to_lower_copy(in.keys["SOLVE_CONTINUUM"]);
	if (key == "y")
		compute_cont = true;
	else if (key == "n")
		compute_cont = false;
	else
		throw std::runtime_error("Invalid argument for SOLVE_CONTINUUM.");

	key = boost::to_lower_copy(in.keys["COMPUTE_ION_STATE"]);
	if (key == "y")
		compute_ion_state = true;
	else if (key == "n")
		compute_ion_state = false;
	else
		throw std::runtime_error("Invalid argument for COMPUTE_ION_STATE.");

	key = in.keys["K"];
	try {
		k_val = std::stod(key);
	} catch (std::exception& e) {
		throw std::runtime_error("Invalid argument for K.");
	}

	key = in.keys["I_STATE_ENERGY"];
	try {
		i_energy = std::stod(key);
	} catch (std::exception& e) {
		throw std::runtime_error("Invalid argument for I_STATE_ENERGY.");
	}

	key = in.keys["F_STATE_ENERGY"];
	try {
		f_energy = std::stod(key);
	} catch (std::exception& e) {
		throw std::runtime_error("Invalid argument for F_STATE_ENERGY.");
	}

	data_path = in.keys["DATA_PATH"];

	key = in.keys["BASIS_FILE"];
	if (boost::to_lower_copy(key) == "auto") {
		std::stringstream ss;
		ss << std::fixed << std::setprecision(3) << k_val;
		basis_file = "He_k" + ss.str() + ".inp";
	}
	else
		basis_file = key;

	key = in.keys["1E_FILE"];
	if (boost::to_lower_copy(key) == "auto") {
		std::stringstream ss;
		ss << std::fixed << std::setprecision(3) << k_val;
		E1_file = "file1E_He_k" + ss.str() + ".F";
	} else
		E1_file = key;

	key = in.keys["2E_FILE"];
	if (boost::to_lower_copy(key) == "auto") {
		std::stringstream ss;
		ss << std::fixed << std::setprecision(3) << k_val;
		E2_file = "file2E_He_k" + ss.str() + ".F";
	} else
		E2_file = key;

	key = in.keys["NORM_FILE"];
	if (boost::to_lower_copy(key) == "auto") {
		std::stringstream ss;
		ss << std::fixed << std::setprecision(3) << k_val;
		norm_file = "norm1E_He_k" + ss.str() + ".F";
	} else
		norm_file = key;

	key = in.keys["HF_VECTORS_FILE"];
	if (boost::to_lower_copy(key) == "auto") {
		HFv_file = "c.dat";
	} else
		HFv_file = key;

	key = in.keys["CI_FILE"];
	if (boost::to_lower_copy(key) == "auto") {
		CI_file = "ci.dat";
	} else
		CI_file = key;

	key = boost::to_lower_copy(in.keys["ION_POTENTIAL"]);
	if (key == "compute")
		potential = f_energy - i_energy;
	else
		try {
			potential = std::stod(key);
		} catch (std::exception& e) {
			throw std::runtime_error("Invalid argument for ION_POTENTIAL.");
		};

	key = boost::to_lower_copy(in.keys["SELECTION_METHOD"]);
	if (key == "energy")
		selection_m = energy;
	else if (key == "norm")
		selection_m = norm;
	else
		std::runtime_error("Invalid argument for SELECTION_METHOD.");

	key = boost::to_lower_copy(in.keys["FORCE_ORTHOGONALITY"]);
	if (key == "y")
		force_orth = true;
	else if (key == "n")
		force_orth = false;
	else
		throw std::runtime_error("Invalid argument for FORCE_ORTHOGONALITY.");

	key = boost::to_lower_copy(in.keys["GAUGE"]);
	if (key == "velocity")
		gauge = velocity;
	else if (key == "dipole")
		gauge = dipole;
	else
		std::runtime_error("Invalid argument for GAUGE.");

	key = boost::to_lower_copy(in.keys["WRITE"]);
	if (key == "y")
		write = true;
	else if (key == "n")
		write = false;
	else
		throw std::runtime_error("Invalid argument for WRITE.");

	res_file = in.keys["WRITE_FILE"];
}












void JobControl::calcualteBasisFunctNumber() {
    l_max = gto_number.size() - 1;
    b_k_l = 0;
    b_nk_l =0;
    for(unsigned i =0; i<gto_number.size(); i++) {
        b_nk_l += crt_siz[i] * gto_number[i];
        b_k_l += crt_siz[i];
    }
    b_l = b_k_l + b_nk_l;
    b_k_l_sqrt = b_k_l * b_k_l;
    b_nk_l_sqrt = b_nk_l * b_nk_l;
    b_l_sqrt = b_l * b_l;
}


void JobControl::print() const {
    std::cout << "\n\n" << " Number of basis set functions: " << b_l << "\n";
    std::cout << " Number of pure GTO's: " << b_nk_l << "\n";
    std::cout << " Number of PW-GTO's: " << b_k_l << "\n";
    std::cout << " Max l value: " << l_max << "\n";
    std::cout << " The binary files loaded: \n";
    std::cout << this->getFile1E() << "\n";
    std::cout << this->getFile2E() << "\n";
    std::cout << this->getFileNorm() << "\n";
    std::cout << " The text files loaded: \n";
    std::cout << this->getFileHF() << "\n";
    std::cout << this->getFileHFEnergy() << "\n";
    std::cout << this->getFileCI() << "\n";
    std::cout << " Ionization potential: " << ionizatoin_pot << "\n";
    std::cout << " Loaded k value: " << k << "\n";
    std::cout << " The He atom energy: " << E_he << "\n\n";

}


void JobControl::writeRes(double sig) const {
    std::ofstream outfile;
    outfile.open(data_path + res_file, std::ios_base::app);
    outfile << std::fixed <<std::setprecision(3) << k <<"\t";
    outfile << std::fixed  <<std::setprecision(3) << photonEeV(k, ionizatoin_pot) << "\t";
    outfile << std::fixed << std::setprecision(4) << sig << "\t\t";
    outfile << ion_pot_int << "\t";
    outfile << forceOrthogonality << "\t";
    outfile << eq_int << "\t";
    outfile << solveIon << "\t";
    outfile << selectionMethod << "\t";
    outfile << gauge << "\n";
    outfile.close();
}


void JobControl::setInFiles( double kk ) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << kk;
    std::string k_ = stream.str();
    file1E = "file1E_He_k" + k_ + ".F";
    file2E = "file2E_He_k" + k_ + ".F";
    norms_file = "norm1E_He_k" + k_ + ".F";
}


void JobControl::readInput( std::string in_file ) {

    basis_name = in_file;
    basis_name.erase(basis_name.end() - 4, basis_name.end());
    std::ifstream file(in_file);
    if(file.is_open()) {

        std::stringstream ss;
        std::string line, key;
        int ask;

        while ( getline ( file, line )) {
            ss.clear();
            ss << line;
            line.clear();
            ss >> key;

            if(key == "EQ") {
                ss >> ask;
                eq_int = ask;
                switch (ask) {
                case 0:
                    check = false;
                    iter = true;
                    eval_1eq = true;
                    eval_2eq = true;
                    break;
                case 1:
                    check = false;
                    iter = true;
                    eval_1eq = true;
                    eval_2eq = false;
                    break;
                case 2:
                    check = false;
                    iter = true;
                    eval_1eq = false;
                    eval_2eq = true;
                    break;
                case 3:
                    check = false;
                    iter = false;
                    eval_1eq = false;
                    eval_2eq = false;
                    break;
                case 4:
                    check = true;
                    iter = false;
                    eval_1eq = false;
                    eval_2eq = false;
                    break;
                default:
                    std::cout << " Using default settings \n";
                    check = false;
                    iter = true;
                    eval_1eq = true;
                    eval_2eq = true;
                    break;
                }
            }
            else if(key == "DATA_PATH")
            {
                ss >> data_path;
                data_path = "/home/mateusz/Documents/cpp/photo/data/" + data_path;
            }
            else if(key == "ION") {
                ss >> ask;
                ion_pot_int = ask;
                switch (ask) {
                case 0:
                    ionizatoin_pot = 0.9035694; //reference value
                    break;
                case 1:
                    ionizatoin_pot = E_ion - E_he; //computed value
                    break;
                case 2: {
                    std::ifstream file(data_path + HF_en_file);
                    if(file.is_open()){
                        double e;
                        file >> e;
                        ionizatoin_pot = -e; //Koopmans theorem value
                    }
                    else {
                        std::cerr << " Cannot open HF energies file! \n";
                        std::cout << " Using default settings. \n";
                        ionizatoin_pot = 0.9035694; //reference value
                    }
                    file.close();
                } break;
                default:
                    std::cout << " Using default settings \n";
                    ionizatoin_pot = 0.9035694; //reference value
                    break;
                }
            }
            else if(key == "K") {
                ss >> k;
                this->setInFiles(k);
            }
            else if(key == "E_HF") ss >> E_he;
            else if(key == "E_HF_ION") ss >> E_ion;
            else if(key == "L") {
                ask =0;
                for(unsigned i=0; i<8; i++){
                    ss >> ask;
                    if(ask == -1) break;
                    gto_number.push_back(ask);
                }
            }
            else if ( key == "FORCE_ORTHOGONALITY" )
            {
                ss >> ask;
                if ( ask == 1 ) forceOrthogonality = true;
                else if ( ask == 0 ) forceOrthogonality = false;
                else
                {
                    std:: cout << " Using default orgthogonality forceing settings. \n";
                    forceOrthogonality = true;
                }
            }
            else if ( key == "WRITE" )
            {
                ss >> ask;
                if ( ask == 1 ) writeResult = true;
                else if ( ask == 0 ) writeResult = false;
                else
                {
                    std:: cout << " Using default writing settings. \n";
                    writeResult = false;
                }
            }
            else if ( key == "SOLVE_ION" )
            {
                ss >> ask;
                if ( ask == 1 ) solveIon = true;
                else if ( ask == 0 ) solveIon = false;
                else
                {
                    std:: cout << " Using default writing settings. \n";
                    solveIon = false;
                }
            }
            else if ( key == "SELECTION_METHOD" )
            {
                ss >> selectionMethod;
            }
            else if (key == "WRITE_TO") {
            	ss >> res_file;
            }
            else if (key == "GAUGE") {
            	ss >> gauge;
            }
        }
    }
    else {
        std::cout << "Cannot open input file! \n";
        std::exit(EXIT_FAILURE);
    }
    file.close();

    if(check)    {
        double photonEV;
        std::cout << " Put photon energy. \n";
        std::cin >> photonEV;
        std::cout << " Length of vector k: ";
        double kk = std::sqrt(2.0 * (photonEV/27.211385 - ionizatoin_pot));
        std::cout << std::fixed << std::setprecision(3)<< kk << "\n\n";
        std::exit(EXIT_SUCCESS);
    }
    std::cout << " All necesary information gathered \n\n";
}



















