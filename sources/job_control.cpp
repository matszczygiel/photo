#include "job_control.h"

Job_control::Job_control()
{
    compute_bound = false;
    compute_cont = false;
    compute_ion_state = false;

    file_HFv = "c.dat";
    file_CI = "ci.dat";
    file_out = "res.out";

    continuum_id = "CONT";

    selection_m = energy;
    gauge = dipole;
    force_orth = false;

    write = false;
}

void Job_control::read(const Input_data &in)
{
    std::string key, arg;

    arg = "COMPUTE_BOUND";
    key = boost::to_lower_copy(in.keys[arg]);
    if (key == "y")
        compute_bound = true;
    else if (key == "n")
        compute_bound = false;
    else
        throw std::runtime_error("Invalid argument for  " + arg + ".");

    arg = "COMPUTE_CONTINUUM";
    key = boost::to_lower_copy(in.keys[arg]);
    if (key == "y")
        compute_cont = true;
    else if (key == "n")
        compute_cont = false;
    else
        throw std::runtime_error("Invalid argument for  " + arg + ".");

    arg = "COMPUTE_ION_STATE";
    key = boost::to_lower_copy(in.keys[arg]);
    if (key == "y")
        compute_ion_state = true;
    else if (key == "n")
        compute_ion_state = false;
    else
        throw std::runtime_error("Invalid argument for  " + arg + ".");

    arg = "K";
    key = in.keys[arg];
    try
    {
        k_val = std::stod(key);
    }
    catch (std::exception &e)
    {
        throw std::runtime_error("Invalid argument for  " + arg + ".");
    }

    arg = "ENERGY_I_STATE";
    key = in.keys[arg];
    try
    {
        energy_i = std::stod(key);
    }
    catch (std::exception &e)
    {
        throw std::runtime_error("Invalid argument for  " + arg + ".");
    }

    arg = "ENERGY_F_STATE";
    key = in.keys[arg];
    try
    {
        energy_f = std::stod(key);
    }
    catch (std::exception &e)
    {
        throw std::runtime_error("Invalid argument for  " + arg + ".");
    }

    arg = "PATH_IN";
    path_in = in.keys[arg];
    path_in += "/";

    arg = "FILE_BASIS";
    key = in.keys[arg];
    if (boost::to_lower_copy(key) == "auto")
    {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << k_val;
        file_basis = "He_k" + ss.str() + ".inp";
    }
    else
        file_basis = key;

    arg = "FILE_1E";
    key = in.keys[arg];
    if (boost::to_lower_copy(key) == "auto")
    {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << k_val;
        file_1E = "file1E_He_k" + ss.str() + ".F";
    }
    else
        file_1E = key;

    arg = "FILE_2E" key = in.keys[arg];
    if (boost::to_lower_copy(key) == "auto")
    {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << k_val;
        file_2E = "file2E_He_k" + ss.str() + ".F";
    }
    else
        file_2E = key;

    arg = "FILE_NORM";
    key = in.keys[arg];
    if (boost::to_lower_copy(key) == "auto")
    {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << k_val;
        file_norm = "norm1E_He_k" + ss.str() + ".F";
    }
    else
        file_norm = key;

    arg = "FILE_HF_VEC";
    key = in.keys[arg];
    if (!(boost::to_lower_copy(key) == "auto"))
        file_HFv = key;

    arg = "FILE_CI";
    key = in.keys[arg];
    if (!(boost::to_lower_copy(key) == "auto"))
        file_CI = key;

    arg = "CONTINUUM_ID";
    key = in.keys[arg];
    if (!(boost::to_lower_copy(key) == "auto"))
        continuum_id = key;

    arg = "ION_POTENTIAL";
    key = boost::to_lower_copy(in.keys[arg]);
    if (key == "compute")
        potential = energy_f - energy_i;
    else
    {
        try
        {
            potential = std::stod(key);
        }
        catch (std::exception &e)
        {
            throw std::runtime_error("Invalid argument for  " + arg + ".");
        }
    }

    arg = "SELECTION_METHOD";
    key = boost::to_lower_copy(in.keys[arg]);
    if (key == "energy")
        selection_m = energy;
    else if (key == "norm")
        selection_m = norm;
    else
        throw std::runtime_error("Invalid argument for  " + arg + ".");

    arg = "FORCE_ORTHOGONALITY";
    key = boost::to_lower_copy(in.keys[arg]);
    if (key == "y")
        force_orth = true;
    else if (key == "n")
        force_orth = false;
    else
        throw std::runtime_error("Invalid argument for  " + arg + ".");

    arg = "GAUGE";
    key = boost::to_lower_copy(in.keys[arg]);
    if (key == "velocity")
        gauge = velocity;
    else if (key == "dipole")
        gauge = dipole;
    else
        throw std::runtime_error("Invalid argument for  " + arg + ".");

    arg = "WRITE";
    key = boost::to_lower_copy(in.keys[arg]);
    if (key == "y")
        write = true;
    else if (key == "n")
        write = false;
    else
        throw std::runtime_error("Invalid argument for  " + arg + ".");

    arg = "PATH_OUT";
    path_out = in.keys[arg];
    path_out += "/";

    arg = "FILE_OUT";
    file_out = in.keys[arg];
}

/*
void JobControl::calcualteBasisFunctNumber()
{
    l_max = gto_number.size() - 1;
    b_k_l = 0;
    b_nk_l = 0;
    for (unsigned i = 0; i < gto_number.size(); i++)
    {
        b_nk_l += crt_siz[i] * gto_number[i];
        b_k_l += crt_siz[i];
    }
    b_l = b_k_l + b_nk_l;
    b_k_l_sqrt = b_k_l * b_k_l;
    b_nk_l_sqrt = b_nk_l * b_nk_l;
    b_l_sqrt = b_l * b_l;
}

void JobControl::print() const
{
    std::cout << "\n\n"
              << " Number of basis set functions: " << b_l << "\n";
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

void JobControl::writeRes(double sig) const
{
    std::ofstream outfile;
    outfile.open(data_path + res_file, std::ios_base::app);
    outfile << std::fixed << std::setprecision(3) << k << "\t";
    outfile << std::fixed << std::setprecision(3) << photonEeV(k, ionizatoin_pot) << "\t";
    outfile << std::fixed << std::setprecision(4) << sig << "\t\t";
    outfile << ion_pot_int << "\t";
    outfile << forceOrthogonality << "\t";
    outfile << eq_int << "\t";
    outfile << solveIon << "\t";
    outfile << selectionMethod << "\t";
    outfile << gauge << "\n";
    outfile.close();
}

void JobControl::setInFiles(double kk)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << kk;
    std::string k_ = stream.str();
    file1E = "file1E_He_k" + k_ + ".F";
    file2E = "file2E_He_k" + k_ + ".F";
    norms_file = "norm1E_He_k" + k_ + ".F";
}

void JobControl::readInput(std::string in_file)
{

    basis_name = in_file;
    basis_name.erase(basis_name.end() - 4, basis_name.end());
    std::ifstream file(in_file);
    if (file.is_open())
    {

        std::stringstream ss;
        std::string line, key;
        int ask;

        while (getline(file, line))
        {
            ss.clear();
            ss << line;
            line.clear();
            ss >> key;

            if (key == "EQ")
            {
                ss >> ask;
                eq_int = ask;
                switch (ask)
                {
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
            else if (key == "DATA_PATH")
            {
                ss >> data_path;
                data_path = "/home/mateusz/Documents/cpp/photo/data/" + data_path;
            }
            else if (key == "ION")
            {
                ss >> ask;
                ion_pot_int = ask;
                switch (ask)
                {
                case 0:
                    ionizatoin_pot = 0.9035694; //reference value
                    break;
                case 1:
                    ionizatoin_pot = E_ion - E_he; //computed value
                    break;
                case 2:
                {
                    std::ifstream file(data_path + HF_en_file);
                    if (file.is_open())
                    {
                        double e;
                        file >> e;
                        ionizatoin_pot = -e; //Koopmans theorem value
                    }
                    else
                    {
                        std::cerr << " Cannot open HF energies file! \n";
                        std::cout << " Using default settings. \n";
                        ionizatoin_pot = 0.9035694; //reference value
                    }
                    file.close();
                }
                break;
                default:
                    std::cout << " Using default settings \n";
                    ionizatoin_pot = 0.9035694; //reference value
                    break;
                }
            }
            else if (key == "K")
            {
                ss >> k;
                this->setInFiles(k);
            }
            else if (key == "E_HF")
                ss >> E_he;
            else if (key == "E_HF_ION")
                ss >> E_ion;
            else if (key == "L")
            {
                ask = 0;
                for (unsigned i = 0; i < 8; i++)
                {
                    ss >> ask;
                    if (ask == -1)
                        break;
                    gto_number.push_back(ask);
                }
            }
            else if (key == "FORCE_ORTHOGONALITY")
            {
                ss >> ask;
                if (ask == 1)
                    forceOrthogonality = true;
                else if (ask == 0)
                    forceOrthogonality = false;
                else
                {
                    std::cout << " Using default orgthogonality forceing settings. \n";
                    forceOrthogonality = true;
                }
            }
            else if (key == "WRITE")
            {
                ss >> ask;
                if (ask == 1)
                    writeResult = true;
                else if (ask == 0)
                    writeResult = false;
                else
                {
                    std::cout << " Using default writing settings. \n";
                    writeResult = false;
                }
            }
            else if (key == "SOLVE_ION")
            {
                ss >> ask;
                if (ask == 1)
                    solveIon = true;
                else if (ask == 0)
                    solveIon = false;
                else
                {
                    std::cout << " Using default writing settings. \n";
                    solveIon = false;
                }
            }
            else if (key == "SELECTION_METHOD")
            {
                ss >> selectionMethod;
            }
            else if (key == "WRITE_TO")
            {
                ss >> res_file;
            }
            else if (key == "GAUGE")
            {
                ss >> gauge;
            }
        }
    }
    else
    {
        std::cout << "Cannot open input file! \n";
        std::exit(EXIT_FAILURE);
    }
    file.close();

    if (check)
    {
        double photonEV;
        std::cout << " Put photon energy. \n";
        std::cin >> photonEV;
        std::cout << " Length of vector k: ";
        double kk = std::sqrt(2.0 * (photonEV / 27.211385 - ionizatoin_pot));
        std::cout << std::fixed << std::setprecision(3) << kk << "\n\n";
        std::exit(EXIT_SUCCESS);
    }
    std::cout << " All necesary information gathered \n\n";
}
*/