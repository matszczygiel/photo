#include "disk_reader.h"

void Disk_reader::initialize(const Job_control &jc)
{
    job = std::make_shared<Job_control>(jc);
    read_file_basis();
}

void Disk_reader::read_file_basis()
{
    std::ifstream bfile(job->get_file_basis());
    if (!bfile.is_open())
        throw std::runtime_error("Cannot open the basis file.");

    basis_l = 0;
    basis_lnk = 0;
    bool n_cont = true;
    std::string line, name;
    std::stringstream ss;
    do
    {
        getline(bfile, line);
    } while (line != "$BASIS");
    getline(bfile, line);
    do
    {
        ss.clear();
        ss << line;
        line.clear();
        ss >> name;
        name == job->get_continuum_id() ? n_cont = false : n_cont = true;
        do
        {
            int end;
            std::string moment;

            ss.clear();
            ss << line;
            line.clear();
            ss >> moment >> end;

            if (moment == "S")
            {
                basis_l += Const_arrays::crt_siz.at(0);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(0);
            }
            else if (moment == "P")
            {
                basis_l += Const_arrays::crt_siz.at(1);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(1);
            }
            else if (moment == "D")
            {
                basis_l += Const_arrays::crt_siz.at(2);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(2);
            }
            else if (moment == "F")
            {
                basis_l += Const_arrays::crt_siz.at(3);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(3);
            }
            else if (moment == "G")
            {
                basis_l += Const_arrays::crt_siz.at(4);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(4);
            }
            else if (moment == "H")
            {
                basis_l += Const_arrays::crt_siz.at(5);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(5);
            }
            else if (moment == "I")
            {
                basis_l += Const_arrays::crt_siz.at(6);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(6);
            }
            else if (moment == "K")
            {
                basis_l += Const_arrays::crt_siz.at(7);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(7);
            }
            else if (moment == "L")
            {
                basis_l += Const_arrays::crt_siz.at(8);
                if (n_cont)
                    basis_lnk += Const_arrays::crt_siz.at(8);
            }
            else
                throw std::runtime_error("The basis file contains unrecognized shell.");

            for (int i = 0; i < end; i++)
                getline(bfile, line);

            getline(bfile, line);
        } while (!line.empty());
        getline(bfile, line);
    } while (line != "$END");
    bfile.close();
}
