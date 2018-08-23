#include "basis.h"

#include <sstream>

const std::map<char, int> Shell::charmap = {{'S', 0},
                                            {'P', 1},
                                            {'D', 2},
                                            {'F', 3},
                                            {'G', 4},
                                            {'H', 5},
                                            {'I', 6}};

const std::array<int, 11> Shell::crt_siz = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55};
const std::array<int, 11> Shell::labels  = {'S', 'P', 'D', 'F', 'G', 'H', 'I'};

Shell char2shell(const char &c) {
    int shlind = Shell::charmap.at(c);
    return Shell(shlind);
}

char shell2char(const Shell &shell) {
    return Shell::labels.at(static_cast<int>(shell.shl));
}

int shell2int(const Shell &shell) {
    return static_cast<int>(shell.shl);
}

std::ostream &operator<<(std::ostream &os, const GTOPW &rhs) {
    os << shell2char(rhs.shl);
    os.width(3);
    os << rhs.size;
    os << "\n";
    std::string spaces = "          ";

    for (int i = 0; i < rhs.size; ++i) {
        os.width(3);
        os << i + 1;
        os << std::scientific;
        os.precision(9);
        os << spaces;
        os.width(18);
        os << rhs.exps[i];
        os << spaces;
        os.width(18);
        os << rhs.coefs[i].real();
        os.width(18);
        os << rhs.coefs[i].imag();
        os << spaces;
        os << std::fixed;
        os.precision(5);
        os.width(10);
        os << rhs.k[0];
        os.width(10);
        os << rhs.k[1];
        os.width(10);
        os << rhs.k[2];
        os << "\n";
    }
    return os;
}

bool GTOPW::read(std::istream &is) {
    std::string line;
    if (!getline(is, line))
        return false;
    if (line.empty())
        return false;

    std::istringstream ss(line);
    char moment = ss.get();
    shl         = char2shell(moment);
    ss >> size;
    exps.clear();
    coefs.clear();
    exps.reserve(size);
    coefs.reserve(size);

    for (int i = 0; i < size; ++i) {
        if (!getline(is, line))
            throw std::runtime_error("Invalind gtopw contraction read - check file.");
        std::istringstream ssl(line);
        double ex, re, im;
        int gnum;
        ssl >> gnum;
        if (gnum != i + 1)
            throw std::runtime_error("Invalind gtopw contraction read.");

        ssl >> ex >> re >> im;
        exps.push_back(ex);
        coefs.push_back(cdouble(re, im));

        double k0, k1, k2;
        ssl >> k0 >> k1 >> k2;
        if (i == 0)
            k = {k0, k1, k2};
        else if (k0 != k[0] || k1 != k[1] || k2 != k[2])
            throw std::runtime_error("Invalind gtopw contraction read - check k.");
    }

    return true;
}

int GTOPW::functions_number() const {
    return Shell::crt_siz.at(shell2int(shl)) * size;
}

bool Basis::read(std::istream &is) {
    Basis();

    std::string line;
    if (!getline(is, line))
        return false;
    if (line.empty())
        return false;

    std::istringstream ss(line);
    std::string token;
    ss >> token;
    if (token == "$END")
        return false;

    label = token;
    ss >> charge;
    ss >> position[0] >> position[1] >> position[2];
    GTOPW g;
    while (g.read(is))
        gtopws.push_back(g);

    return true;
}

std::ostream &operator<<(std::ostream &os, const Basis &rhs) {
    std::string spaces = "          ";
    os << rhs.label;
    os.width(7);
    os.precision(2);
    os << std::fixed;
    os << rhs.charge;
    os << spaces;
    os.precision(5);
    os.width(10);
    os << rhs.position[0];
    os.width(10);
    os << rhs.position[1];
    os.width(10);
    os << rhs.position[2];
    os << "\n";

    for (const auto &x : rhs.gtopws)
        os << x;

    os << "\n";
    return os;
}

int Basis::functions_number() const {
    int num = 0;
    for (const auto &x : gtopws)
        num += x.functions_number();

    return num;
}

Shell Basis::get_max_shell() const {
    int max = 0;
    for (auto it = gtopws.begin(); it != gtopws.end(); it++) {
        auto shl = shell2int(it->get_shell());
        if (shl > max)
            max = shl;
    }
    return Shell(max);
}

void Basis::truncate_at(const Shell &shl) {
    for (auto it = gtopws.begin(); it != gtopws.end();) {
        if (shell2int(it->get_shell()) > shell2int(shl))
            gtopws.erase(it);
        else
            it++;
    }
}

void Basis::set_kvec(const Vec3d &kvec) {
    for (auto &x : gtopws)
        x.set_kvec(kvec);
}

Vec3d Basis::get_kvec() const {
    return gtopws.at(0).get_kvec();
}