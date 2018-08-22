#include "input_data.h"

#include <regex>
#include <stdexcept>

void Input_data::read_input(std::ifstream &input_file) {
    if (!input_file.is_open())
        throw std::runtime_error("Input file is not open!");

    std::string line;
    while (std::getline(input_file, line)) {
        if (line.empty())
            continue;

        std::regex reg("\\s+");

        std::sregex_token_iterator beg(line.begin(), line.end(), reg, -1);
        std::sregex_token_iterator end;

        std::vector<std::string> vec(beg, end);
        keys.insert(std::make_pair(vec.at(0), vec.at(1)));
    }
}
