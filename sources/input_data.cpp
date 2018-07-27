#include "input_data.h"

void Input_data::read_input(std::ifstream& input_file) {
	if (!input_file.is_open())
		throw std::runtime_error("Input file is not open!");

	std::string line;
	while (std::getline(input_file, line)) {
		boost::tokenizer<> tok(line);
		boost::tokenizer<>::iterator beg = tok.begin();
		keys[*beg] = *beg++;
	}
}

