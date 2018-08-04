#include "input_data.h"


void Input_data::read_input(std::ifstream &input_file)
{
	if (!input_file.is_open())
		throw std::runtime_error("Input file is not open!");

	std::string line;
	while (std::getline(input_file, line))
	{
		if(line.empty()) continue;
		boost::tokenizer<boost::char_separator<char>> tok(line, boost::char_separator<char>(" \t\n"));
		auto beg = tok.begin();
		std::string key = *beg;
		beg++;
		std::string val = *beg;

		keys.insert(std::make_pair(key, val));
	}
}
