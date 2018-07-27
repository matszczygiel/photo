#ifndef SOURCES_INPUT_DATA_H_
#define SOURCES_INPUT_DATA_H_

#include <boost/tokenizer.hpp>
#include <string>
#include <map>
#include <fstream>
#include <stdexcept>

class Input_data {
public:
	std::map<std::string, std::string> keys;
	void read_input(std::ifstream& input_file);
};

#endif /* SOURCES_INPUT_DATA_H_ */
