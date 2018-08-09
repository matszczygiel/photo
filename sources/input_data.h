#ifndef SOURCES_INPUT_DATA_H_
#define SOURCES_INPUT_DATA_H_

#include <boost/tokenizer.hpp>
#include <string>
#include <map>
#include <fstream>
#include <stdexcept>
#include <regex>

class Input_data
{
  public:
	void read_input(std::ifstream &input_file);

	std::map<std::string, std::string> keys;
};

#endif /* SOURCES_INPUT_DATA_H_ */
