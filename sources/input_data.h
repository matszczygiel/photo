#ifndef SOURCES_INPUT_DATA_H
#define SOURCES_INPUT_DATA_H

#include <string>
#include <map>
#include <fstream>

class Input_data
{
  public:
	void read_input(std::ifstream &input_file);

	std::map<std::string, std::string> keys;
};

#endif /* SOURCES_INPUT_DATA_H_ */
