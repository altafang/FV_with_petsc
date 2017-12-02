// Most of this was copy-pasted from Ryan Davis's PADI software:
// https://github.com/rsdavis/ParallelDiffuseInterface-PADI
// It provides utilities for reading input parameters from a text file.
// The logging code has been removed.

#ifndef IO_TOOLS_H
#define IO_TOOLS_H

#include <map>
#include <string>
#include "bc.hpp"

// Process boundary condition types input from the input file.
BC_type convert_to_BCtype(std::string name);

// utility function for appending numbers to hdf5 filenames, if there is a series
std::string number_filename(std::string base_filename, int counter);

template <typename T> // primary template
void unpack(std::map<std::string, std::string> params, std::string name, T &parameter);

template <> // explicit specialization for T = double
void unpack(std::map<std::string, std::string> hash, std::string name, double &parameter);

template <> // explicit specialization for T = int
void unpack(std::map<std::string, std::string> hash, std::string name, int &parameter);

template <> // explicit specialization for T = BC
void unpack(std::map<std::string, std::string> hash, std::string name, BC &parameter);

void unpack(std::map<std::string, int> name_index, std::string name, int &index);

void read_parameters(std::map<std::string, std::string> &params, 
                     std::string input_file="input.txt");

#endif
