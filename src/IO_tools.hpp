// Most of this was copy-pasted from Ryan Davis's PADI software:
// https://github.com/rsdavis/ParallelDiffuseInterface-PADI
// It provides utilities for reading input parameters from a text file.
// The logging code has been removed.

#ifndef IO_TOOLS_H
#define IO_TOOLS_H

#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>
#include <algorithm> // std::remove

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

// utility function for appending numbers to hdf5 filenames, if there is a series
std::string number_filename(std::string base_filename, int counter);

template <typename T> // primary template
void unpack(std::map<std::string, std::string> params, std::string name, T &parameter);

template <> // explicit specialization for T = double
void unpack(std::map<std::string, std::string> hash, std::string name, double & parameter);

template <> // explicit specialization for T = int
void unpack(std::map<std::string, std::string> hash, std::string name, int & parameter);

template <> // explicit specialization for T = float
void unpack(std::map<std::string, std::string> hash, std::string name, float & parameter);

void unpack(std::map<std::string, int> name_index, std::string name, int &index);

void read_parameters(std::map<std::string, std::string> &params, std::string input_file="input.txt");

#endif
