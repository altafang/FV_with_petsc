// This was copy-pasted from Ryan Davis's PADI software: https://github.com/rsdavis/ParallelDiffuseInterface-PADI
// It provides utilities for reading input parameters from a text file.
// The logging code has been removed.

#ifndef INPUT_TOOLS_H
#define INPUT_TOOLS_H

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

template <typename T> // primary template
inline void unpack(std::map<std::string, std::string> params, std::string name, T &parameter)
{
    MPI_Abort(PETSC_COMM_WORLD, 0);
}

template <> // explicit specialization for T = double
inline void unpack(std::map<std::string, std::string> hash, std::string name, double & parameter)
{
    std::map<std::string, std::string>::iterator it;
    it = hash.find(name);
    if (it==hash.end()) { // parameter not found
        PetscPrintf(PETSC_COMM_WORLD, "Parameter not found\n");
        MPI_Abort(PETSC_COMM_WORLD, 0);
    } else { // parameter found
        parameter = std::stod(it->second);
    }
}

template <> // explicit specialization for T = int
inline void unpack(std::map<std::string, std::string> hash, std::string name, int & parameter)
{
    std::map<std::string, std::string>::iterator it;
    it = hash.find(name);
    if (it==hash.end()) { // parameter not found
        PetscPrintf(PETSC_COMM_WORLD, "Parameter not found\n");
        MPI_Abort(PETSC_COMM_WORLD, 0);
    } else { // parameter found
        parameter = std::stoi(it->second);
    }
}

template <> // explicit specialization for T = float
inline void unpack(std::map<std::string, std::string> hash, std::string name, float & parameter)
{
    std::map<std::string, std::string>::iterator it;
    it = hash.find(name);
    if (it==hash.end()) { // parameter not found
        PetscPrintf(PETSC_COMM_WORLD, "Parameter not found\n");
        MPI_Abort(PETSC_COMM_WORLD, 0);
    } else { // parameter found
        parameter = std::stof(it->second);
    }
}

inline void unpack(std::map<std::string, int> name_index, std::string name, int &index)
{
    std::map<std::string, int>::iterator it;
    it = name_index.find(name);
    if (it==name_index.end()) { // parameter not found
        PetscPrintf(PETSC_COMM_WORLD, "Parameter not found\n");
        MPI_Abort(PETSC_COMM_WORLD, 0);
    } else { // parameter found
        index = it->second;
    }
}

void readParameters(std::map<std::string, std::string> &params)
{
    
    // This function parses the input file.
    // It reads a key-value pair from any line involving an "=" sign.
    // All key-value are stored as string in params for use later.
    // Comment lines begin with "#"
    
    std::ifstream input("input.txt");
    std::string line;
    
    while (std::getline(input, line))
    {
        size_t pos;
        std::string key, value;
        
        // remove whitespace before parsing
        line.erase(std::remove(line.begin(),line.end(),' '), line.end());
        line.erase(std::remove(line.begin(),line.end(),'\t'), line.end());
        line.erase(std::remove(line.begin(),line.end(),'\n'), line.end());
        line.erase(std::remove(line.begin(),line.end(),'\r'), line.end());
        line.erase(std::remove(line.begin(),line.end(),'\v'), line.end());
        
        // skip empty lines
        if (line.empty()) continue;
        
        // skip comments
        if (line[0] == '#') continue;
        
        // parse line using the equal sign
        pos = line.find("=");
        key = line.substr(0, pos);
        value = line.substr(pos+1, line.size());
        
        params[key] = value;
    }
    
    input.close();
}

#endif
