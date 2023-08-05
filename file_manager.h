#ifndef FILE_MANAGER_H
#define FILE_MANAGER_H

#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "progress.h"
#include "globals.h"

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}

class FileManager {
public:

    //reading data from input pdb file.
    int readTrajectory() {
        const std::string& filename = config.trajectoryFilename;
        if (DEBUG) {
            std::cout << "Reading file: " << filename << std::endl;
        }
        std::string line;
        std::ifstream file1(filename);
        int lines_count = 0;
        if (file1.is_open()) {
            while (getline(file1, line)) {
                lines_count++;
            }
            file1.close();
        } else {
            if (DEBUG) {
                std::cout << "Cannot find file: " << filename << std::endl;
            }
            return 1;
        }
        Progress p1(lines_count);
        std::ifstream file(filename);
        if (file.is_open()) {
            int frame = 0;
            int atom;
            SPHERES = 0;
            FRAMES = 0;
            ATOMS = 0;
            A = {};
            sphereCA = {};
            // sphereSize = {};
            while (getline(file, line)) {
                p1.improve();
                if (line[0] == 'M') {
                    frame = stoi(line.substr(9, 5));
                    frame--;
                    A.push_back({});
                    FRAMES++;
                } else if (line[0] == 'A') {
                    atom = stoi(line.substr(6, 5));
                    atom--;
                    A[frame].push_back({});
                    A[frame][atom].push_back(stod(line.substr(30, 8)));
                    A[frame][atom].push_back(stod(line.substr(38, 8)));
                    A[frame][atom].push_back(stod(line.substr(46, 8)));
                    if (frame == 0) {
                        ATOMS++;
                        if (line[14] == 'A' and line[13] == 'C') {
                            sphereCA.push_back(atom);
                            // sphereSize.push_back(0);
                            SPHERES++;
                        }
                    }
                }
            }
            p1.end();
            file.close();
            if (DEBUG) {
                std::cout << "File parsed" << std::endl;
            }
            return 0;
        } else {
            if (DEBUG) {
                std::cout << "Cannot find file: " << filename << std::endl;
            }
            return 1;
        }
    }

    bool readConfig(const std::string& filename) {
        std::ifstream file(filename);
        if (file.is_open()) {
            std::string line;
            std::unordered_map<std::string, std::string> configMap;
            while (getline(file, line)) {
                if (line[0] == '#') {
                    continue;
                }
                std::istringstream iss(line);
                std::string key, value;
                std::getline(iss, key, ':');
                std::getline(iss, value);
                trim(key);
                trim(value);
                configMap[key] = value;
            }

            if (configMap.find("trajectoryFilename") != configMap.end()) {
                config.trajectoryFilename = configMap["trajectoryFilename"];
            }
            if (configMap.find("matrixSize") != configMap.end()) {
                config.matrixSize = std::stoi(configMap["matrixSize"]);
            }
            if (configMap.find("timeLimitMinutes") != configMap.end()) {
                config.timeLimitMinutes = std::stod(configMap["timeLimitMinutes"]);
            }
            if (configMap.find("showDebugCurrentBest") != configMap.end()) {
                config.showDebugCurrentBest = configMap["showDebugCurrentBest"] == "true" ? true : false;
            }
            if (configMap.find("showDebugRouteBest") != configMap.end()) {
                config.showDebugRouteBest = configMap["showDebugRouteBest"] == "true" ? true : false;
            }
            if (configMap.find("jumpFromLocalAreaChance") != configMap.end()) {
                config.jumpFromLocalAreaChance = std::stod(configMap["jumpFromLocalAreaChance"]);
            }
            if (configMap.find("randomFrameWhileSwappingChance") != configMap.end()) {
                config.randomFrameWhileSwappingChance = std::stod(configMap["randomFrameWhileSwappingChance"]);
            }
            if (configMap.find("randomSeed") != configMap.end()) {
                config.randomSeed = configMap["randomSeed"] == "true" ? true : false;
            }
            if (configMap.find("ompThreadsPerCore") != configMap.end()) {
                config.ompThreadsPerCore = std::stod(configMap["ompThreadsPerCore"]);
            }
            if (configMap.find("memorySize") != configMap.end()) {
                config.memorySize = std::stod(configMap["memorySize"]);
            }
            if (configMap.find("writeAsCSV") != configMap.end()) {
                config.writeAsCSV = configMap["writeAsCSV"] == "true" ? true : false;
            }
            if (configMap.find("showLogs") != configMap.end()) {
                config.showLogs = configMap["showLogs"] == "true" ? true : false;
            }
            if (configMap.find("showRMSDCounter") != configMap.end()) {
                config.showRMSDCounter = configMap["showRMSDCounter"] == "true" ? true : false;
            }
            if (configMap.find("runRepetitions") != configMap.end()) {
                config.runRepetitions = std::stoi(configMap["runRepetitions"]);
            }
            
            DEBUG = config.showLogs;
            DEBUG_RMSD = config.showRMSDCounter;

            if (DEBUG) {
                std::cout << "Reading file: " << filename << std::endl;
            }

            file.close();
            if (DEBUG) {
                std::cout << "Config file loaded" << std::endl;
            }
            return true;
        } else {
            if (DEBUG) {
                std::cout << "Cannot find config file: " << filename << std::endl;
            }
            return false;
        }
    }

    static void writeResultsAsCSV(int bestI, int bestJ, double bestValue, double elapsedTime) {

        std::cout
            << "local_search;" 
            << config.trajectoryFilename << ";"
            << config.timeLimitMinutes << ";"
            << config.ompThreadsPerCore << ";"
            << config.runRepetitions << ";"

            << config.jumpFromLocalAreaChance << ";"
            << config.randomFrameWhileSwappingChance << ";"
            << config.memorySize << ";"

            << config.matrixSize << ";"

            << bestI << ";"
            << bestJ << ";"
            << bestValue << ";"
            << elapsedTime
        << std::endl;
    }
};

#endif // FILE_MANAGER_H
