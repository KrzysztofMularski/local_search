#ifndef GLOBALS_H
#define GLOBALS_H

#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>
#include <omp.h>

extern bool DEBUG;
extern bool DEBUG_RMSD;
extern int RMSDCalculationCount;
extern int AllocationsCount;
extern bool AlreadyShowedRMSDCalculationCount;
extern int omp_thread_id;

extern double sphereRadius;
extern int SPHERES;
extern int ATOMS;
extern int FRAMES;
extern int FRAMEONE;
extern int FRAMETWO;

#pragma omp threadprivate(\
    RMSDCalculationCount,\
    AllocationsCount,\
    AlreadyShowedRMSDCalculationCount,\
    omp_thread_id,\
    FRAMEONE,\
    FRAMETWO)

// Atoms[<frame>][<atom>][<coordinate>]
extern std::vector<std::vector<std::vector<double>>> A;

// Maps sphere to CA; CAAtomNumber[<sphere>]
extern std::vector<int> sphereCA;

// List of atoms in [<sphere>]
extern std::vector<std::vector<int>>* sphereAtoms;

struct Config {
    std::string trajectoryFilename;             // trajectory filename
    int matrixSize;                             // analysing first [matrixSize] frames of pairs matrix
    double timeLimitMinutes;                    // max time for whole local search to finish
    bool showDebugCurrentBest;                  // showing current best value
    bool showDebugRouteBest;                    // showing current route best value
    double jumpFromLocalAreaChance;             // probability of jumping from local area
    double randomFrameWhileSwappingChance;      // probability of choosing random frame while swapping allocations
    bool randomSeed;                            // random seed for srand()
    double ompThreadsPerCore;                   // omp threads number per one cpu core
    double memorySize;                          // [0, 1] where 0 is no memory, and 1 is remembering whole matrix
    bool writeAsCSV;                            // each run of a program generates one line in CSV format
    bool showLogs;                              // show logs in the console
    bool showRMSDCounter;                       // show rsmd counter in the console
    int runRepetitions;                         // program execution repetition number

    void print() {
        if (!DEBUG) {
            return;
        }
        std::cout << "Config: " << std::endl;
        std::cout << " - " << "trajectoryFilename = " << trajectoryFilename << std::endl;
        std::cout << " - " << "matrixSize = " << matrixSize << std::endl;
        std::cout << " - " << "timeLimitMinutes = " << timeLimitMinutes << std::endl;
        std::cout << " - " << "showDebugCurrentBest = " << (showDebugCurrentBest ? "true" : "false") << std::endl;
        std::cout << " - " << "showDebugRouteBest = " << (showDebugRouteBest ? "true" : "false") << std::endl;
        std::cout << " - " << "jumpFromLocalAreaChance = " << jumpFromLocalAreaChance << std::endl;
        std::cout << " - " << "randomFrameWhileSwappingChance = " << randomFrameWhileSwappingChance << std::endl;
        std::cout << " - " << "randomSeed = " << (randomSeed ? "true" : "false") << std::endl;
        std::cout << " - " << "ompThreadsPerCore = " << ompThreadsPerCore << std::endl;
        std::cout << " - " << "memorySize = " << memorySize << std::endl;
        std::cout << " - " << "writeAsCSV = " << (writeAsCSV ? "true" : "false") << std::endl;
        std::cout << " - " << "showLogs = " << (showLogs ? "true" : "false") << std::endl;
        std::cout << " - " << "showRMSDCounter = " << (showRMSDCounter ? "true" : "false") << std::endl;
        std::cout << " - " << "runRepetitions = " << runRepetitions << std::endl;
    }

    void initDefault() {
        trajectoryFilename = "";
        timeLimitMinutes = 0.5;
        ompThreadsPerCore = 0;
        writeAsCSV = false;
        runRepetitions = 1;

        jumpFromLocalAreaChance = 0.1;
        randomFrameWhileSwappingChance = 0.01;
        memorySize = 0.1;

        randomSeed = true;
        matrixSize = -1;
        showLogs = true;
        showRMSDCounter = false;
        showDebugCurrentBest = true;
        showDebugRouteBest = false;
    }
};

extern Config config;

struct PairHash {
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U>& p) const {
        // Use std::hash for each element and combine the hashes
        return std::hash<int>()(p.first) * config.matrixSize + std::hash<int>()(p.second);
    }
};

extern std::unordered_set<std::pair<int, int>, PairHash> memorySet;

extern omp_lock_t memoryMutex;

extern int MEMORY_SIZE;

inline extern int getRandom(int offset, int range) {
    return offset + (rand() % (range + 1));
}

template <class T> inline extern void _debug(T t) {
    std::cout << t << std::endl;
}

template <class T, class... Args> inline extern void _debug(T t, Args... args) {
    std::cout << t;
    _debug(args...);
}

template <class T> inline extern void debug(T t) {
    if (DEBUG) {
        if (AlreadyShowedRMSDCalculationCount) {
            std::cout << std::endl;
        }
        AlreadyShowedRMSDCalculationCount = false;
        std::cout << "[DEBUG] " << t << std::endl;
    }
}

template <class T, class... Args> inline extern void debug(T t, Args... args) {
    if (DEBUG) {
        if (AlreadyShowedRMSDCalculationCount) {
            std::cout << std::endl;
        }
        AlreadyShowedRMSDCalculationCount = false;
        std::cout << "[DEBUG] " << t;
        _debug(args...);
    }
}

inline extern void debugRMSD() {
    if (DEBUG_RMSD) {
        if (AlreadyShowedRMSDCalculationCount) {
            std::cout << "\r";
            std::cout.flush();
        }
        std::cout << "[RMSD]: " << RMSDCalculationCount;
        AlreadyShowedRMSDCalculationCount = true;
    }
}

template <class T> inline extern void print(T t) {
    if (DEBUG) {
        std::cout << t << std::endl;
    }
}

template <class T, class... Args> inline extern void print(T t, Args... args) {
    if (DEBUG) {
        std::cout << t;
        print(args...);
    }
}

#endif // GLOBALS_H
