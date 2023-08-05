#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <stdexcept>

#include "RMSD_calculation.h"
#include "file_manager.h"
#include "globals.h"

bool DEBUG = true;
bool DEBUG_RMSD = false;
int RMSDCalculationCount = 0;
int AllocationsCount = 0;
bool AlreadyShowedRMSDCalculationCount = false;
int omp_thread_id;

double sphereRadius = 8;
int SPHERES;
int ATOMS;
int FRAMES;
int FRAMEONE;
int FRAMETWO;

// Atoms[<frame>][<atom>][<coordinate>]
std::vector<std::vector<std::vector<double>>> A;

// Maps sphere to CA; CAAtomNumber[<sphere>]
std::vector<int> sphereCA;

// List of atoms in [<sphere>]
std::vector<std::vector<int>> *sphereAtoms;

Config config;

std::unordered_set<std::pair<int, int>, PairHash> memorySet;

omp_lock_t memoryMutex;

int MEMORY_SIZE;

class LocalSearch {
  public:
    struct LocalSearchResult {
        double rmsdValue;
        int i;
        int j;
        LocalSearchResult() : rmsdValue(-1), i(-1), j(-1) {}
        LocalSearchResult(double rmsdValue, int i, int j) : rmsdValue(rmsdValue), i(i), j(j) {}
    };
    LocalSearchResult bestResult;
    RMSDCalculation rmsd;
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds> start;

    LocalSearch() {
        if (config.matrixSize == -1) {
            config.matrixSize = FRAMES;
        }
    }

    void choosePairRandom(int &i, int &j) {
        i = getRandom(0, config.matrixSize - 1);
        j = getRandom(0, config.matrixSize - 1);
        while (i == j) {
            j = getRandom(0, config.matrixSize - 1);
        }
    }

    inline bool saveIfRouteBest(LocalSearchResult &routeBest, double value, int i, int j) {
        if (value > routeBest.rmsdValue) {
            routeBest = {
                value,
                i,
                j,
            };
            if (config.showDebugRouteBest) {
                debug("[Current route best]: [", i, ", ", j, "] = ", value);
            }
            return true;
        }
        return false;
    }

    inline void saveIfBest(double value, int i, int j) {
        if (value > bestResult.rmsdValue) {
            bestResult = {
                value,
                i,
                j,
            };
            bestResult.j = j;
            bestResult.rmsdValue = value;
            if (config.showDebugCurrentBest) {
                debug("[Current best]: [", i, ", ", j, "] = ", value);
            }
        }
    }

    inline bool insideMatrixBoundaries(int &i) {
        return i >= 0 && i < config.matrixSize;
    }

    inline double changeAllocationsAndCalculate(int &allocatedOnFrame, int &changingFrame) {
        if (getRandom(1, 100) <= config.randomFrameWhileSwappingChance * 100) {
            // (A, B) -> (C, D)
            int new_i = getRandom(0, config.matrixSize - 1);
            int new_j = getRandom(0, config.matrixSize - 1);
            while (new_i == allocatedOnFrame || new_i == changingFrame) {
                new_i = getRandom(0, config.matrixSize - 1);
            }
            while (new_j == allocatedOnFrame || new_j == changingFrame || new_j == new_i) {
                new_j = getRandom(0, config.matrixSize - 1);
            }
            allocatedOnFrame = new_i;
            changingFrame = new_j;
            rmsd.atomsAllocation(allocatedOnFrame);

            return rmsd.calculateRMSDSuperpose(changingFrame);
        }
        // (A, B) -> (B, C)
        int new_j = getRandom(0, config.matrixSize - 1);
        while (new_j == allocatedOnFrame || new_j == changingFrame) {
            new_j = getRandom(0, config.matrixSize - 1);
        }
        allocatedOnFrame = changingFrame;
        changingFrame = new_j;
        rmsd.atomsAllocation(allocatedOnFrame);
        return rmsd.calculateRMSDSuperpose(changingFrame);
    }

    bool jump(int &allocatedOnFrame, int &changingFrame, LocalSearchResult &routeBest) {
        int new_j = getRandom(0, config.matrixSize - 1);

        while (new_j == allocatedOnFrame || new_j == changingFrame) {
            new_j = getRandom(0, config.matrixSize - 1);
        }

        double newValue = rmsd.calculateRMSDSuperpose(new_j);

        if (saveIfRouteBest(routeBest, newValue, allocatedOnFrame, new_j)) {
            changingFrame = new_j;
            return true;
        }
        return false;
    }

    inline bool identifiersGood(int &i, int &j) {
        return insideMatrixBoundaries(j) && i != j;
    }

    LocalSearchResult traverse(int i, int j) {
        // traversing route: i != j checked
        int allocatedOnFrame = i;
        int changingFrame = j;
        rmsd.atomsAllocation(allocatedOnFrame);
        LocalSearchResult routeBest = {
            rmsd.calculateRMSDSuperpose(changingFrame),
            allocatedOnFrame,
            changingFrame,
        };

        bool changedSidesAlready = false;
        int step = getRandom(0, 1) * 2 - 1; // -1 or +1

        while (true) {
            auto stop = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = stop - start;
            if (elapsed.count() > config.timeLimitMinutes * 60) {
                return routeBest;
            }

            int newChangingFrame = changingFrame + step;
            if (!identifiersGood(allocatedOnFrame, newChangingFrame)) {
                // changing direction
                if (!changedSidesAlready) {
                    step *= -1;
                    changedSidesAlready = true;
                    continue;
                }
                step = getRandom(0, 1) * 2 - 1; // -1 or +1
                changedSidesAlready = false;
                // otherwise, try to jump
                if (getRandom(1, 100) <= config.jumpFromLocalAreaChance * 100 && jump(allocatedOnFrame, changingFrame, routeBest)) {
                    // jump
                    continue;
                }
                // otherwise, no jump, change allocations
                double newValue = changeAllocationsAndCalculate(allocatedOnFrame, changingFrame);

                if (saveIfRouteBest(routeBest, newValue, allocatedOnFrame, changingFrame)) {
                    continue;
                } else {
                    return routeBest;
                }
            }

            double newValue = rmsd.calculateRMSDSuperpose(newChangingFrame);
            if (saveIfRouteBest(routeBest, newValue, allocatedOnFrame, newChangingFrame)) {
                // better than current best
                // going in straight line from now on
                changingFrame = newChangingFrame;
                while (true) {
                    newChangingFrame = changingFrame + step;
                    if (!identifiersGood(allocatedOnFrame, newChangingFrame)) {
                        break;
                    }
                    newValue = rmsd.calculateRMSDSuperpose(newChangingFrame);

                    if (saveIfRouteBest(routeBest, newValue, allocatedOnFrame, newChangingFrame)) {
                        changingFrame = newChangingFrame;
                        continue;
                    } else {
                        break;
                    }
                }
                // cannot go straight line anymore
                // cannot change direction because we came from there
                // trying to jump
                step = getRandom(0, 1) * 2 - 1; // -1 or +1
                changedSidesAlready = false;
                if (getRandom(1, 100) <= config.jumpFromLocalAreaChance * 100 && jump(allocatedOnFrame, changingFrame, routeBest)) {
                    // jump
                    continue;
                }
                // otherwise, no jump, change allocations
                double newValue = changeAllocationsAndCalculate(allocatedOnFrame, changingFrame);

                if (saveIfRouteBest(routeBest, newValue, allocatedOnFrame, changingFrame)) {
                    continue;
                } else {
                    return routeBest;
                }
            }

            // not better
            // trying to change sides
            if (!changedSidesAlready) {
                step *= -1;
                changedSidesAlready = true;
                continue;
            }

            // cannot change directions
            step = getRandom(0, 1) * 2 - 1; // -1 or +1
            changedSidesAlready = false;
            // otherwise, try to jump
            if (getRandom(1, 100) <= config.jumpFromLocalAreaChance * 100 && jump(allocatedOnFrame, changingFrame, routeBest)) {
                // jump
                continue;
            }
            // otherwise, no jump, change allocations
            newValue = changeAllocationsAndCalculate(allocatedOnFrame, changingFrame);

            if (saveIfRouteBest(routeBest, newValue, allocatedOnFrame, changingFrame)) {
                continue;
            } else {
                return routeBest;
            }
        }
    }

    void run() {
        omp_set_num_threads(omp_get_num_procs() * config.ompThreadsPerCore);
        int AllocationsCountGlobal = 0;
        int RMSDCalculationCountGlobal = 0;

        start = std::chrono::steady_clock::now();
        std::atomic<bool> time_exceeded(false);

#pragma omp parallel
        {
            omp_thread_id = omp_get_thread_num();
            if (omp_thread_id == 0) {
                debug("[OMP] [Number of threads]: ", omp_get_num_threads());
                sphereAtoms = new std::vector<std::vector<int>>[omp_get_num_threads()];
                // testSpecificPairs();
            }

#pragma omp barrier

            while (!time_exceeded) {
                // one route
                int i, j;
                choosePairRandom(i, j);
                // debug("[New route] (", i, ", ", j, ")");

                LocalSearchResult routeBest = traverse(i, j);
                saveIfBest(routeBest.rmsdValue, routeBest.i, routeBest.j);

                if (omp_thread_id == 0) {
                    auto stop = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed = stop - start;
                    if (elapsed.count() > config.timeLimitMinutes * 60) {
                        time_exceeded.store(true, std::memory_order_relaxed);
                        break;
                    }
                }
            }
#pragma omp atomic
            RMSDCalculationCountGlobal += RMSDCalculationCount;
#pragma omp atomic
            AllocationsCountGlobal += AllocationsCount;
        }

        delete[] sphereAtoms;

        auto stop = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = stop - start;
        print("Local Search Results:");
        print(" - Computation time: ", elapsed.count(), "s");
        print(" - RMSD counted: ", RMSDCalculationCountGlobal, " times.");
        print(" - Atoms allocated: ", AllocationsCountGlobal, " times.");

        if (config.writeAsCSV) {
            FileManager::writeResultsAsCSV(bestResult.i, bestResult.j, bestResult.rmsdValue, elapsed.count());
        }

        return;
    }
};

void resetGlobals() {
    RMSDCalculationCount = 0;
    AllocationsCount = 0;
    AlreadyShowedRMSDCalculationCount = false;
    memorySet.clear();
}

// Function to parse a value of type T from a string
template <typename T> T parseValue(const std::string &str) {
    T value;
    std::istringstream iss(str);
    if (!(iss >> value)) {
        throw std::runtime_error("Failed to parse value from string: " + str);
    }
    return value;
}

// Function to parse a boolean value from a string
bool parseBoolean(const std::string &str) {
    std::string lowerStr = str;
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), ::tolower);
    if (lowerStr == "true" || lowerStr == "t" || lowerStr == "1" || lowerStr == "yes" || lowerStr == "y" || lowerStr == "on" || lowerStr == "") {
        return true;
    } else if (lowerStr == "false" || lowerStr == "f" || lowerStr == "0" || lowerStr == "no" || lowerStr == "n" || lowerStr == "off") {
        return false;
    } else {
        throw std::runtime_error("Failed to parse boolean value from string: " + str);
    }
}

int readArgs(int argc, char *argv[], FileManager &fileManager) {

    if (argc == 1) {
        std::cout << "local_search: too few arguments" << std::endl;
        std::cout << "Try 'local_search --help' for more information." << std::endl;
        return 1;
    } else if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        std::cout << "Usage: local_search [OPTION]..." << std::endl;
        std::cout << std::endl;
        std::cout << "Performing local search algorithm to find the greatest deviation of provided trajectory frames by calculating RMSD value between them."
                  << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -c CONFIG                           provide config parameters via CONFIG file" << std::endl;
        std::cout << "  -h, --help                          display this message and exit" << std::endl;
        std::cout << std::endl;
        std::cout << "Parameters if no config provided. In descriptions: [type:default] format is used," << std::endl;
        std::cout << "where type is the type of parameter, and default is its default value." << std::endl;
        std::cout << "  --trajectory=TRAJECTORY             [string:] [mandatory] trajectory filename in .pdb format" << std::endl;
        std::cout << "  --time-limit=TIME                   [double:1.0] max time in minutes for whole local search to finish" << std::endl;
        std::cout << "  --omp-threads=NUM                   [double:0] omp threads number per one cpu core" << std::endl;
        std::cout << "  --write-as-csv=[true/false]         [bool:false] each run of a program generates one line in CSV format" << std::endl;
        std::cout << "  --repetitions=REPS                  [int:2] number of program executions" << std::endl;

        std::cout << "  --jump-chance=PROB                  [double:0.1] probability of jumping from local area" << std::endl;
        std::cout << "  --random-frame-chance=PROB          [double:0.01] probability of choosing random frame while swapping allocations" << std::endl;
        std::cout << "  --memory-size=SIZE                  [double:0.1] [0, 1] where 0 is no memory, and 1 is remembering whole matrix" << std::endl;

        std::cout << "  --random-seed=[true/false]          [bool:true] random seed for srand()" << std::endl;
        std::cout << "  --matrix-size=SIZE                  [int:-1] limiting matrix to SIZE by SIZE, if -1 then SIZE is max for current trajectory file" << std::endl;
        std::cout << "  --show-logs=[true/false]            [bool:true] show any logs in the console" << std::endl;
        std::cout << "  --show-rmsd-counter=[true/false]    [bool:false] show rsmd counter in the console" << std::endl;
        std::cout << "  --show-current-best=[true/false]    [bool:true] show current best value, works only if --show-logs is set" << std::endl;
        std::cout << "  --show-route-best=[true/false]      [bool:false] show current route best value, works only if --show-logs is set" << std::endl;
        std::cout << std::endl;
        std::cout << "Examples:" << std::endl;
        std::cout << "  local_search -c config.yml" << std::endl;
        std::cout << "  local_search --trajectory=traj.pdb --time-limit=0.5 --repetitions=5" << std::endl;
        std::cout << std::endl;
        std::cout << "All bool possible values:" << std::endl;
        std::cout << "  maps to true:  [true]  [t] [1] [yes] [y] [on]  []" << std::endl;
        std::cout << "  maps to false: [false] [f] [0] [no]  [n] [off]" << std::endl;
        return 1;
        
    } else if (argc == 3 && strcmp(argv[1], "-c") == 0) {
        std::string configFilename = std::string(argv[2]);

        if (!fileManager.readConfig(configFilename)) {
            return 1;
        }
        return 0;
    } else {
        std::unordered_map<std::string, std::string> argMap;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg.size() > 2 && arg.substr(0, 2) == "--") {
                arg = arg.substr(2);
                size_t equalsPos = arg.find('=');
                std::string argName = arg.substr(0, equalsPos);
                std::string argValue;
                if (equalsPos != std::string::npos) {
                    argValue = arg.substr(equalsPos + 1);
                }
                argMap[argName] = argValue;
            }
        }

        config.initDefault();

        if (argMap.count("trajectory")) {
            config.trajectoryFilename = argMap["trajectory"];
        } else {
            std::cout << "Trajectory file is mandatory." << std::endl;
            std::cout << "Try 'local_search --help' for more information." << std::endl;
            return 1;
        }
        if (argMap.count("time-limit")) {
            config.timeLimitMinutes = parseValue<double>(argMap["time-limit"]);
        }
        if (argMap.count("omp-threads")) {
            config.ompThreadsPerCore = parseValue<double>(argMap["omp-threads"]);
        }
        if (argMap.count("write-as-csv")) {
            config.writeAsCSV = parseBoolean(argMap["write-as-csv"]);
        }
        if (argMap.count("repetitions")) {
            config.runRepetitions = parseValue<int>(argMap["repetitions"]);
        }
        if (argMap.count("jump-chance")) {
            config.jumpFromLocalAreaChance = parseValue<double>(argMap["jump-chance"]);
        }
        if (argMap.count("random-frame-chance")) {
            config.randomFrameWhileSwappingChance = parseValue<double>(argMap["random-frame-chance"]);
        }
        if (argMap.count("memory-size")) {
            config.memorySize = parseValue<double>(argMap["memory-size"]);
        }
        if (argMap.count("random-seed")) {
            config.randomSeed = parseBoolean(argMap["random-seed"]);
        }
        if (argMap.count("matrix-size")) {
            config.matrixSize = parseValue<int>(argMap["matrix-size"]);
        }
        if (argMap.count("show-logs")) {
            config.showLogs = parseBoolean(argMap["show-logs"]);
        }
        if (argMap.count("show-rmsd-counter")) {
            config.showRMSDCounter = parseBoolean(argMap["show-rmsd-counter"]);
        }
        if (argMap.count("show-current-best")) {
            config.showDebugCurrentBest = parseBoolean(argMap["show-current-best"]);
        }
        if (argMap.count("show-route-best")) {
            config.showDebugRouteBest = parseBoolean(argMap["show-route-best"]);
        }

        DEBUG = config.showLogs;
        DEBUG_RMSD = config.showRMSDCounter;

        if (DEBUG) {
            for (const auto &kv : argMap) {
                std::cout << "(arg) [" << kv.first << "]: [" << kv.second << "]" << std::endl;
            }
        }
        return 0;
    }

    return 1;
}

int main(int argc, char *argv[]) {

    omp_init_lock(&memoryMutex);
    FileManager fileManager;
    int result = readArgs(argc, argv, fileManager);
    if (result != 0) {
        return result;
    }

    result = fileManager.readTrajectory();
    if (result != 0) {
        return result;
    }

    MEMORY_SIZE = config.matrixSize * config.matrixSize * config.memorySize;

    if (config.randomSeed) {
        srand((unsigned)time(NULL));
    } else {
        srand((unsigned)NULL);
    }

    for (int i = 0; i < config.runRepetitions; i++) {
        LocalSearch localSearch;
        resetGlobals();
        if (i == 0) {
            config.print();
        }
        localSearch.run();
    }

    omp_destroy_lock(&memoryMutex);
}
