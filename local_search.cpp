#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <iomanip>

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

// Numer of atoms in [<sphere>]
// std::vector<int> sphereSize;

// List of atoms in [<sphere>]
// std::vector<std::vector<int>> sphereAtoms;
std::vector<std::vector<int>> *sphereAtoms;

Config config;

std::unordered_set<std::pair<int, int>, PairHash> memorySet;

omp_lock_t memoryMutex;

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

    inline void swappingAllocations(int &allocatedOnFrame, int &changingFrame) {
        int temp = changingFrame;
        changingFrame = allocatedOnFrame;
        allocatedOnFrame = temp;
        // debug("[swapping] (", changingFrame, ", ", allocatedOnFrame, ") -> (",
        // allocatedOnFrame, ", ", changingFrame, ")");

        rmsd.atomsAllocation(allocatedOnFrame);
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

    LocalSearchResult traverse_old(int i, int j, bool allocationOnInit = true) {
        // traversing route: i != j checked
        int allocatedOnFrame = i;
        int changingFrame = j;
        if (allocationOnInit) {
            rmsd.atomsAllocation(allocatedOnFrame);
        }
        LocalSearchResult routeBest = {
            rmsd.calculateRMSDSuperpose(changingFrame),
            allocatedOnFrame,
            changingFrame,
        };
        // debug("[_New route] (", routeBest.i, ", ", routeBest.j, ") = ",
        // routeBest.rmsdValue); debug("(", routeBest.i, ", ", routeBest.j, ")=",
        // routeBest.rmsdValue);
        bool changedSidesAlready = false;
        bool swappedAlready = false;

        // choosing direction
        int step = getRandom(0, 1) * 2 - 1; // -1 or +1

        // debug("[_1] step=", step, ", changedSidesAlready=", changedSidesAlready,
        // ", swappedAlready=", swappedAlready);

        while (true) {
            // traversing route part - straight line

            int newChangingFrame = changingFrame + step;
            // debug("[_Cheking cell] (", allocatedOnFrame, ", ", newChangingFrame,
            // ")"); debug("(", allocatedOnFrame, ", ", newChangingFrame, ")");
            // debug("[_2] newChangingFrame=", newChangingFrame);
            if (!insideMatrixBoundaries(newChangingFrame) || newChangingFrame == allocatedOnFrame) {
                // debug("[_Failed]");
                // debug("[_3] !insideMatrixBoundaries=",
                // !insideMatrixBoundaries(newChangingFrame)); debug("[_4] frames has
                // same id: ", newChangingFrame == allocatedOnFrame); try changing sides
                if (!changedSidesAlready) {
                    // debug("[_5] !changedSidesAlready ", !changedSidesAlready);
                    step *= -1;
                    changedSidesAlready = true;
                    continue;
                } else if (!swappedAlready) {
                    // debug("[_6] !swappedAlready ", !swappedAlready);

                    swappedAlready = true;
                    step = getRandom(0, 1) * 2 - 1; // -1 or +1
                    changedSidesAlready = false;
                    swappingAllocations(allocatedOnFrame, changingFrame);
                    continue;
                } else {
                    // debug("[_7] !swappedAlready ", !swappedAlready);
                    return routeBest;
                }
            }

            double newValue = rmsd.calculateRMSDSuperpose(newChangingFrame);

            // debug("[_New value] (", allocatedOnFrame, ", ", newChangingFrame, ") =
            // ", newValue); debug("(", allocatedOnFrame, ", ", newChangingFrame, ") =
            // ", newValue); debug("[_8] newValue=", newValue);
            if (newValue > routeBest.rmsdValue) {
                saveIfRouteBest(routeBest, newValue, allocatedOnFrame, newChangingFrame);
                changingFrame = newChangingFrame;
                while (true) {
                    // following direction
                    newChangingFrame = changingFrame + step;
                    // debug("[_9] newChangingFrame=", newChangingFrame);
                    // debug("[_Cheking cell] (", allocatedOnFrame, ", ",
                    // newChangingFrame, ")"); debug("(", allocatedOnFrame, ", ",
                    // newChangingFrame, ") ");
                    if (!insideMatrixBoundaries(newChangingFrame) || newChangingFrame == allocatedOnFrame) {
                        // debug("[_Failed] ", "!inside=",
                        // !insideMatrixBoundaries(newChangingFrame), ", ",
                        //       newChangingFrame == allocatedOnFrame);
                        // debug("[_10]");
                        break;
                    }
                    newValue = rmsd.calculateRMSDSuperpose(newChangingFrame);
                    // debug("[_New value] (", allocatedOnFrame, ", ", newChangingFrame,
                    // ") = ", newValue); debug("(", allocatedOnFrame, ", ",
                    // newChangingFrame, ") = ", newValue); debug("[_11] newValue=",
                    // newValue);
                    if (newValue > routeBest.rmsdValue) {
                        // debug("[_12]");
                        saveIfRouteBest(routeBest, newValue, allocatedOnFrame, newChangingFrame);
                        changingFrame = newChangingFrame;
                        // debug("[_13] changingFrame=", changingFrame);
                        continue;
                    } else {
                        // debug("[_14]");
                        break;
                    }
                }
                // cannot change direction because we came from it, so we swap
                swappedAlready = true;
                changedSidesAlready = false;
                swappingAllocations(allocatedOnFrame, changingFrame);
                step = getRandom(0, 1) * 2 - 1; // -1 or +1
                // debug("[_15] step=", step);
                continue;
            } else if (!changedSidesAlready) {
                step *= -1;
                changedSidesAlready = true;
                // debug("[_17] step=", step);
                continue;
            } else {
                // debug("[_18]");
                if (!swappedAlready) {
                    // debug("[_19]");

                    swappedAlready = true;
                    changedSidesAlready = false;
                    swappingAllocations(allocatedOnFrame, changingFrame);
                    step = getRandom(0, 1) * 2 - 1; // -1 or +1
                    continue;
                } else {
                    // debug("[_20]");
                    return routeBest;
                }
            }
            auto stop = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = stop - start;
            if (elapsed.count() > config.timeLimitMinutes * 60) {
                return routeBest;
            }
        }
    }

    // void calculateBestPairs() {
    //     for (const std::pair<int, int> &pair : bestResultsArray) {
    //         rmsd.atomsAllocation(pair.first);
    //         double newValue = rmsd.calculateRMSDSuperpose(pair.second);
    //         debug("[No sqrt -> sqrt]: [", pair.first, ", ", pair.second, "] = ", newValue);
    //     }
    // }

    // void testSpecificPairs() {
    //     // [DEBUG] [Current best]: [803, 220] = 15.3397  94.911
    //     // [DEBUG] [Current best]: [887, 106] = 18.3842  113.509
    //     // [DEBUG] [Current best]: [872, 115] = 18.5634  114.105
    //     // [DEBUG] [Current best]: [965, 57] = 21.369    118.398
    //     // [DEBUG] [Current best]: [23, 949] = 22.6254   124.4
    //     //
    //     // [DEBUG] [Current best]: [0, 934] = 21.5337
    //     // [DEBUG] [Current best]: [995, 40] = 21.9388
    //     // [DEBUG] [Current best]: [45, 986] = 22.179653
    //
    //     const std::vector<std::pair<int, int>> arr = {
    //         { 803, 220 }, { 887, 106 }, { 872, 115 }, { 965, 57 }, { 23, 949 }, { 0, 934 }, { 995, 40 }, { 45, 986 },
    //     };
    //
    //     for (const std::pair<int, int> &pair : arr) {
    //         rmsd.atomsAllocation(pair.first);
    //         double newValue = rmsd.calculateRMSDSuperpose(pair.second);
    //         debug("[custom]: [", pair.first, ", ", pair.second, "] = ", newValue);
    //     }
    // }

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

        // przez 5(parametr) procent przez ileś iteracji(parametr) nie poprawia
        // zawsze sfery są z niższego indeksu
        // wylosować 1 i 2 klatki
        // pierwsza zostaje 1, druga podróżuje
        // potem druga zostaje i 1 podróżuje
        // potem losujemy inną parę
        //
        // while(1) {
        // czy koniec będzie zależny od czasu od ostatniego rozwiązania, od
        // różnic rmsd, czy jako parametr czy przejmować się już znalezionymi
        // wartościami czy uznać że rmsd pomiędzy 1 a 2 jest takie samo jak
        // pomiędzy 2 a 1?
        // }

        auto stop = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = stop - start;
        print("Local Search Results:");
        print(" - Computation time: ", elapsed.count(), "s");
        print(" - RMSD counted: ", RMSDCalculationCountGlobal, " times.");
        print(" - Atoms allocated: ", AllocationsCountGlobal, " times.");
        return;
    }
};

int main(int argc, char *argv[]) {

    omp_init_lock(&memoryMutex);
    FileManager fileManager;
    // if (argc <= 1) {
    //     std::cout << "Too few arguments" << std::endl;
    //     return 1;
    // }

    // std::string configFilename = argv[1];
    std::string configFilename = "./config.yml";
    if (!fileManager.readConfig(configFilename)) {
        return 1;
    }

    fileManager.readTrajectory();

    if (config.randomSeed) {
        srand((unsigned)time(NULL));
    } else {
        srand((unsigned)NULL);
    }

    LocalSearch localSearch;

    config.print();

    localSearch.run();
    omp_destroy_lock(&memoryMutex);
}
