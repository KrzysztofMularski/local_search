#ifndef PROGRESS_H
#define PROGRESS_H

#include <iostream>
#include "globals.h"

class Progress {
private:
    int barWidth;
    int currentSteps;
    int allStepsCount;

public:
    Progress(int stepsCount) {
        allStepsCount = stepsCount;
        barWidth = 70;
        currentSteps = 0;
        if (!DEBUG) {
            return;
        }
        std::cout << "[>";
        for (int i = 1; i < barWidth; ++i) {
            std::cout << " ";
        }
        printf("] 0 %%\r");
        std::cout.flush();
    }

    void improve() {
        if (!DEBUG) {
            return;
        }
        if (currentSteps % 100000) {
            currentSteps++;
            return;
        }
        std::cout << "[";
        int pos = barWidth * float(currentSteps)/float(allStepsCount);
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        printf("] %.2f%%\r", float(currentSteps)/float(allStepsCount)*100.0);
        std::cout.flush();
        currentSteps++;
    }

    void end() {
        if (!DEBUG) {
            return;
        }
        std::cout << "[";
        int pos = barWidth * 1;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] 100%      ";
        std::cout.flush();
        std::cout << std::endl;
        std::cout << "Done!" << std::endl;
    }
};

#endif // PROGRESS_H
