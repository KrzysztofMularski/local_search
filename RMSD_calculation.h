#ifndef RMSD_CALCULATION_H
#define RMSD_CALCULATION_H

#include <cmath>
#include <eigen3/Eigen/Geometry>

#include "globals.h"

class RMSDCalculation {
  private:

    // calculating distance between 2 atoms, used when allocating atoms to spheres.
    double atomsDistanceCalc(int atom1, int atom2) {

        double dx = A[FRAMEONE][atom1][0] - A[FRAMEONE][atom2][0];
        double dy = A[FRAMEONE][atom1][1] - A[FRAMEONE][atom2][1];
        double dz = A[FRAMEONE][atom1][2] - A[FRAMEONE][atom2][2];

        double result = dx*dx + dy*dy + dz*dz;

        return sqrt(result);
    }

    // Find3DAffineTransform is from oleg-alexandrov repository on github, available here
    // https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp [as of 27.01.2022]
    // Given two sets of 3D points, find the rotation + translation + scale
    // which best maps the first set to the second.
    // Source: http://en.wikipedia.org/wiki/Kabsch_algorithm

    // The input 3D points are stored as columns.
    Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out) {

        // Default output
        Eigen::Affine3d A;
        A.linear() = Eigen::Matrix3d::Identity(3, 3);
        A.translation() = Eigen::Vector3d::Zero();

        if (in.cols() != out.cols())
            throw "Find3DAffineTransform(): input data mis-match";

        // First find the scale, by finding the ratio of sums of some distances,
        // then bring the datasets to the same scale.
        double dist_in = 0, dist_out = 0;
        for (int col = 0; col < in.cols() - 1; col++) {
            dist_in += (in.col(col + 1) - in.col(col)).norm();
            dist_out += (out.col(col + 1) - out.col(col)).norm();
        }
        if (dist_in <= 0 || dist_out <= 0)
            return A;
        double scale = dist_out / dist_in;
        out /= scale;

        // Find the centroids then shift to the origin
        Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
        Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
        for (int col = 0; col < in.cols(); col++) {
            in_ctr += in.col(col);
            out_ctr += out.col(col);
        }
        in_ctr /= in.cols();
        out_ctr /= out.cols();
        for (int col = 0; col < in.cols(); col++) {
            in.col(col) -= in_ctr;
            out.col(col) -= out_ctr;
        }

        // SVD
        Eigen::MatrixXd Cov = in * out.transpose();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // Find the rotation
        double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
        if (d > 0)
            d = 1.0;
        else
            d = -1.0;
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
        I(2, 2) = d;
        Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

        // The final transform
        A.linear() = scale * R;
        A.translation() = scale * (out_ctr - R * in_ctr);

        return A;
    }

    // superpose changes std::vector of atoms frame 2 to map atoms from frame 1 in the way to
    // minimise RMSD between both frames
    void superpose(const std::vector<std::vector<double>> &frame1, std::vector<std::vector<double>> &frame2) {
        int atomsInSphere = frame1.size();
        Eigen::Matrix3Xd S1(3, atomsInSphere);
        Eigen::Matrix3Xd S2(3, atomsInSphere);
        for (int j = 0; j < atomsInSphere; j++) {
            for (int k = 0; k < 3; k++) {
                S1(k, j) = frame1[j][k];
                S2(k, j) = frame2[j][k];
            }
        }
        Eigen::Affine3d RT = Find3DAffineTransform(S2, S1);
        S2 = RT.linear() * S2;
        for (int j = 0; j < atomsInSphere; j++) {
            S2.block<3, 1>(0, j) += RT.translation();
        }
        for (int j = 0; j < atomsInSphere; j++) {
            for (int k = 0; k < 3; k++) {
                frame2[j][k] = S2(k, j);
            }
        }
    }

    bool pairInMemory(int f1, int f2) {
        omp_set_lock(&memoryMutex);
        std::pair<int, int> newPair = std::make_pair(f1, f2);
        if (memorySet.find(newPair) == memorySet.end()) {
            // not in memory
            memorySet.insert(newPair);
            if (memorySet.size() > MEMORY_SIZE) {
                memorySet.erase(memorySet.begin());
            }
            omp_unset_lock(&memoryMutex);
            return false;
        } else {
            // in memory
            omp_unset_lock(&memoryMutex);
            return true;
        }
    }

  public:
    // calculating RMSD on spheres, on choosen frames
    double calculateRMSDSuperpose(int secondFrame) {
        if (config.memorySize != 0 && pairInMemory(FRAMEONE, FRAMETWO)) {
            return -1.0;
        }
        // else calculate rmsd
        FRAMETWO = secondFrame;
        RMSDCalculationCount++;
        debugRMSD();
        double result = 0;
        std::vector<std::vector<std::vector<double>>> sphereMatrix;
        for (int s = 0; s < SPHERES; s++) {
            int atomsInSphere = sphereAtoms[omp_thread_id][s].size();
            sphereMatrix.assign(2, {});
            sphereMatrix[0].assign(atomsInSphere, {});
            sphereMatrix[1].assign(atomsInSphere, {});

            for (int j = 0; j < atomsInSphere; j++) {
                sphereMatrix[0][j] = A[FRAMEONE][sphereAtoms[omp_thread_id][s][j]];
                sphereMatrix[1][j] = A[FRAMETWO][sphereAtoms[omp_thread_id][s][j]];
            }
            double tempResult = 0;
            superpose(sphereMatrix[0], sphereMatrix[1]);
            for (int j = 0; j < atomsInSphere; j++) {
                for (int k = 0; k < 3; k++ ) {
                    double tempRMSD = sphereMatrix[1][j][k] - sphereMatrix[0][j][k];
                    tempResult += tempRMSD * tempRMSD;
                }
            }
            tempResult /= atomsInSphere * 3.0;
            tempResult = sqrt(tempResult);
            result += tempResult;
        }
        return result;
    }

    // allocating atoms into spheres, based on sphereRadius
    void atomsAllocation(int firstFrame) {
        FRAMEONE = firstFrame;
        AllocationsCount++;
        std::vector<int> temp;
        temp.reserve(ATOMS);
        sphereAtoms[omp_thread_id].assign(SPHERES, temp);
        for (int i = 0; i < ATOMS; i++) {
            for (int j = 0; j < SPHERES; j++) {
                if (atomsDistanceCalc(i, sphereCA[j]) <= sphereRadius) {
                    sphereAtoms[omp_thread_id][j].push_back(i);
                }
            }
        }
    }
};

#endif // RMSD_CALCULATION_H
