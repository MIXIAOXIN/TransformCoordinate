#include <iostream>
#include "Trajectory.h"

int main() {
    std::string input_file = "/data/mxx/data/guangxi1/1.pos";
    std::string output_file = "/data/mxx/data/guangxi1/1-3.txt";

    utility::Trajectory* traj_io = new utility::Trajectory;
    traj_io->readPosTFile_guangxi(input_file);
    traj_io->writeTraj_XYZT(output_file);
    delete traj_io;

    return 0;
}
