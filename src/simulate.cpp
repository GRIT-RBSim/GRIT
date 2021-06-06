/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#define LOGURU_IMPLEMENTATION 1
#include "loguru.hpp"
#include "system.hpp"
#include "utils.hpp"

using namespace std;
using namespace rb_sim;

int main(int argc, char *argv[]) {

  // Check the command "./simulate dir".
  assert(argc == 2 && "Usage: ./simulate dir");
  string dir = string(argv[1]) + "/";

  System system;
  // Initialize the system from the files in the directory `dir`.
  system.InitFromFile(dir);
  // Store a copy of init_system as a new file.
  system.OutputInitSysPosVel(dir);

  // Search for files: "target_time.json".
  auto tt_file = GetFileContents((dir + "/target_time.json").c_str());
  auto tt_json = json::parse(tt_file);
  // Get the target time of simulation.
  ld target_time = tt_json["target_time"];

  // Set the target time of simulation.
  system.SetTargetTime(target_time);
  // Simulate.
  system.Simulate();

  return 0;
}
