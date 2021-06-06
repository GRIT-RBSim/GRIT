/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "body.hpp"
#include "loguru.hpp"
#include "math_func.hpp"
#include "utils.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <random>
#include <vector>

namespace rb_sim {

// Struct of output format.
struct OutputFormat {
  bool orbital_elements = true;
  bool position = true;
  bool velocity = true;
  bool axis = true;
  bool spin = true;
  bool obliquity = true;
  bool axial_tilt = true;
  bool Hamiltonian = true;
  bool momentum = true;
  bool customize = true;
};

enum Coordinates { CENTRAL, BARYCENTRIC, JACOBI };
enum Scheme { DEFAULT, M42_, M642_, T2_, T4_, T6_, K2_, RK4_ };

inline std::string NameOfScheme(Scheme s) {
  switch (s) {
  case DEFAULT:
    return "default";
  case M42_:
    return "M42";
  case M642_:
    return "M642";
  case T2_:
    return "T2";
  case T4_:
    return "T4";
  case T6_:
    return "T6";
  case K2_:
    return "K2";
  case RK4_:
    return "RK4";
  }
  return "default";
}

inline std::string NameOfCoordinates(Coordinates c) {
  switch (c) {
  case CENTRAL:
    return "central";
  case BARYCENTRIC:
    return "barycentric";
  case JACOBI:
    return "Jacobi";
  }
  return "central";
}

// Struct of parameters.
struct Parameters {
  // Step size of the numerical integration.
  ld step_size = 1e-3;
  // Save the system every `save_interval` time.
  ld save_interval = 1;
  // Output data every `output_gap` steps.
  uint64_t output_gaps = 100;
  // Include tidal force if `tide` is true.
  bool tide = true;
  // Include GR force if `GR` is true.
  bool GR = true;
  // The order of potential considered for the rigid body.
  uint8_t potential_order = 2;
  // Coordinates in calculating the orbital elements.
  Coordinates coordinates = CENTRAL;
  // Numerical scheme.
  Scheme scheme = DEFAULT;
};

class System {
private:
  // Directory of the system.
  std::string sys_dir;
  // true if the system is initialized from "current_system.json"
  bool is_continued_sys;

  // Specify whether compensated summation is used for q, p dynamics 
  bool compensated_summation;

  // Vector of all the bodies in the system.
  std::vector<Body> bodies;

  // total mass of all the bodies
  ld tot_mass;

  // Map the name of each body to their indices in `bodies`.
  std::map<std::string, int> name_to_id;

  // Current time of the numerical integration.
  ld current_time;
  // Target time of the numerical integration.
  ld target_time;

  // Name of the system, e.g. solar_system, Earth_Moon_system.
  std::string system_name;
  // Parameters for the integration.
  Parameters parameters;
  OutputFormat output_format;
  // Tital pairs need to be considered.
  std::set<std::pair<std::string, std::string>> tidal_pairs;

  // whether some external forces (tidal forces, GR) are specified
  bool external_force;
  // whether all the bodies are symmetric
  bool symmetric;

  std::vector<Matrix<ld>> RIR;

  // Add a body `body` to the system
  void AddBody(Body body);

  // Simualte `steps` using the Verlet integrator with order h^2.
  void T2(int steps);
  // Simualte `steps` using the Suzuki Triple Jump for Verlet with order h^4.
  void T4(int steps);
  // Simualte `steps` using a 6th order method.
  void T6(int steps);
  // Simulate `steps` using the default integrator of order h^4 + eps h^2.
  void M42(int steps);
  // Simulate `steps` using the integrator of order h^6 + eps h^4 + eps^2 h^2.
  void M642(int steps);
  void M642_S6(ld h, bool is_middle);
  void M642_Verlet(ld h);
  // Keperian splitting
  void K2(int steps);

  // traditional splitting
  void ABSchemeA(ld h);
  void ABSchemeB(ld h);
  // Tailored splitting: 4 parts
  void ABCDSchemeA(ld h);
  void ABCDSchemeB(ld h);
  void ABCDSchemeC(ld h);
  void ABCDSchemeD(ld h);
  // Tailored splitting: Keplerian
  // K1: HKeplerian
  // K2: ABCDSchemeC and HSun
  // K3: ABCDSchemeD and HInteraction
  void KABCSchemeA(ld h);
  void KABCSchemeB(ld h);
  void KABCSchemeC(ld h);
  // non symmetric H
  void SchemeE(ld h);


  // Keplerian splitting: H = H_Kep + H_Sun + H_interaction (Duncan et. al. 1998)
  void HKeplerian(ld h);
  void HSun(ld h);
  // The interaction among all the bodies except the central body with index 0
  void HInteraction(ld h);

  void UpdateRIR();

  // Simulate `steps` according to the Equations of Motion using RK4.
  void RK4(int steps);

  Vec3<ld> PVPq(unsigned int i,
                bool V2 = false) const;  // $\frac{\partial{V}}{\partial{q}}$
  Vec3<ld> PVKPq(unsigned int i) const;  // $\frac{\partial{V}}{\partial{q}}$
  Matrix<ld> PVPR(unsigned int i) const; // $\frac{\partial{V}}{\partial{R}}$
  Vec3<ld> PVPq4thOrder(unsigned int i) const;
  Matrix<ld> PVPR4thOrder(unsigned int i) const;
  // V1: point mass potential
  Vec3<ld>
  PV1Pq(unsigned int i) const; // $\frac{\partial{V1}}{\partial{q}}$
  // V2: potential from the rigid part
  Vec3<ld> PV2Pq(unsigned int i) const;   // $\frac{\partial{V2}}{\partial{q}}$
  Matrix<ld> PV2PR(unsigned int i) const; // $\frac{\partial{V2}}{\partial{R}}$

  void ClearSystem();

  // Init files in "data_in_mat/" as an empty file, then output information
  // to those files from "data_tmp.json" if "data_tmp.json" exists.
  void InitMatFile() const;
  void AppendMatFile(const std::vector<json> &j) const;
  void OutputBodyData(const json &j, std::vector<std::ofstream> &fouts) const;
  void OutputHamiltonian(const json &j, std::ofstream &fout) const;
  void OutputMomentum(const json &j, std::ofstream &fout) const;
  void OutputCustomizeData(const json &j, std::ofstream &fout) const;

  // Sanity check: serialization and deserialization.
  void SanityCheck();

  // Numerical integrate `steps` steps.
  void Simulate(int steps);
  void StepwiseOperations();

  // From input json file to `parameters`.
  void Json2Parameters(const json &j);
  // From `parameters` to output json file.
  json Parameters2Json() const;

  // From input json file to `output_format`.
  void Json2OutputFormat(const json &j);
  // From `output_format` to output json file.
  json OutputFormat2Json() const;

  // From input json file to `tidal_pairs`.
  void Json2TidalPairs(const json &j);
  // From `tidal_pairs` to output json file.
  json TidalPairs2Json() const;

  // Print the system summary to the screen.
  void PrintSystemSummary() const;

  static constexpr const OutputFormat kDefaultOutputFormat{.orbital_elements =
                                                               true,
                                                           .position = false,
                                                           .velocity = false,
                                                           .axis = false,
                                                           .spin = false,
                                                           .obliquity = true,
                                                           .axial_tilt = false,
                                                           .Hamiltonian = false,
                                                           .momentum = false,
                                                           .customize = false};

  static constexpr const Parameters kDefaultParameters{.step_size = 1e-3,
                                                       .save_interval = 1.0,
                                                       .output_gaps = 100,
                                                       .tide = false,
                                                       .GR = false,
                                                       .potential_order = 2,
                                                       .coordinates = CENTRAL,
                                                       .scheme = DEFAULT};

  // Calculating the general relativity force.
  Vec3<ld> GRForce(int id_planet) const;
  // Apply the general relativity.
  void EnforceGR(ld h);

public:
  // Return the number of bodies in the system.
  auto GetNum() const { return bodies.size(); }
  // Return the name of the `i`th body.
  auto GetName(int i) const { return bodies.at(i).name; }
  // Return the mass of the `i`th body.
  auto GetMass(int i) const { return bodies.at(i).mass; }
  // Return the current of the system.
  auto GetCurrentTime() const { return current_time; }
  // Return the step size of the integration.
  ld GetStepSize() const { return parameters.step_size; }
  // Return the output gaps of the system.
  int GetOutputGaps() const { return parameters.output_gaps; }
  // Return the save interval of the system.
  ld GetSaveInterval() const { return parameters.save_interval; }
  // Return the `sys_dir`.
  std::string GetDir() const { return sys_dir; }
  // Return current time.
  ld GetTime() const { return current_time; }

  // Initialize the system from files in `dir`
  void InitFromFile(const std::string &dir);
  // Initialize the system according to the json.
  void InitSystem(const json &j);
  // Set the target time of the numerical simulation as `t`.
  void SetTargetTime(const ld t);
  // Set output directory to be `dir`.
  void SetOutputDir(const std::string &dir);

  void Normalize();
  void Simulate();

  // Return the Hamiltonian of the system.
  ld Hamiltonian() const;
  // Return the potential energy of the system.
  ld PotentialEnergy() const;
  // Return the angular momentum of the system.
  Vec3<ld> AngularMomentum() const;
  // Return the linear momentum of the system.
  Vec3<ld> LinearMomentum() const;

  // Get the axial tilt angle of the `idx`th body.
  ld GetAxialTilt(int idx) const;
  // Get the Eular angle of the `idx`th body.
  Matrix<ld> GetEulerAngles(int idx) const;
  // Get the position of the `idx`th body.
  Vec3<ld> GetPos(int idx) const;
  // Get the velocity of the `idx`th body.
  Vec3<ld> GetVel(int idx) const;

  // Equations of the motion of the system: for the RK4 integrator.
  void EoM(std::vector<Vec3<ld>> &kq, std::vector<Vec3<ld>> &kp,
           std::vector<Vec3<ld>> &kpi, std::vector<Matrix<ld>> &kR);

  // Return the tidal acceleration.
  // host: rigid body, guest: point mass.
  // `id_host`: the id of the host body.
  // `id_guest`: the id of the guest body.
  Vec3<ld> TidalForceEquilibrium(int id_host, int id_guest) const;
  // Apply the tidal force at each time step of the numerical integration.
  void EnforceTide(ld h);

  // Serialization.
  std::ostream &operator>>(std::ostream &os) const;
  // Deserialization.
  std::istream &operator<<(std::istream &os);

  void Dump(std::ostream &&os) const;
  void DumpDataTmp(std::vector<json> &data_json) const;
  json Dump() const;
  void OutputData(std::ostream &os) const;
  json OutputDataJson() const;
  void StoreDataInMat(const std::string &dir, const json &data_json) const;
  void OutputInitSysPosVel(const std::string &dir) const;
  };

inline std::istream &operator>>(std::istream &is, System &sys) {
  return sys << is;
}

inline std::ostream &operator<<(std::ostream &os, const System &sys) {
  return sys >> os;
}

} // namespace rb_sim

std::vector<ld> CustomizeData(const rb_sim::System &sys);

#endif
