/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#include "system.hpp"
#include <omp.h>

namespace rb_sim {

constexpr const OutputFormat System::kDefaultOutputFormat;

constexpr const Parameters System::kDefaultParameters;

void System::InitFromFile(const std::string &dir) {

  sys_dir = dir;
  // whether Normalize the system such that the linear momentum is (0,0,0).
  bool normalize = false;

  // Search for file as the initialization of the system: "current_system.json".
  auto init_file = dir + "/current_system.json";

  // If "current_system.json" exists, set it as `init_file`, otherwise set
  // "init_system.json" as `init_file`.
  if (!IsReadable(init_file)) {
    // If the "current_system.json" doesn't exist, search for "init_system.json"
    // and erase the "data_tmp.json" file if it exists.
    auto data_tmp_file = dir + "/data_tmp.json";
    if (std::ifstream(data_tmp_file)) {
      remove(data_tmp_file.c_str());
    }
    init_file = dir + "/init_system.json";
    normalize = true;
    is_continued_sys = false;
  } else
    is_continued_sys = true;
  // If none of "current_system.json", "init_system.json" exist. Error.
  if (!IsReadable(init_file)) {
    std::cerr << "file not found" << std::endl;
    exit(1);
  }
  std::cout << "System initialized from " << init_file << "\n";

  // If the directory "data_in_mat" doesn't exist, create one.
  if (!loguru::create_directories((sys_dir + "/data_in_mat/").c_str())) {
    std::cerr << "create " << sys_dir << "/data_in_mat failed" << std::endl;
    return;
  }

  auto f = GetFileContents(init_file.c_str());
  json j;

  try {
    j = json::parse(f);
  } catch (json::parse_error & ex) {
    std::cerr << "Invalid json format. Check your input file." << ex.byte << std::endl;
  }

  InitSystem(j);
  if (normalize) Normalize();

  InitMatFile();

  // SanityCheck();
}

void System::ClearSystem() {
  bodies.clear();
  name_to_id.clear();
  tidal_pairs.clear();
  RIR.clear();
}

void System::InitSystem(const json &j) {

  compensated_summation = false;

  ClearSystem();

  assert(j.is_object() && "Not a valid json object. Check your json format.");

  // Get the name of the system.
  if (j.contains("system_name"))
    system_name = j["system_name"].get<std::string>();
  else {
    system_name = "a_poor_system_with_no_name";
    std::cerr << "Warning: You didn't provide a name for the system, it is set "
                 "to be the default name.\n";
  }

  // Get the current_time of the system.
  if (j.contains("current_time"))
    current_time = j["current_time"];
  else {
    current_time = 0;
    std::cerr << "Warning: You didn't provide a valid current_time, it is set "
                 "to be 0.0 as default.\n";
  }

  // Set the parameters to be the default one first.
  parameters = kDefaultParameters;
  assert(j.contains("parameters") &&
         "Need parameters in your system file. See \"sample/init_system.json\" "
         "for example.");
  Json2Parameters(j["parameters"]);

  // Add data of all the bodies to the system.
  assert(j.contains("body") &&
         "Need information of bodies: add json object array in \"body\".");
  auto jbody = j["body"];
  symmetric = true;
  for (json::iterator it = jbody.begin(); it != jbody.end(); it++) {
    Body new_body(it.value());
    AddBody(new_body);
  }

  // Set the output format: specify the values included in the output files.
  output_format = kDefaultOutputFormat;
  if (j.contains("output_format")) {
    Json2OutputFormat(j["output_format"]);
  } else {
    std::cerr
        << "Warning: Your output format is set to be the default one. To "
           "customize your output, add \"output_format\" in your input file\n";
  }

  // Set the tidal pairs.
  if (parameters.tide) {
    if (j.contains("tidal_pairs")) {
      Json2TidalPairs(j["tidal_pairs"]);
    } else {
      for (unsigned int i = 0; i < bodies.size(); i++)
        for (unsigned int j = 0; j < bodies.size(); j++)
          if (i != j && bodies[i].rigid) {
            assert(bodies[i].CheckTidalParameters() &&
                   "Need time_lag and Q to evaulate tidal forces.");
            tidal_pairs.insert(std::make_pair(GetName(i), GetName(j)));
            std::cerr
                << "Warning: You didn't specify the tidal pairs. It is set to "
                   "be the default one: all the rigid objects are affected by "
                   "the tidal forces raised by all other objects.\n";
          }
    }
  }

  external_force = parameters.tide || parameters.GR;

  PrintSystemSummary();

  // Initialize `target_time` to be the `current_time`.
  target_time = current_time;

  // Calculate totol mass
  tot_mass = 0;
  for (unsigned int i = 0; i < bodies.size(); i++)
    tot_mass += bodies[i].mass;

  // Initialize the RIR vector which stores R(t)*I*R(t)^T for each object.
  RIR.resize(bodies.size());
  UpdateRIR();
}

void System::Json2Parameters(const json &j) {
  // Parameters that must be provided by users.
  assert(j.contains("step_size") && "Need step size.");
  parameters.step_size = j["step_size"];
  assert(j.contains("save_interval") && "Need save_interval.");
  parameters.save_interval = j["save_interval"];
  assert(j.contains("output_gaps") && "Need output_gaps.");
  parameters.output_gaps = j["output_gaps"];
  // Tide.
  if (j.contains("tide"))
    parameters.tide = j["tide"];
  // GR.
  if (j.contains("GR"))
    parameters.GR = j["GR"];
  // Coordinates
  if (j.contains("coordinates")) {
    std::string cds = j["coordinates"];
    if (cds == "barycentric")
      parameters.coordinates = BARYCENTRIC;
    else if (cds == "Jacobi")
      parameters.coordinates = JACOBI;
    else if (cds == "central")
      parameters.coordinates = CENTRAL;
    else {
      parameters.coordinates = kDefaultParameters.coordinates;
      std::cerr
          << "Warning: Your coordinates is not supported. The coordinates "
             "is set to be "
          << kDefaultParameters.coordinates << " coordinates.\n";
    }
  }

  // Potential order.
  if (j.contains("potential_order")) {
    parameters.potential_order = j["potential_order"];
    assert(
        (parameters.potential_order == 2 || parameters.potential_order == 4) &&
        "Invalid potential order.");
  }

  // Scheme.
  if (j.contains("scheme")) {
    if (j["scheme"] == "default")
      parameters.scheme = DEFAULT;
    else if (j["scheme"] == "T2")
      parameters.scheme = T2_;
    else if (j["scheme"] == "T4")
      parameters.scheme = T4_;
    else if (j["scheme"] == "T6")
      parameters.scheme = T6_;
    else if (j["scheme"] == "M42")
      parameters.scheme = M42_;
    else if (j["scheme"] == "M642")
      parameters.scheme = M642_;
    else if (j["scheme"] == "RK4")
      parameters.scheme = RK4_;
    else if (j["scheme"] == "K2")
      parameters.scheme = K2_;
    else
      std::cerr << "Your scheme is invalid. It is set as the default scheme.\n";
  }
}

json System::Parameters2Json() const {
  json j;
  j["step_size"] = parameters.step_size;
  j["save_interval"] = parameters.save_interval;
  j["output_gaps"] = parameters.output_gaps;
  j["tide"] = parameters.tide;
  j["GR"] = parameters.GR;
  j["potential_order"] = parameters.potential_order;
  j["coordinates"] = NameOfCoordinates(parameters.coordinates);
  j["scheme"] = NameOfScheme(parameters.scheme);
  return j;
}

void System::Json2OutputFormat(const json &j) {
  if (j.contains("orbital_elements"))
    output_format.orbital_elements = j["orbital_elements"];
  if (j.contains("position"))
    output_format.position = j["position"];
  if (j.contains("velocity"))
    output_format.velocity = j["velocity"];
  if (j.contains("axis"))
    output_format.axis = j["axis"];
  if (j.contains("spin"))
    output_format.spin = j["spin"];
  if (j.contains("obliquity"))
    output_format.obliquity = j["obliquity"];
  if (j.contains("axial_tilt"))
    output_format.axial_tilt = j["axial_tilt"];
  if (j.contains("Ftf"))
      output_format.Ftf = j["Ftf"];
  if (j.contains("Eheat"))
      output_format.Eheat = j["Eheat"];
  if (j.contains("Hamiltonian"))
    output_format.Hamiltonian = j["Hamiltonian"];
  if (j.contains("momentum"))
    output_format.momentum = j["momentum"];
  if (j.contains("customize"))
    output_format.customize = j["customize"];
}

json System::OutputFormat2Json() const {
  json j;
  j["orbital_elements"] = output_format.orbital_elements;
  j["position"] = output_format.position;
  j["velocity"] = output_format.velocity;
  j["axis"] = output_format.axis;
  j["spin"] = output_format.spin;
  j["obliquity"] = output_format.obliquity;
  j["axial_tilt"] = output_format.axial_tilt;
  j["Ftf"] = output_format.Ftf;
  j["Eheat"] = output_format.Eheat;
  j["Hamiltonian"] = output_format.Hamiltonian;
  j["momentum"] = output_format.momentum;
  j["customize"] = output_format.customize;
  return j;
}

void System::Json2TidalPairs(const json &j) {
  for (auto it = j.begin(); it != j.end(); ++it) {
    std::vector<std::string> name_list = it.value();
    int idx = name_to_id[it.key()];
    bodies[idx].CheckTidalParameters();
    for (auto name : name_list) {
      tidal_pairs.insert(std::make_pair(it.key(), name));
    }
  }
}

json System::TidalPairs2Json() const {
  json j;
  for (unsigned int i = 0; i < bodies.size(); i++)
    j[GetName(i)] = json::array();
  for (auto tp : tidal_pairs)
    j[tp.first].push_back(tp.second);
  return j;
}

void System::EnforceGR(ld h) {

  for (auto i = 0u; i < bodies.size(); i++) {
    auto tmp = GRForce(i);
    auto mu = bodies[i].mass;
    auto F_GR = tmp;
    bodies[i].p = bodies[i].p + F_GR * h * mu;
  }
}

Vec3<ld> System::GRForce(int id_planet) const {

  const int id_star = 0;
  if (id_planet == id_star)
    return {0, 0, 0};
  ld cconst = 173.1445988 * kDaysOfYear; // in unit of AU/yr
  const ld m_star = bodies[id_star].mass;
  const auto m_planet = bodies[id_planet].mass;
  const auto rad_vec = bodies[id_planet].q - bodies[id_star].q;
  const auto vel_vec =
      bodies[id_planet].p / m_planet - bodies[id_star].p / m_star;
  const auto r = rad_vec.Norm();
  const auto r2 = r * r, r3 = r2 * r;
  const auto vm2 = vel_vec.DotProduct(vel_vec);
  const ld vr = rad_vec.DotProduct(vel_vec);
  const ld cons = ld(4.0) * kG * m_star / r - vm2;
  const ld kk = kG * m_star / (r3 * cconst * cconst);
  auto F_GR = (rad_vec * cons + vel_vec * vr * 4.0) * kk;

  return F_GR;
}

void System::PrintSystemSummary() const {
  std::cout << std::setfill('=') << std::setw(50) << "System Summary"
            << std::setw(40) << "=\n";

  std::cout << "# System name: " << system_name << "\n";
  std::cout << "# Bodies in the system:"
            << "\n";

  std::cout << "# " << std::setfill(' ') << std::setw(32) << "name"
            << std::setw(16) << "rigidity\n";
  for (unsigned int i = 0; i < bodies.size(); i++) {
    std::cout << "# " << std::setfill(' ') << std::setw(10) << i
              << "th body: " << std::setw(15) << GetName(i);
    if (bodies[i].rigid)
      std::cout << std::setw(15) << "rigid";
    else
      std::cout << std::setw(15) << "point mass";
    std::cout << "\n";
  }
  std::cout << "# Coordinates of orbital elements: "
            << NameOfCoordinates(parameters.coordinates) << "\n";

  std::cout << "\n";
  std::cout << std::setfill('=') << std::setw(54) << "Numerical Integrator"
            << std::setw(36) << "=\n";
  std::cout << "# Scheme: " << NameOfScheme(parameters.scheme) << "\n";
  std::cout << "# Step size: " << GetStepSize() << " year \n";
  std::cout << "# System will be saved every " << GetSaveInterval()
            << " year in \"current_system.json\"\n";
  std::cout << "# Data output scale: " << GetOutputGaps() * GetStepSize()
            << " year\n";
  std::cout << "# Tidal effects: " << (parameters.tide ? "" : "not ")
            << "considered\n";
  std::cout << "# General relativity: " << (parameters.GR ? "" : "not ")
            << "considered\n";

  std::cout << "# Output variables: ";
  auto j = OutputFormat2Json();
  for (auto it = j.begin(); it != j.end(); it++)
    if (it.value())
      std::cout << it.key() << " ";
  std::cout << "\n";
}

void System::SetTargetTime(const ld t) { target_time = t; }

void System::SetOutputDir(const std::string &dir) { sys_dir = dir; }

void System::SanityCheck() {

  // Dump the current system to the `j`.
  auto j = Dump();

  System x, y;

  std::string dir_x = sys_dir + "/.tmp_x";
  std::string dir_y = sys_dir + "/.tmp_y";
  if (!loguru::create_directories((dir_x).c_str())) {
    std::cerr << "create " << dir_x << "failed" << std::endl;
  }
  if (!loguru::create_directories((dir_y).c_str())) {
    std::cerr << "create " << dir_y << "failed" << std::endl;
  }

  // Evolve x for 100 steps.
  x.InitSystem(j);
  x.SetTargetTime(100 * GetStepSize());
  x.SetOutputDir(dir_x);
  x.Simulate();

  // Evolve y for 50 steps + 50 steps.
  y.InitSystem(j);
  y.SetOutputDir(dir_y);
  y.SetTargetTime(50 * GetStepSize());
  y.Simulate();
  j = y.Dump();
  y.InitSystem(j);
  y.SetTargetTime(100 * GetStepSize());
  y.Simulate();

  // Compare two systems.
  assert(x.Dump() == y.Dump() && "serialization && deserialization failed");
}

void System::Normalize() {
  auto lm = LinearMomentum();
  ld mass_sum = 0;
  for (unsigned int i = 0; i < bodies.size(); i++) {
    mass_sum += bodies[i].mass;
  }

  for (unsigned int i = 0; i < bodies.size(); i++) {
    bodies[i].p -= lm * bodies[i].mass / mass_sum;
  }
}

json System::Dump() const {
  json j = {{"body", json::array()}};
  j["system_name"] = system_name;
  j["current_time"] = current_time;
  for (unsigned int idx = 0; idx < bodies.size(); idx++) {
    j["body"].push_back(bodies[idx]);
  }
  j["parameters"] = Parameters2Json();
  j["tidal_pairs"] = TidalPairs2Json();
  j["output_format"] = OutputFormat2Json();
  return j;
}

void System::DumpDataTmp(std::vector<json> &data_json) const {
  auto fout = std::ofstream(sys_dir + "/data_tmp.json", std::ios::app);
  for (auto &d : data_json) {
    fout << d.dump() << std::endl;
  }
  AppendMatFile(data_json);
  data_json.clear();
  fout.close();
}

void System::Dump(std::ostream &&os) const {
  os << Dump().dump(4) << std::endl;
}

void System::AddBody(Body body) {
  bodies.push_back(body);
  symmetric &= (body.symmetric || (!body.rigid));
}

void System::Simulate() {

  // Fetch step size `h`, simulation target time `target_time` and the gaps for
  // saving data `gaps`.
  const ld h = GetStepSize();
  const int gaps = GetOutputGaps();

  // Initialize the `last_save_time` to be the current time for further counting
  // the next saving time.
  auto last_save_time = current_time;

  // `data_json` is the vector storing temp file.
  std::vector<json> data_json;

  // Timer.
  ld tot_time = target_time - current_time;
  auto start = std::chrono::high_resolution_clock::now();
  auto prev = start;
  std::chrono::duration<double> elapsed_in_s;

  // Iterate over all time steps of the integration.
  std::cout << "\n"
            << std::setfill('=') << std::setw(48) << "Simulation"
            << std::setw(42) << "\n";
  std::cout << "# From t0=" << std::to_string(current_time)
            << " to T=" << std::to_string(target_time) + ": \n\n";
  ProgressBar(std::cout, 0, 0, current_time);
  if (!is_continued_sys)
    data_json.push_back(OutputDataJson());
  while (current_time < target_time) {
    // Print progress bar.
    elapsed_in_s = std::chrono::high_resolution_clock::now() - start;
    std::chrono::duration<double> time_gap =
        std::chrono::high_resolution_clock::now() - prev;
    if (double(time_gap.count()) >= 0.2) {
      ProgressBar(std::cout, float(1 - (target_time - current_time) / tot_time),
                  elapsed_in_s.count(), current_time);
      prev = std::chrono::high_resolution_clock::now();
    }
    // simulate
    if (current_time + h * gaps < target_time)
      Simulate(gaps);
    else
      while (current_time < target_time)
        Simulate(1);

    // Append output json data of current step to data_json.
    data_json.push_back(OutputDataJson());
    // Output file after 1000 steps in case it takes too much memory.
    // To store less frequently, the threshold can be set larger, but must be
    // greater than the steps in SanityCheck().
    if (data_json.size() > 1) //GL change
      DumpDataTmp(data_json);

    // Once hit the save_interval set by the user. Save the current system in
    // the file "current_system.json".
    if (current_time - last_save_time + h >= GetSaveInterval()) {
      last_save_time = current_time;
      Dump(std::ofstream(sys_dir + "/current_system.json"));
    }
  }
  // The final progress bar.
  elapsed_in_s = std::chrono::high_resolution_clock::now() - start;
  ProgressBar(std::cout, 1, elapsed_in_s.count(), current_time);
  std::cout << "Finished.\n";
  // Dump rest data.
  DumpDataTmp(data_json);
  // Dump the final system after numerical simulation.
  Dump(std::ofstream(sys_dir + "/current_system.json"));
}

void System::Simulate(int steps) {
  switch (parameters.scheme) {
  case DEFAULT:
    M42(steps);
    break;
  case M42_:
    M42(steps);
    break;
  case M642_:
    M642(steps);
    break;
  case T2_:
    T2(steps);
    break;
  case T4_:
    T4(steps);
    break;
  case T6_:
    T6(steps);
    break;
  case K2_:
    K2(steps);
    break;
  case RK4_:
    RK4(steps);
    break;
  }
}

void System::StepwiseOperations() {
  if (parameters.tide)
    EnforceTide(GetStepSize());
  if (parameters.GR)
    EnforceGR(GetStepSize());
}

void System::T2(int steps) {
  const ld h = GetStepSize();
  for (int current_step = 0; current_step < steps; current_step++) {

    current_time += h;

    if (current_step == 0)
      ABSchemeA(h / ld(2.0));
    else
      ABSchemeA(h);

    if (!external_force && symmetric)
      ABSchemeB(h);
    else if (external_force && symmetric) {
      ABSchemeB(h / ld(2.0));
      StepwiseOperations();
      ABSchemeB(h / ld(2.0));
    } else if (!external_force && !symmetric) {
      ABSchemeB(h / ld(2.0));
      SchemeE(h);
      ABSchemeB(h / ld(2.0));
    } else {
      ABSchemeB(h / ld(2.0));
      SchemeE(h / ld(2.0));
      StepwiseOperations();
      SchemeE(h / ld(2.0));
      ABSchemeB(h / ld(2.0));
    }
  }
  ABSchemeA(h / ld(2.0));
}

void System::T4(int steps) {

  // order 4 Suzuki Triple Jump
  // verlet: phi, then phi(gamma3 * h)•phi(gamma2 * h)•phi(gamma1 * h)

  const ld h = GetStepSize();
  ld gamma1 = ld(1.0) / (ld(2.0) - pow(ld(2.0), ld(1.0) / ld(3.0)));
  ld gamma2 = ld(1.0) - ld(2.0) * gamma1;
  for (int current_step = 0; current_step < steps; current_step++) {

    current_time += h;

    if (current_step == 0)
      ABSchemeA(h * gamma1 / ld(2.0));
    else
      ABSchemeA(h * gamma1);
    ABSchemeB(h * gamma1);
    ABSchemeA(h * (gamma1 + gamma2) / ld(2.0));
    if (!external_force && symmetric)
      ABSchemeB(h * gamma2);
    else if (external_force && symmetric) {
      ABSchemeB(h * gamma2 / ld(2.0));
      StepwiseOperations();
      ABSchemeB(h * gamma2 / ld(2.0));
    } else if (!external_force && !symmetric) {
      ABSchemeB(h * gamma2 / ld(2.0));
      SchemeE(h);
      ABSchemeB(h * gamma2 / ld(2.0));
    } else {
      ABSchemeB(h * gamma2 / ld(2.0));
      SchemeE(h / ld(2.0));
      StepwiseOperations();
      SchemeE(h / ld(2.0));
      ABSchemeB(h * gamma2 / ld(2.0));
    }
    ABSchemeA(h * (gamma2 + gamma1) / ld(2.0));
    ABSchemeB(h * gamma1);
  }
  ABSchemeA(h * gamma1 / ld(2.0));
}

void System::T6(int steps) {

  const ld h = GetStepSize();

  // Yoshida 1990 order 6 Solution A
  const ld a1 = h * 0.784513610477560;
  const ld a2 = h * 0.235573213359357;
  const ld a3 = h * (-0.117767998417887e1);
  const ld a4 = h * 1.31518632068390628476;

  for (int current_step = 0; current_step < steps; current_step++) {

    current_time += h;

    if (current_step == 0)
      ABSchemeA(a1 / 2);
    else
      ABSchemeA(a1);
    ABSchemeB(a1);
    ABSchemeA((a1 + a2) / 2);
    ABSchemeB(a2);
    ABSchemeA((a2 + a3) / 2);
    ABSchemeB(a3);
    ABSchemeA((a3 + a4) / 2);
    ABSchemeB(a4 / 2);
    SchemeE(h / ld(2.0));
    StepwiseOperations();
    SchemeE(h / ld(2.0));
    ABSchemeB(a4 / 2);
    ABSchemeA((a4 + a3) / 2);
    ABSchemeB(a3);
    ABSchemeA((a3 + a2) / 2);
    ABSchemeB(a2);
    ABSchemeA((a2 + a1) / 2);
    ABSchemeB(a1);
  }
  ABSchemeA(a1 / 2);
}

void System::M42(int steps) {
  const ld h = GetStepSize();
  ld gamma1 = ld(1.0) / (ld(2.0) - pow(ld(2.0), ld(1.0) / ld(3.0)));
  ld gamma2 = ld(1.0) - ld(2.0) * gamma1;
  for (int current_step = 0; current_step < steps; current_step++) {
    current_time += h;

    if (current_step == 0) {
      ABCDSchemeA(h * gamma1 / ld(2.0));
    } else
      ABCDSchemeA(h * gamma1);
    ABCDSchemeB(h * gamma1);
    ABCDSchemeA(h * (gamma1 + gamma2) / ld(2.0));
    ABCDSchemeB(h * gamma2 / ld(2.0));

    ABCDSchemeC(h / ld(2.0));
    if (!external_force && symmetric)
      ABCDSchemeD(h);
    else if (external_force && symmetric) {
      ABCDSchemeD(h / ld(2.0));
      StepwiseOperations();
      ABCDSchemeD(h / ld(2.0));
    } else if (!external_force && !symmetric) {
      ABCDSchemeD(h / ld(2.0));
      SchemeE(h);
      ABCDSchemeD(h / ld(2.0));
    } else {
      ABCDSchemeD(h / ld(2.0));
      SchemeE(h / ld(2.0));
      StepwiseOperations();
      SchemeE(h / ld(2.0));
      ABCDSchemeD(h / ld(2.0));
    }
    ABCDSchemeC(h / ld(2.0));

    ABCDSchemeB(h * gamma2 / ld(2.0));
    ABCDSchemeA(h * (gamma2 + gamma1) / ld(2.0));
    ABCDSchemeB(h * gamma1);
  }
  ABCDSchemeA(h * gamma1 / ld(2.0));
}

void System::M642(int steps) {

  const ld h = GetStepSize();
  for (int current_step = 0; current_step < steps; current_step++) {

    current_time += h;

    // ABA42
    M642_S6((3 - sqrt(3)) * h / 6, false);
    M642_Verlet(h / 2);
    M642_S6(h / sqrt(3), true);
    M642_Verlet(h / 2);
    M642_S6((3 - sqrt(3)) * h / 6, false);
  }
}

void System::M642_S6(ld h, bool is_middle_step) {

  // Yoshida 1990
  const ld a1 = h * 0.784513610477560;
  const ld a2 = h * 0.235573213359357;
  const ld a3 = h * (-0.117767998417887e1);
  const ld a4 = h * 1.31518632068390628476;

  ABCDSchemeA(a1 / 2);
  ABCDSchemeB(a1);
  ABCDSchemeA((a1 + a2) / 2);
  ABCDSchemeB(a2);
  ABCDSchemeA((a2 + a3) / 2);
  ABCDSchemeB(a3);
  ABCDSchemeA((a3 + a4) / 2);
  if (is_middle_step) {
    ABCDSchemeB(a4 / 2);
    SchemeE(h / ld(2.0));
    StepwiseOperations();
    SchemeE(h / ld(2.0));
    ABCDSchemeB(a4 / 2);
  } else
    ABCDSchemeB(a4);
  ABCDSchemeA((a4 + a3) / 2);
  ABCDSchemeB(a3);
  ABCDSchemeA((a3 + a2) / 2);
  ABCDSchemeB(a2);
  ABCDSchemeA((a2 + a1) / 2);
  ABCDSchemeB(a1);
  ABCDSchemeA(a1 / 2);
}

void System::M642_Verlet(ld h) {
  ABCDSchemeC(h / ld(2.0));
  ABCDSchemeD(h);
  ABCDSchemeC(h / ld(2.0));
}

void System::K2(int steps) {
  const ld h = GetStepSize();
  for (int current_step = 0; current_step < steps; current_step++) {

    current_time += h;

    if (current_step == 0)
      KABCSchemeC(h / 2.0);
    else
      KABCSchemeC(h);

    KABCSchemeB(h / 2.0);

    if (!external_force && symmetric) {
      KABCSchemeA(h);
    } else if (external_force && symmetric) {
      KABCSchemeA(h / 2.0);
      StepwiseOperations();
      KABCSchemeA(h / 2.0);
    } else if (!external_force && !symmetric) {
      KABCSchemeA(h / 2.0);
      SchemeE(h);
      KABCSchemeA(h / 2.0);
    } else {
      KABCSchemeA(h / 2.0);
      SchemeE(h / 2.0);
      StepwiseOperations();
      SchemeE(h / 2.0);
      KABCSchemeA(h / 2.0);
    }

    KABCSchemeB(h / 2.0);
  }
  KABCSchemeC(h / 2.0);
}

void System::ABSchemeA(ld h) {

  for (auto i = 0u; i < bodies.size(); i++) {
    if (compensated_summation) {
      auto y = bodies[i].p * (h / bodies[i].mass) - bodies[i].qc;
      auto t = bodies[i].q + y;
      bodies[i].qc = (t - bodies[i].q) - y;
      bodies[i].q = t;
    } else
      bodies[i].q += bodies[i].p * (h / bodies[i].mass);
  }

#pragma omp parallel for if (bodies.size()>50)
  for (auto i = 0u; i < bodies.size(); i++) {
    if (bodies[i].rigid) {
      ld tmp = bodies[i].pi.Norm();
      ld theta = h * (ld(1.0) / bodies[i].I3 - ld(1.0) / bodies[i].I2) *
                 bodies[i].pi.GetElement<3>();
      ld gamma = h * tmp / bodies[i].I2;

      auto Rz = RotZ(theta);
      auto Rpi = RotAngle2Mat(bodies[i].pi, gamma);
      
      bodies[i].pi = Rz.Transpose() * bodies[i].pi;
      bodies[i].R = bodies[i].R * Rpi * Rz;
    }
  }

  UpdateRIR();
}

void System::ABSchemeB(ld h) {
// scheme B
#pragma omp parallel for if (bodies.size()>50)
  for (auto i = 0u; i < bodies.size(); i++) {
    if (compensated_summation) {
      auto y = PVPq(i) * (-h) - bodies[i].pc;
      auto t = bodies[i].p + y;
      bodies[i].pc = (t - bodies[i].p) - y;
      bodies[i].p = t;
    } else
      bodies[i].p -= PVPq(i) * h;
    if (bodies[i].rigid)
      bodies[i].pi -= Rot(bodies[i].R.Transpose() * PVPR(i)) * h;
  }
}

void System::ABCDSchemeA(ld h) {
#pragma omp parallel for if (bodies.size()>50)
  for (auto i = 0u; i < bodies.size(); i++) {
    auto add = bodies[i].p * (h / bodies[i].mass);
    if (compensated_summation) {
      auto y = add - bodies[i].qc;
      auto t = bodies[i].q + y;
      bodies[i].qc = (t - bodies[i].q) - y;
      bodies[i].q = t;
    } else
      bodies[i].q += add;
  }
}

void System::ABCDSchemeB(ld h) {
#pragma omp parallel for if (bodies.size()>50)
  for (auto i = 0u; i < bodies.size(); i++) {
    auto add = PV1Pq(i) * (-h);
    if (compensated_summation) {
      auto y = add - bodies[i].pc;
      auto t = bodies[i].p + y;
      bodies[i].pc = (t - bodies[i].p) - y;
      bodies[i].p = t;
    } else
      bodies[i].p += add;
  }
}

void System::ABCDSchemeC(ld h) {
#pragma omp parallel for if (bodies.size()>50)
  for (auto i = 0u; i < bodies.size(); i++) {
    if (bodies[i].rigid) {
      ld tmp = bodies[i].pi.Norm();
      ld theta = h * (ld(1.0) / bodies[i].I3 - ld(1.0) / bodies[i].I2) *
                 bodies[i].pi.GetElement<3>();
      ld gamma = h * tmp / bodies[i].I2;

      auto Rz= RotAngle2Mat(kE3, theta);
      auto Rpi = RotAngle2Mat(bodies[i].pi, gamma);
      
      bodies[i].pi = Rz.Transpose() * bodies[i].pi;
      bodies[i].R = bodies[i].R * Rpi * Rz;
    }
  }

  UpdateRIR();
}

void System::ABCDSchemeD(ld h) {
#pragma omp parallel for if (bodies.size()>50)
  for (auto i = 0u; i < bodies.size(); i++) {
    if (compensated_summation) {
      auto y = PV2Pq(i) * (-h) - bodies[i].pc;
      auto t = bodies[i].p + y;
      bodies[i].pc = (t - bodies[i].p) - y;
      bodies[i].p = t;
    } else
      bodies[i].p -= PV2Pq(i) * h;

    if (bodies[i].rigid)
      bodies[i].pi -= Rot(bodies[i].R.Transpose() * PV2PR(i)) * h;
  }
}

void System::SchemeE(ld h) {
  for (auto i = 0u; i < bodies.size(); i++)
    if (bodies[i].rigid)
      if (!bodies[i].symmetric) {
        ld theta = bodies[i].delta * bodies[i].pi[1] * h;
        auto rot = RotY(theta);
        bodies[i].R = rot * bodies[i].R;
        bodies[i].pi = rot.Transpose() * bodies[i].pi;
      }
}

void System::KABCSchemeA(ld h) { HKeplerian(h); }

void System::KABCSchemeB(ld h) {
  HSun(h);
  ABCDSchemeC(h);
}

void System::KABCSchemeC(ld h) {
  HInteraction(h);
  ABCDSchemeD(h);
}

void System::HKeplerian(ld h) {

  // switch to democratic heliocentric variables Q, P
  // Duncan et. al. 1998
  std::vector<Vec3<ld>> Q;
  std::vector<Vec3<ld>> P;
  Q.push_back(bodies[0].q * bodies[0].mass);
  P.push_back(bodies[0].p);
  for (unsigned int i = 1; i < bodies.size(); i++) {
    Q.push_back(bodies[i].q - bodies[0].q);
    Q[0] = Q[0] + bodies[i].q * bodies[i].mass;
    P[0] = P[0] + bodies[i].p;
  }
  Q[0] = Q[0] / tot_mass;
  for (unsigned int i = 1; i < bodies.size(); i++) {
    P.push_back(bodies[i].p - P[0] * bodies[i].mass / tot_mass);
  }

  // Drift here
  for (unsigned int i = 1; i < bodies.size(); i++) {
    auto p = Q[i];
    auto v = P[i] / bodies[i].mass;
    KeplerSlvDanby(h, kG * GetMass(0), p, v);
    Q[i] = p;
    P[i] = v * bodies[i].mass;
  }

  // switch back to Cartesian coordinates in the inertia frame
  bodies[0].q = Q[0];
  for (unsigned int i = 1; i < bodies.size(); i++) {
    bodies[0].q -= Q[i] * bodies[i].mass / tot_mass;
  }
  for (unsigned int i = 1; i < bodies.size(); i++) {
    bodies[i].q = Q[i] + bodies[0].q;
  }
  bodies[0].p = P[0];
  for (unsigned int i = 1; i < bodies.size(); i++) {
    bodies[i].p = P[i] + P[0] * bodies[i].mass / tot_mass;
    bodies[0].p -= bodies[i].p;
  }
}

void System::HSun(ld h) {
  auto tmp = bodies[0].p;
  ld m0 = bodies[0].mass;
  for (auto i = 0u; i < bodies.size(); i++) {
    tmp -= bodies[i].p * (m0 / tot_mass);
  }
  for (auto i = 0u; i < bodies.size(); i++) {
    auto add = tmp * h;
    if (i == 0)
      add *= (1 - m0 / tot_mass) / m0;
    else
      add *= -1 / tot_mass;
    if (compensated_summation) {
      auto y = add - bodies[i].qc;
      auto t = bodies[i].q + y;
      bodies[i].qc = (t - bodies[i].q) - y;
      bodies[i].q = t;
    } else
      bodies[i].q += add;
  }
}

void System::HInteraction(ld h) {
#pragma omp parallel for if (bodies.size()>50)
  for (auto i = 1u; i < bodies.size(); i++) {
    auto add = PVKPq(i) * (-h);
    if (compensated_summation) {
      auto y = add - bodies[i].pc;
      auto t = bodies[i].p + y;
      bodies[i].pc = (t - bodies[i].p) - y;
      bodies[i].p = t;
    } else
      bodies[i].p += add;
  }
}

void System::UpdateRIR() {
  for (auto i = 0u; i < bodies.size(); i++) {
    if (bodies[i].rigid) {
      RIR[i] = bodies[i].R * bodies[i].I * bodies[i].R.Transpose();
    }
  }
}

Vec3<ld> System::PVPq(unsigned int i, bool V2) const {
  Vec3<ld> ret(ld(0.0));
  // potential order equals 2
  if (bodies[i].rigid) {
    for (unsigned int j = 0u; j < bodies.size(); j++) {
      if (j != i) {
        Vec3<ld> qi_qj = bodies[i].q - bodies[j].q;
        ld div2 = qi_qj.NormSquared();
        ld div = sqrt(div2);
        ld div3 = div2 * div;
        ld div5 = div3 * div2;
        ld div7 = div5 * div2;

        ld tmp = 0.0;
        if (!V2) tmp = bodies[i].mass * bodies[j].mass / div3;
        if (bodies[j].rigid) {
          tmp += ld(1.5) / div5 *
                 (bodies[j].mass * bodies[i].traceI +
                  bodies[i].mass * bodies[j].traceI);
          tmp -= qi_qj.DotProduct(
                     (bodies[j].mass * RIR[i] + bodies[i].mass * RIR[j]) *
                     (qi_qj)) *
                 (ld(7.5) / div7);
          ret += (bodies[j].mass * RIR[i] + bodies[i].mass * RIR[j]) * qi_qj *
                 (ld(3.0) / div5);
        } else {
          tmp += ld(1.5) / div5 * (bodies[j].mass * bodies[i].traceI);
          tmp -= qi_qj.DotProduct(bodies[j].mass * RIR[i] * (qi_qj)) *
                 (ld(7.5) / div7);
          ret += bodies[j].mass * RIR[i] * qi_qj * (ld(3.0) / div5);
        }
        ret += qi_qj * tmp;
      }
    }
  } else {
    for (auto j = 0u; j < bodies.size(); j++) {
      if (j != i) {
        Vec3<ld> qi_qj = bodies[i].q - bodies[j].q;
        ld div2 = qi_qj.NormSquared();
        ld div = sqrt(div2);
        ld div3 = div2 * div;

        ld tmp = 0.0;
        if (!V2) tmp = bodies[i].mass * bodies[j].mass / div3;

        if (bodies[j].rigid) {
          ld div5 = div3 * div2;
          ld div7 = div5 * div2;
          tmp += ld(1.5) / div5 * (bodies[i].mass * bodies[j].traceI);
          tmp -= qi_qj.DotProduct(bodies[i].mass * RIR[j] * (qi_qj)) *
                 (ld(7.5) / div7);
          ret += bodies[i].mass * RIR[j] * qi_qj * (ld(3.0) / div5);
        }
        ret += qi_qj * tmp;
      }
    }
  }
  ret = ret * kG;
  if (parameters.potential_order == 4) {
    ret += PVPq4thOrder(i);
  }
  return ret;
}

Vec3<ld> System::PV1Pq(unsigned int i) const {
  Vec3<ld> ret(ld(0.0));
  for (auto j = 0u; j < bodies.size(); j++) {
    if (j != i) {
      Vec3<ld> qi_qj = bodies[i].q - bodies[j].q;
      ld div2 = qi_qj.NormSquared();
      ld div = sqrt(div2);
      ld div3 = div2 * div;
      ld tmp = bodies[i].mass * bodies[j].mass / div3;

      ret += qi_qj * tmp;
    }
  }
  ret = ret * kG;
  return ret;
}

Vec3<ld> System::PV2Pq(unsigned int i) const {
  return PVPq(i, true);
}

Vec3<ld> System::PVKPq(unsigned int i) const {
  Vec3<ld> ret(ld(0.0));
  // potential order equals 2
    for (auto j = 1u; j < bodies.size(); j++) {
      if (j != i) {
        Vec3<ld> qi_qj = bodies[i].q - bodies[j].q;
        ld div2 = qi_qj.NormSquared();
        ld div = sqrt(div2);
        ld div3 = div2 * div;

        ret += qi_qj * (bodies[i].mass * bodies[j].mass / div3);
      }
  }
  ret = ret * kG;
  return ret;
}

Vec3<ld> System::PVPq4thOrder(unsigned int i) const {

  Vec3<ld> ret(ld(0.0));

  for (auto j = 0u; j < bodies.size(); j++) {
    if (j != i && (bodies[i].rigid || bodies[j].rigid)) {
      Vec3<ld> qi_qj = bodies[i].q - bodies[j].q;
      ld div2 = qi_qj.NormSquared();
      ld div = sqrt(div2);
      ld div3 = div2 * div;
      ld div5 = div3 * div2;
      ld div7 = div5 * div2;
      ld div9 = div7 * div2;
      ld div11 = div9 * div2;

      ld rai = 0, rbi = 0, rci = 0, rai2 = 0, rbi2 = 0, rci2 = 0;
      if (bodies[i].rigid) {
        bodies[i].GetSemiAxes(rai, rbi, rci);
        rai2 = rai * rai;
        rbi2 = rbi * rbi;
        rci2 = rci * rci;
      }
      ld raj = 0, rbj = 0, rcj = 0, raj2 = 0, rbj2 = 0, rcj2 = 0;
      if (bodies[j].rigid) {
        bodies[j].GetSemiAxes(raj, rbj, rcj);
        raj2 = raj * raj;
        rbj2 = rbj * rbj;
        rcj2 = rcj * rcj;
      }

      auto Ri = kIdentity;
      if (bodies[i].rigid)
        Ri = bodies[i].R;
      auto Rj = kIdentity;
      if (bodies[j].rigid)
        Rj = bodies[j].R;
      ld x = qi_qj[0], y = qi_qj[1], z = qi_qj[2];
      Matrix<ld> qi_qj2 =
          Matrix<ld>(std::vector<ld>{x * x, x * y, x * z, y * x, y * y, y * z,
                                     z * x, z * y, z * z},
                     3, 3);

      ld tmp = 0;
      tmp += 15.0 / 8 / div7 *
             (3 * rai2 * rai2 + 3 * rbi2 * rbi2 + 3 * rci2 * rci2 +
              2 * (rai2 * rbi2 + rai2 * rci2 + rbi2 * rci2)) /
             35;
      tmp += 15.0 / 8 / div7 *
             (3 * raj2 * raj2 + 3 * rbj2 * rbj2 + 3 * rcj2 * rcj2 +
              2 * (raj2 * rbj2 + raj2 * rcj2 + rbj2 * rcj2)) /
             35;
      auto tmp_mat = Ri.Transpose() * Rj *
                     Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5, 0,
                                                0, 0, rcj2 / 5},
                                3, 3) *
                     Rj.Transpose() * Ri *
                     Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5, 0,
                                                0, 0, rci2 / 5},
                                3, 3);
      tmp += 15.0 / 4 / div7 * tmp_mat.Trace();
      tmp_mat =
          (Ri *
           Matrix<ld>(std::vector<ld>{rai2 * (3 * rai2 + rbi2 + rci2), 0, 0, 0,
                                      rbi2 * (rai2 + 3 * rbi2 + rci2), 0, 0, 0,
                                      rci2 * (rai2 + rbi2 + 3 * rci2)},
                      3, 3) *
           Ri.Transpose()) /
          35;
      tmp -= 15.0 * 7 / 4 / div9 * (qi_qj.DotProduct(tmp_mat * (qi_qj)));
      ret +=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (15.0 / 2 / div7);
      tmp_mat =
          Rj *
          Matrix<ld>(std::vector<ld>{raj2 * (3 * raj2 + rbj2 + rcj2), 0, 0, 0,
                                     rbj2 * (raj2 + 3 * rbj2 + rcj2), 0, 0, 0,
                                     rcj2 * (raj2 + rbj2 + 3 * rcj2)},
                     3, 3) *
          Rj.Transpose() / 35;
      tmp -= 15.0 * 7 / 4 / div9 * (qi_qj.DotProduct(tmp_mat * (qi_qj)));
      ret +=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (15.0 / 2 / div7);

      tmp_mat =
          (rai2 + rbi2 + rci2) / 5 * Rj *
              Matrix<ld>(std::vector<ld>{raj2, 0, 0, 0, rbj2, 0, 0, 0, rcj2}, 3,
                         3) /
              5 * Rj.Transpose() +
          (raj2 + rbj2 + rcj2) / 5 * Ri *
              Matrix<ld>(std::vector<ld>{rai2, 0, 0, 0, rbi2, 0, 0, 0, rci2}, 3,
                         3) /
              5 * Ri.Transpose();
      tmp -= 15.0 * 7 / 4 / div9 * (qi_qj.DotProduct(tmp_mat * qi_qj));
      ret +=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (15.0 / 2 / div7);

      tmp_mat = Ri *
                Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5, 0, 0, 0,
                                           rci2 / 5},
                           3, 3) *
                Ri.Transpose() * Rj *
                Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5, 0, 0, 0,
                                           rcj2 / 5},
                           3, 3) *
                Rj.Transpose();
      tmp -= 15 * 7 / div9 * (qi_qj.DotProduct(tmp_mat * qi_qj));
      ret +=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (15 * 2 / div7);

      tmp_mat = qi_qj2 * Ri *
                Matrix<ld>(std::vector<ld>{rai2 * rai2 * 3 / 35, 0, 0, 0,
                                           rbi2 * rbi2 * 3 / 35, 0, 0, 0,
                                           rci2 * rci2 * 3 / 35},
                           3, 3) *
                Ri.Transpose() * qi_qj2;
      tmp += 35.0 * 9 / 8 / div11 * tmp_mat.Trace();
      tmp_mat = qi_qj2 * Ri *
                Matrix<ld>(std::vector<ld>{rai2 * rai2 * 3 / 35, 0, 0, 0,
                                           rbi2 * rbi2 * 3 / 35, 0, 0, 0,
                                           rci2 * rci2 * 3 / 35},
                           3, 3) *
                Ri.Transpose();
      ret -=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (35.0 / 4 / div9);
      tmp_mat = Ri *
                Matrix<ld>(std::vector<ld>{rai2 * rai2 * 3 / 35, 0, 0, 0,
                                           rbi2 * rbi2 * 3 / 35, 0, 0, 0,
                                           rci2 * rci2 * 3 / 35},
                           3, 3) *
                Ri.Transpose() * qi_qj2;
      ret -=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (35.0 / 4 / div9);

      tmp_mat = qi_qj2 * Rj *
                Matrix<ld>(std::vector<ld>{raj2 * raj2 * 3 / 35, 0, 0, 0,
                                           rbj2 * rbj2 * 3 / 35, 0, 0, 0,
                                           rcj2 * rcj2 * 3 / 35},
                           3, 3) *
                Rj.Transpose() * qi_qj2;
      tmp += 35.0 * 9 / 8 / div11 * tmp_mat.Trace();
      tmp_mat = qi_qj2 * Rj *
                Matrix<ld>(std::vector<ld>{raj2 * raj2 * 3 / 35, 0, 0, 0,
                                           rbj2 * rbj2 * 3 / 35, 0, 0, 0,
                                           rcj2 * rcj2 * 3 / 35},
                           3, 3) *
                Rj.Transpose();
      ret -=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (35.0 / 4 / div9);
      tmp_mat = Rj *
                Matrix<ld>(std::vector<ld>{raj2 * raj2 * 3 / 35, 0, 0, 0,
                                           rbj2 * rbj2 * 3 / 35, 0, 0, 0,
                                           rcj2 * rcj2 * 3 / 35},
                           3, 3) *
                Rj.Transpose() * qi_qj2;
      ret -=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (35.0 / 4 / div9);

      tmp += 105.0 * 9 / 4 / div11 *
             qi_qj.DotProduct(
                 Ri *
                 Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5, 0, 0,
                                            0, rci2 / 5},
                            3, 3) *
                 Ri.Transpose() * qi_qj2 * Rj *
                 Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5, 0, 0,
                                            0, rcj2 / 5},
                            3, 3) *
                 Rj.Transpose() * qi_qj);
      tmp_mat = Ri *
                Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5, 0, 0, 0,
                                           rci2 / 5},
                           3, 3) *
                Ri.Transpose() * qi_qj2 * Rj *
                Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5, 0, 0, 0,
                                           rcj2 / 5},
                           3, 3) *
                Rj.Transpose();
      ret -=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (105.0 / 2 / div9);
      tmp_mat = Rj *
                Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5, 0, 0, 0,
                                           rcj2 / 5},
                           3, 3) *
                Rj.Transpose() * qi_qj2 * Ri *
                Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5, 0, 0, 0,
                                           rci2 / 5},
                           3, 3) *
                Ri.Transpose();
      ret -=
          bodies[i].mass * bodies[j].mass * tmp_mat * qi_qj * (105.0 / 2 / div9);

      tmp *= bodies[i].mass * bodies[j].mass;
      ret += qi_qj * tmp;
    }
  }
  return ret * kG;
}

Matrix<ld> System::PVPR(unsigned int i) const {
  Matrix<ld> ret(ld(0.0), 3, 3);
  if (bodies[i].rigid) {
    for (auto j = 0u; j < bodies.size(); j++) {
      if (j != i) {
        Vec3<ld> qi_qj_m = bodies[i].q - bodies[j].q;
        ld div2 = qi_qj_m.NormSquared();
        ld div = sqrt(div2);
        ld div5 = div2 * div2 * div;
        std::vector<ld> v{qi_qj_m.GetElement<1>(), qi_qj_m.GetElement<2>(),
                          qi_qj_m.GetElement<3>()};
        Matrix<ld> qi_qj(v, 3, 1);
        ret += bodies[j].mass * ld(3.0) / div5 * qi_qj * (qi_qj.Transpose());
      }
    }
    return kG * ret * bodies[i].R * bodies[i].I;
    if (parameters.potential_order == 4) {
      ret += PVPR4thOrder(i);
    }
  }
  return ret;
}

Matrix<ld> System::PVPR4thOrder(unsigned int i) const {

  Matrix<ld> ret(ld(0.0), 3, 3);

  for (auto j = 0u; j < bodies.size(); j++) {
    if (j != i && (bodies[i].rigid || bodies[j].rigid)) {
      Vec3<ld> qi_qj = bodies[i].q - bodies[j].q;
      ld div2 = qi_qj.NormSquared();
      ld div = sqrt(div2);
      ld div3 = div2 * div;
      ld div5 = div3 * div2;
      ld div7 = div5 * div2;
      ld div9 = div7 * div2;

      ld rai = 0, rbi = 0, rci = 0, rai2 = 0, rbi2 = 0, rci2 = 0;
      if (bodies[i].rigid) {
        bodies[i].GetSemiAxes(rai, rbi, rci);
        rai2 = rai * rai;
        rbi2 = rbi * rbi;
        rci2 = rci * rci;
      }
      ld raj = 0, rbj = 0, rcj = 0, raj2 = 0, rbj2 = 0, rcj2 = 0;
      if (bodies[j].rigid) {
        bodies[j].GetSemiAxes(raj, rbj, rcj);
        raj2 = raj * raj;
        rbj2 = rbj * rbj;
        rcj2 = rcj * rcj;
      }

      auto Ri = kIdentity;
      if (bodies[i].rigid)
        Ri = bodies[i].R;
      auto Rj = kIdentity;
      if (bodies[j].rigid)
        Rj = bodies[j].R;

      auto v = std::vector<ld>{qi_qj[0], qi_qj[1], qi_qj[2]};
      Matrix<ld> qi_qj_mat(v, 3, 1);
      Matrix<ld> qi_qj2 = qi_qj_mat * (qi_qj_mat.Transpose());

      auto tmp_mat =
          ld(2.0) * (Rj *
                     Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5, 0,
                                                0, 0, rcj2 / 5},
                                3, 3) *
                     Rj.Transpose() * Ri *
                     Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5, 0,
                                                0, 0, rci2 / 5},
                                3, 3));
      ret -= bodies[j].mass * tmp_mat * 3 / 4 / div5;

      tmp_mat =
          ld(2.0) * qi_qj2 * Ri *
          Matrix<ld>(std::vector<ld>{rai2 * (3 * rai2 + rbi2 + rci2), 0, 0, 0,
                                     rbi2 * (rai2 + 3 * rbi2 + rci2), 0, 0, 0,
                                     rci2 * (rai2 + rbi2 + 3 * rci2)},
                     3, 3) /
          35;
      ret += bodies[j].mass * tmp_mat * 15 / 4 / div7;

      tmp_mat = ld(2.0) * qi_qj2 * Ri * (raj2 + rbj2 + rcj2) / 5 * Ri *
                Matrix<ld>(std::vector<ld>{rai2, 0, 0, 0, rbi2, 0, 0, 0, rci2},
                           3, 3) /
                5;
      ret += bodies[j].mass * tmp_mat * 15 / 4 / div7;

      tmp_mat = ld(2.0) * qi_qj2 * Rj *
                Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5, 0, 0, 0,
                                           rcj2 / 5},
                           3, 3) *
                Rj.Transpose() * Ri *
                Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5, 0, 0, 0,
                                           rci2 / 5},
                           3, 3);
      ret += bodies[j].mass * tmp_mat * 15 / div7;

      tmp_mat = ld(2.0) * qi_qj2 * qi_qj2 * Ri *
                Matrix<ld>(std::vector<ld>{rai2 * rai2 * 3 / 35, 0, 0, 0,
                                           rbi2 * rbi2 * 3 / 35, 0, 0, 0,
                                           rci2 * rci2 * 3 / 35},
                           3, 3);
      ret -= bodies[j].mass * tmp_mat * 35 / 8 / div9;

      tmp_mat = ld(2.0) * qi_qj2 * Rj *
                Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5, 0, 0, 0,
                                           rcj2 / 5},
                           3, 3) *
                Rj.Transpose() * qi_qj2 * Ri *
                Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5, 0, 0, 0,
                                           rci2 / 5},
                           3, 3);
      ret -= bodies[j].mass * tmp_mat * 105 / 4 / div9;
    }
  }
  return kG * bodies[i].mass * ret;
}

Matrix<ld> System::PV2PR(unsigned int i) const { return PVPR(i); }

ld System::Hamiltonian() const {
  ld H = 0.0;
  for (auto i = 0u; i < bodies.size(); i++) {
    H += bodies[i].p.NormSquared() / bodies[i].mass / 2.0;
    if (bodies[i].rigid) {
      H += bodies[i].pi.DotProduct((bodies[i].Iinv * bodies[i].pi)) / 2.0;
    }
  }
  return H + PotentialEnergy();
}

ld System::PotentialEnergy() const {
  ld E = 0;
  for (auto i = 0u; i < bodies.size(); i++) {
    for (auto j = 0u; j < i; j++) {
        Vec3<ld> qi_qj = bodies[i].q - bodies[j].q;
        ld div2 = qi_qj.NormSquared();
        ld div = sqrt(div2);

        ld tmp = -bodies[i].mass * bodies[j].mass / div;
        if (bodies[i].rigid || bodies[j].rigid) {
          ld div3 = div2 * div;
          ld div5 = div3 * div2;

          Matrix<ld> Mat(ld(0.0), 3, 3);
          if (bodies[i].rigid) {
              tmp -= (bodies[i].traceI) * bodies[j].mass / div3 / 2.0;
              Mat += RIR[i] * bodies[j].mass;
          }
          if (bodies[j].rigid) {
              tmp -= (bodies[j].traceI) * bodies[i].mass / div3 / 2.0;
              Mat += RIR[j] * bodies[i].mass;
          }

          tmp += (ld)(1.5) / div5 * qi_qj.DotProduct(Mat*qi_qj);
        }
        E += tmp;
    }
  }
  if (parameters.potential_order == 4) {

    for (auto i = 0u; i < bodies.size(); i++) {
      for (auto j = 0u; j < bodies.size(); j++) {
        if (j != i && (bodies[i].rigid || bodies[j].rigid)) {
          Vec3<ld> qi_qj = bodies[i].q - bodies[j].q;
          ld div2 = qi_qj.NormSquared();
          ld div = sqrt(div2);
          ld div3 = div2 * div;
          ld div5 = div3 * div2;
          ld div7 = div5 * div2;
          ld div9 = div7 * div2;

          ld rai = 0, rbi = 0, rci = 0, rai2 = 0, rbi2 = 0, rci2 = 0;
          if (bodies[i].rigid) {
            bodies[i].GetSemiAxes(rai, rbi, rci);
            rai2 = rai * rai;
            rbi2 = rbi * rbi;
            rci2 = rci * rci;
          }
          ld raj = 0, rbj = 0, rcj = 0, raj2 = 0, rbj2 = 0, rcj2 = 0;
          if (bodies[j].rigid) {
            bodies[j].GetSemiAxes(raj, rbj, rcj);
            raj2 = raj * raj;
            rbj2 = rbj * rbj;
            rcj2 = rcj * rcj;
          }

          auto Ri = kIdentity;
          if (bodies[i].rigid)
            Ri = bodies[i].R;
          auto Rj = kIdentity;
          if (bodies[j].rigid)
            Rj = bodies[j].R;
          ld x = qi_qj[0], y = qi_qj[1], z = qi_qj[2];
          Matrix<ld> qi_qj2 =
              Matrix<ld>(std::vector<ld>{x * x, x * y, x * z, y * x, y * y,
                                         y * z, z * x, z * y, z * z},
                         3, 3);
          E -= bodies[i].mass * bodies[j].mass * 3.0 / 8.0 / div5 *
               (3 * rai2 * rai2 + 3 * rbi2 * rbi2 + 3 * rci2 * rci2 +
                2 * rai2 * rbi2 + 2 * rai2 * rci2 + 2 * rbi2 * rci2 +
                3 * raj2 * raj2 + 3 * rbj2 * rbj2 + 3 * rcj2 * rcj2 +
                2 * raj2 * rbj2 + 2 * raj2 * rcj2 + 2 * rbj2 * rcj2) /
               35;

          auto tmp_mat = Ri.Transpose() * Rj *
                         Matrix<ld>(std::vector<ld>{raj2 / 5, 0, 0, 0, rbj2 / 5,
                                                    0, 0, 0, rcj2 / 5},
                                    3, 3) *
                         Rj.Transpose() * Ri *
                         Matrix<ld>(std::vector<ld>{rai2 / 5, 0, 0, 0, rbi2 / 5,
                                                    0, 0, 0, rci2 / 5},
                                    3, 3);
          E -= bodies[i].mass * bodies[j].mass * 3.0 / 4.0 / div5 * tmp_mat.Trace();

          tmp_mat = (Ri *
                     Matrix<ld>(
                         std::vector<ld>{rai2 * (3 * rai2 + rbi2 + rci2), 0, 0,
                                         0, rbi2 * (rai2 + 3 * rbi2 + rci2), 0,
                                         0, 0, rci2 * (rai2 + rbi2 + 3 * rci2)},
                         3, 3) *
                     Ri.Transpose()) /
                    35;
          E += bodies[i].mass * bodies[j].mass * 15.0 / 4.0 / div7 *
               (qi_qj.DotProduct(tmp_mat * (qi_qj)));

          tmp_mat = (Rj *
                     Matrix<ld>(
                         std::vector<ld>{raj2 * (3 * raj2 + rbj2 + rcj2), 0, 0,
                                         0, rbj2 * (raj2 + 3 * rbj2 + rcj2), 0,
                                         0, 0, rcj2 * (raj2 + rbj2 + 3 * rcj2)},
                         3, 3) *
                     Rj.Transpose()) /
                    35;
          E += bodies[i].mass * bodies[j].mass * 15.0 / 4.0 / div7 *
               (qi_qj.DotProduct(tmp_mat * (qi_qj)));

          tmp_mat = (rai2 + rbi2 + rci2) / 5 * Rj *
                        Matrix<ld>(
                            std::vector<ld>{raj2, 0, 0, 0, rbj2, 0, 0, 0, rcj2},
                            3, 3) /
                        5 * Rj.Transpose() +
                    (raj2 + rbj2 + rcj2) / 5 * Ri *
                        Matrix<ld>(
                            std::vector<ld>{rai2, 0, 0, 0, rbi2, 0, 0, 0, rci2},
                            3, 3) /
                        5 * Ri.Transpose();
          E += bodies[i].mass * bodies[j].mass * 15.0 / div7 *
               (qi_qj.DotProduct(tmp_mat * (qi_qj)));

          tmp_mat = qi_qj2 * Ri *
                    Matrix<ld>(std::vector<ld>{rai2 * rai2 * 3 / 35, 0, 0, 0,
                                               rbi2 * rbi2 * 3 / 35, 0, 0, 0,
                                               rci2 * rci2 * 3 / 35},
                               3, 3) *
                    Ri.Transpose() * qi_qj2;
          E -= bodies[i].mass * bodies[j].mass * 35.0 / 8.0 / div9 * tmp_mat.Trace();

          tmp_mat = qi_qj2 * Rj *
                    Matrix<ld>(std::vector<ld>{raj2 * raj2 * 3 / 35, 0, 0, 0,
                                               rbj2 * rbj2 * 3 / 35, 0, 0, 0,
                                               rcj2 * rcj2 * 3 / 35},
                               3, 3) *
                    Rj.Transpose() * qi_qj2;
          E -= bodies[i].mass * bodies[j].mass * 35.0 / 8.0 / div9 * tmp_mat.Trace();

          tmp_mat = Ri *
                    Matrix<ld>(std::vector<ld>{rai2 / 5.0, 0, 0, 0, rbi2 / 5.0,
                                               0, 0, 0, rci2 / 5.0},
                               3, 3) *
                    Ri.Transpose() * qi_qj2 * Rj *
                    Matrix<ld>(std::vector<ld>{raj2 / 5.0, 0, 0, 0, rbj2 / 5.0,
                                               0, 0, 0, rcj2 / 5.0},
                               3, 3) *
                    Rj.Transpose();
          E -= bodies[i].mass * bodies[j].mass * 105.0 / 4.0 / div9 *
               qi_qj.DotProduct(tmp_mat * qi_qj);
        }
      }
    }
  }
  return kG * E;
}

Vec3<ld> System::LinearMomentum() const {
  Vec3<ld> ret(0, 0, 0);
  for (auto i = 0u; i < bodies.size(); i++) {
    ret += bodies[i].p;
  }
  return ret;
}

Vec3<ld> System::AngularMomentum() const {
  Vec3<ld> ret(0, 0, 0);
  for (auto i = 0u; i < bodies.size(); i++) {
    ret += bodies[i].q.CrossProduct(bodies[i].p);
    if (bodies[i].rigid)
      ret += bodies[i].R * bodies[i].pi;
  }
  return ret;
}

ld System::GetAxialTilt(int idx) const {
  return acos(bodies[idx].R.GetElement(3, 3));
}

Matrix<ld> System::GetEulerAngles(int idx) const {
  ld phi, theta, psi;
  RotMat2EulerAngles(bodies[idx].R, phi, theta, psi);
  return Matrix<ld>(std::vector<ld>{phi, theta, psi}, 3, 1);
}

Vec3<ld> System::GetPos(int idx) const { return bodies[idx].q; }

Vec3<ld> System::GetVel(int idx) const {
  return bodies[idx].p / bodies[idx].mass;
}

Vec3<ld> System::GetFtf(int idx) const {
  return TidalForceEquilibrium(idx, 0);
}

ld System::GetEheat(int idx) const {
    Vec3<ld> spin1(0, 0, 0);
    spin1 += bodies[idx].R * bodies[idx].Iinv * bodies[idx].pi;
    Vec3<ld> vel1(0, 0, 0);
    vel1 += bodies[idx].p / bodies[idx].mass;
    Vec3<ld> pos1(0, 0, 0);
    pos1 += bodies[idx].q;
    Vec3<ld> ddot(0, 0, 0);
    ddot=vel1-(spin1.CrossProduct(pos1));
    Vec3<ld> F1(0, 0, 0);
    F1 += TidalForceEquilibrium(idx, 0);
    ld ret=ddot.DotProduct(F1);
    
  return ret;
}

std::ostream &System::operator>>(std::ostream &os) const {
  for (const auto &s : bodies)
    os << s;
  return os;
}

std::istream &System::operator<<(std::istream &is) {
  for (auto &s : bodies)
    is >> s;
  return is;
}

void System::RK4(int steps) {
  // K: delta_q, delta_p, delta_pi, delta_R

  ld h = GetStepSize();
  std::vector<Vec3<ld>> kq1, kq2, kq3, kq4;
  std::vector<Vec3<ld>> kp1, kp2, kp3, kp4;
  std::vector<Vec3<ld>> kpi1, kpi2, kpi3, kpi4;
  std::vector<Matrix<ld>> kR1, kR2, kR3, kR4;
  kq1.resize(bodies.size());
  kp1.resize(bodies.size());
  kpi1.resize(bodies.size());
  kR1.resize(bodies.size());
  kq2.resize(bodies.size());
  kp2.resize(bodies.size());
  kpi2.resize(bodies.size());
  kR2.resize(bodies.size());
  kq3.resize(bodies.size());
  kp3.resize(bodies.size());
  kpi3.resize(bodies.size());
  kR3.resize(bodies.size());
  kq4.resize(bodies.size());
  kp4.resize(bodies.size());
  kpi4.resize(bodies.size());
  kR4.resize(bodies.size());

  for (int step_now = 0; step_now < steps; step_now++) {

    StepwiseOperations();
    current_time += h;

    EoM(kq1, kp1, kpi1, kR1);
    for (auto i = 0u; i < bodies.size(); i++) {
      bodies[i].q += kq1[i] / 2.0 * h;
      bodies[i].p += kp1[i] / 2.0 * h;
      if (bodies[i].rigid) {
        bodies[i].pi += kpi1[i] / 2.0 * h;
        bodies[i].R += kR1[i] / 2.0 * h;
      }
    }

    EoM(kq2, kp2, kpi2, kR2);
    for (auto i = 0u; i < bodies.size(); i++) {
      bodies[i].q += (kq2[i] / 2.0 - kq1[i] / 2.0) * h;
      bodies[i].p += (kp2[i] / 2.0 - kp1[i] / 2.0) * h;
      if (bodies[i].rigid) {
        bodies[i].pi += (kpi2[i] / 2.0 - kpi1[i] / 2.0) * h;
        bodies[i].R += (kR2[i] / 2.0 - kR1[i] / 2.0) * h;
      }
    }

    EoM(kq3, kp3, kpi3, kR3);
    for (auto i = 0u; i < bodies.size(); i++) {
      bodies[i].q += (kq3[i] - kq2[i] / 2.0) * h;
      bodies[i].p += (kp3[i] - kp2[i] / 2.0) * h;
      if (bodies[i].rigid) {
        bodies[i].pi += (kpi3[i] - kpi2[i] / 2.0) * h;
        bodies[i].R += (kR3[i] - kR2[i] / 2.0) * h;
      }
    }

    EoM(kq4, kp4, kpi4, kR4);
    for (auto i = 0u; i < bodies.size(); i++) {
      bodies[i].q -= kq3[i] * h;
      bodies[i].p -= kp3[i] * h;
      if (bodies[i].rigid) {
        bodies[i].pi -= kpi3[i] * h;
        bodies[i].R -= kR3[i] * h;
      }
    }

    for (auto i = 0u; i < bodies.size(); i++) {
      bodies[i].q +=
          (kq1[i] / 6.0 + kq2[i] / 3.0 + kq3[i] / 3.0 + kq4[i] / 6.0) * h;
      bodies[i].p +=
          (kp1[i] / 6.0 + kp2[i] / 3.0 + kp3[i] / 3.0 + kp4[i] / 6.0) * h;
      if (bodies[i].rigid) {
        bodies[i].pi +=
            (kpi1[i] / 6.0 + kpi2[i] / 3.0 + kpi3[i] / 3.0 + kpi4[i] / 6.0) * h;
        bodies[i].R +=
            (kR1[i] / 6.0 + kR2[i] / 3.0 + kR3[i] / 3.0 + kR4[i] / 6.0) * h;
      }
    }
  }
}

void System::EoM(std::vector<Vec3<ld>> &kq, std::vector<Vec3<ld>> &kp,
                 std::vector<Vec3<ld>> &kpi, std::vector<Matrix<ld>> &kR) {

  for (auto i = 0u; i < bodies.size(); i++) {
    kq[i] = bodies[i].p / bodies[i].mass;
    kp[i] = PVPq(i) * (ld(-1));
    if (bodies[i].rigid) {
      kpi[i] = bodies[i].pi.CrossProduct(bodies[i].Iinv * bodies[i].pi) -
               Rot(bodies[i].R.Transpose() * PVPR(i));
      kR[i] = bodies[i].R * (bodies[i].Iinv * bodies[i].pi).Skew();
    }
  }
}

void System::EnforceTide(ld h) {
  // for (auto i = 0u; i < bodies.size(); i++)
  //     if (bodies[i].rigid) for (auto j = 0u; j < bodies.size(); j++)
  //         if  (i!=j) {
  //             auto tmp = tidalForceEquilibrium(i, j);
  //             auto mu = (bodies[i].mass *
  //             bodies[j].mass)/(bodies[i].mass + bodies[j].mass); auto
  //             F_tides = tmp; // * bodies[i].mass; auto N_tides =
  //             bodies[i].R.transpose() *
  //             ((bodies[j].q-bodies[i].q).crossProduct(tmp));
  //             bodies[i].p = bodies[i].p - F_tides * h * mu; bodies[j].p
  //             = bodies[j].p + F_tides * h * mu; bodies[i].pi =
  //             bodies[i].pi - N_tides * mu * h;
  //         }

  for (auto i = 0u; i < bodies.size(); i++)
    for (auto j = 0u; j < bodies.size(); j++)
      if (i != j && bodies[i].rigid &&
          tidal_pairs.find(std::make_pair(bodies[i].name, bodies[j].name)) !=
              tidal_pairs.end()) {
        auto tmp = TidalForceEquilibrium(i, j);
        auto mu = (bodies[i].mass * bodies[j].mass) /
                  (bodies[i].mass + bodies[j].mass);
        auto F_tides = tmp; // * bodies[i].mass;
        auto N_tides = bodies[i].R.Transpose() *
                       ((bodies[j].q - bodies[i].q).CrossProduct(tmp));
        bodies[i].p = bodies[i].p - F_tides * h * mu;
        bodies[j].p = bodies[j].p + F_tides * h * mu;
        bodies[i].pi = bodies[i].pi - N_tides * mu * h;
      }
}

Vec3<ld> System::TidalForceEquilibrium(int id_host, int id_guest) const {

  const auto m_host = bodies[id_host].mass;
  const auto m_guest = bodies[id_guest].mass;
  const auto rad_vec = bodies[id_guest].q - bodies[id_host].q;
  const auto vel_vec =
      bodies[id_guest].p / m_guest - bodies[id_host].p / m_host;
  const auto r = rad_vec.Norm();
  const auto r2 = r * r, r5 = r2 * r2 * r, r10 = r5 * r5;
  const auto radius_host = bodies[id_host].radius,
             rh2 = radius_host * radius_host, rh5 = rh2 * rh2 * radius_host;
  ld Q = bodies[id_host].Q;
  const auto sigma =
      bodies[id_host].time_lag * 4.0 * kG / 3.0 / rh5 * (1 - Q) / Q;
  // `Omega_host`: the angular velocity of host
  const auto Omega_host =
      bodies[id_host].R * bodies[id_host].Iinv * bodies[id_host].pi;
  const auto A = rh5 * Q / (1 - Q);
  const auto mu = m_host * m_guest / (m_host + m_guest);
  auto h = rad_vec.CrossProduct(vel_vec);

  auto F_tides = (rad_vec * (rad_vec.DotProduct(vel_vec) * 3.0) +
                  (h - Omega_host * r2).CrossProduct(rad_vec)) *
                 (-9.0L * sigma * m_guest * m_guest * A * A / 2.0 / mu / r10);

  return F_tides;
//    std::cout << "hahaha you are including tides\n";
}

json System::OutputDataJson() const {
  json dynamics;
  dynamics["time"] = current_time;
  if (output_format.Hamiltonian)
    dynamics["Hamiltonian"] = Hamiltonian();
  if (output_format.momentum) {
    dynamics["linear_momentum"] = LinearMomentum();
    dynamics["angular_momentum"] = AngularMomentum();
  }
  if (output_format.customize) {
    dynamics["customize"] = CustomizeData(*this);
  }

  std::vector<std::vector<ld>> orbs;
  // Calculate orbital elements.
  // Central coordinates.
  if (parameters.coordinates == CENTRAL) {
    // The first body is set to be the central body.
    unsigned int center_id = 0;
    ld center_mass = GetMass(center_id);

    for (unsigned int idx = 0; idx < bodies.size(); idx++)
      if (idx != center_id) {
        json j;
        std::vector<ld> orb;
        orb = PosVel2OrbElements(GetPos(idx) - GetPos(center_id),
                                 GetVel(idx) - GetVel(center_id), center_mass,
                                 GetMass(idx));
        orbs.push_back(orb);
      } else
        orbs.push_back({0, 0, 0, 0, 0, 0});
  } else if (parameters.coordinates == BARYCENTRIC) {
    // center of mass
    Vec3<ld> CoM = {0, 0, 0};
    Vec3<ld> avg_vel = {0, 0, 0};
    ld mass_sum = 0;
    for (unsigned int idx = 0; idx < bodies.size(); idx++) {
      CoM = CoM + GetPos(idx) * GetMass(idx);
      avg_vel = avg_vel + GetVel(idx) * GetMass(idx);
      mass_sum = mass_sum + GetMass(idx);
    }
    CoM = CoM / mass_sum;
    avg_vel = avg_vel / mass_sum;
    for (unsigned int idx = 0; idx < bodies.size(); idx++) {
      auto orb = PosVel2OrbElements(GetPos(idx) - CoM, GetVel(idx) - avg_vel,
                                    GetMass(0), GetMass(idx));
      orbs.push_back(orb);
    }
  } else if (parameters.coordinates == JACOBI) {
    // calc orbital elements using Jacobi coordinates
    // first body
    orbs.push_back({0, 0, 0, 0, 0, 0});

    unsigned int idx = 0;
    ld mass = GetMass(idx);
    auto r = GetPos(idx) * GetMass(idx);
    auto mom = GetVel(idx) * GetMass(idx);

    // other bodies
    for (idx = 1; idx < bodies.size(); idx++) {
      json j;

      std::vector<ld> orb;
      orb = PosVel2OrbElements((r * (1 / mass) - GetPos(idx)) * ld(-1.0),
                               (mom * (1 / mass) - GetVel(idx)) * ld(-1.0),
                               mass, GetMass(idx));
      orbs.push_back(orb);

      r = r + GetPos(idx) * GetMass(idx);
      mom = mom + GetVel(idx) * GetMass(idx);
      mass = mass + GetMass(idx);
    }
  } else {
    std::cerr << "Please input a valid coordinate system.";
    exit(1);
  }

  for (unsigned int idx = 0; idx < bodies.size(); idx++) {
    json j;
    auto i = orbs[idx][2], Omega = orbs[idx][3];
    if (output_format.orbital_elements)
      j["orbital_elements"] = orbs[idx];
    if (output_format.position)
      j["pos"] = GetPos(idx);
    if (output_format.velocity)
      j["vel"] = GetVel(idx);
    if (bodies[idx].rigid) {
      Vec3<ld> v1 = bodies[idx].R * kE3;
      Vec3<ld> v2(sin(i) * sin(Omega), -sin(i) * cos(Omega), cos(i));
      Vec3<ld> v3 = (bodies[idx].R * bodies[idx].Iinv * bodies[idx].pi);
      Vec3<ld> v3N = v3;
      v3N.Normalize();
      Vec3<ld> v4 = GetPos(idx);
      v4.Normalize();

      if (output_format.obliquity)
        j["obliquity"] = acos(v2.DotProduct(v3N));
      if (output_format.axial_tilt)
        j["axial_tilt"] = acos(v1.DotProduct(v4));
      if (output_format.Ftf)
        j["Ftf"] = GetFtf(idx);
      if (output_format.Eheat)
        j["Eheat"] = GetEheat(idx);
      if (output_format.spin)
        j["spin"] = bodies[idx].R * bodies[idx].Iinv * bodies[idx].pi;
      if (output_format.axis)
        j["axis"] = bodies[idx].R * kE3;
    }
    dynamics["body"][GetName(idx)] = j;
  }

  return dynamics;
}

void System::OutputInitSysPosVel(const std::string &dir) const {
  json j = Dump();
  for (unsigned int idx = 0; idx < bodies.size(); idx++) {
    auto &tmp = j["body"][idx];
    tmp.erase("_R");
    tmp.erase("_pi");
    tmp.erase("_p");
  }
  auto file = std::ofstream(dir + "init_system_posvel.json");
  file << j.dump(4);
}

void System::InitMatFile() const {

  // Create the file to store pos,vel,etc of each body.
  std::vector<std::ofstream> body_files;
  for (unsigned int idx = 0; idx < GetNum(); idx++) {
    body_files.push_back(std::ofstream(sys_dir + "/data_in_mat/" + GetName(idx),
                                       std::ios::trunc));
    // Print the header line.
    body_files[idx] << std::left << std::setw(30) << "time";
    if (output_format.orbital_elements)
      body_files[idx] << std::setw(30) << "a" << std::setw(30) << "e"
                      << std::setw(30) << "i" << std::setw(30) << "Omega"
                      << std::setw(30) << "omega" << std::setw(30)
                      << "true_anomaly";
    if (output_format.position)
      body_files[idx] << std::setw(30) << "x" << std::setw(30) << "y"
                      << std::setw(30) << "z";
    if (output_format.velocity)
      body_files[idx] << std::setw(30) << "u" << std::setw(30) << "v"
                      << std::setw(30) << "w";
    if (bodies[idx].rigid) {
      if (output_format.axis)
        body_files[idx] << std::setw(30) << "axis_x" << std::setw(30)
                        << "axis_y" << std::setw(30) << "axis_z";
      if (output_format.spin)
        body_files[idx] << std::setw(30) << "spin_x" << std::setw(30)
                        << "spin_y" << std::setw(30) << "spin_z";
      if (output_format.obliquity)
        body_files[idx] << std::setw(30) << "obliquity";
      if (output_format.axial_tilt)
        body_files[idx] << std::setw(30) << "axial_tilt";
      if (output_format.Ftf)
          body_files[idx] << std::setw(30)  << "Fx" << std::setw(30) << "Fy"
                        << std::setw(30) << "Fz";
      if (output_format.Eheat)
            body_files[idx] << std::setw(30) << "Eheat";
    }
    body_files[idx] << std::endl;
  }

  // Create the file to store the Hamiltonian.
  std::ofstream Hamiltonian_file;
  if (output_format.Hamiltonian) {
    Hamiltonian_file =
        std::ofstream(sys_dir + "/data_in_mat/Hamiltonian", std::ios::trunc);
    // Print the header line.
    Hamiltonian_file << std::setw(30) << "time" << std::setw(30) << "Hamiltonian" << std::endl;
  }

  // Create the file to store the momentum.
  std::ofstream momentum_file;
  if (output_format.momentum) {
    momentum_file =
        std::ofstream(sys_dir + "/data_in_mat/momentum", std::ios::trunc);
    // Print the header line.
    momentum_file << std::setw(30) << "time" << std::setw(30)
                  << "linear_momentum_x" << std::setw(30) << "linear_momentum_y"
                  << std::setw(30) << "linear_momentum_z" << std::setw(30)
                  << "angular_momentum_x" << std::setw(30)
                  << "angular_momentum_y" << std::setw(30)
                  << "angular_momentum_z" << std::endl;
  }

  // Create the file to store the customized output data.
  std::ofstream customize_file;
  if (output_format.customize) {
    customize_file =
        std::ofstream(sys_dir + "/data_in_mat/customize", std::ios::trunc);
    customize_file << std::setw(30) << "time" << std::setw(30)
                   << "your_customize_data_array" << std::endl;
  }

  // Output information to `body_files` line by line from "data_tmp.json"
  // if "data_tmp.json" exists.
  std::ifstream fin(sys_dir + "/data_tmp.json");
  if (fin) {
    std::string str;
    while (std::getline(fin, str, '\n')) {
      auto data = json::parse(str);
      OutputBodyData(data, body_files);
      if (output_format.Hamiltonian)
        OutputHamiltonian(data, Hamiltonian_file);
      if (output_format.momentum)
        OutputMomentum(data, momentum_file);
      if (output_format.customize)
        OutputCustomizeData(data, customize_file);
    }
    fin.close();
  }

  for (auto &file : body_files) {
    file.close();
  }
  if (output_format.Hamiltonian)
    Hamiltonian_file.close();
  if (output_format.momentum)
    momentum_file.close();
  if (output_format.customize)
    customize_file.close();
}

void System::AppendMatFile(const std::vector<json> &j) const {

  // Open the file to append pos,vel,etc of each body.
  std::vector<std::ofstream> body_files;
  for (unsigned int idx = 0; idx < GetNum(); idx++) {
    body_files.push_back(
        std::ofstream(sys_dir + "/data_in_mat/" + GetName(idx), std::ios::app));
  }
  for (auto &e : j)
    OutputBodyData(e, body_files);
  for (auto &file : body_files)
    file.close();

  // Append the Hamiltonian.
  if (output_format.Hamiltonian) {
    std::ofstream Hamiltonian_file;
    Hamiltonian_file =
        std::ofstream(sys_dir + "/data_in_mat/Hamiltonian", std::ios::app);
    for (auto &e : j)
      OutputHamiltonian(e, Hamiltonian_file);
    Hamiltonian_file.close();
  }

  // Append the momentum.
  if (output_format.momentum) {
    std::ofstream momentum_file;
    momentum_file =
        std::ofstream(sys_dir + "/data_in_mat/momentum", std::ios::app);
    for (auto &e : j)
      OutputMomentum(e, momentum_file);
    momentum_file.close();
  }

  // Append the customized output data.
  if (output_format.customize) {
    std::ofstream customize_file;
    customize_file =
        std::ofstream(sys_dir + "/data_in_mat/customize", std::ios::app);
    for (auto &e : j)
      OutputCustomizeData(e, customize_file);
    customize_file.close();
  }
}

void System::OutputBodyData(const json &j,
                            std::vector<std::ofstream> &fouts) const {
  for (unsigned int idx = 0; idx < fouts.size(); idx++) {
    fouts[idx] << std::setprecision(std::numeric_limits<double>::digits10 + 1)
               << std::scientific << std::left << std::setw(30) << double(j["time"]);
    auto body_json = j["body"][GetName(idx)];
    if (output_format.orbital_elements) {
      auto orb = body_json["orbital_elements"];
      for (auto e : orb) {
        fouts[idx] << std::setw(30) << double(e);
      }
    }
    if (output_format.position) {
      auto pos = body_json["pos"];
      for (json::iterator it = pos.begin(); it != pos.end(); it++) {
        fouts[idx] << std::setw(30) << double(it.value());
      }
    }
    if (output_format.velocity) {
      auto vel = body_json["vel"];
      for (json::iterator it = vel.begin(); it != vel.end(); it++) {
        fouts[idx] << std::setw(30) << double(it.value());
      }
    }
    if (bodies[idx].rigid) {
      if (output_format.axis) {
        auto axis = body_json["axis"];
        for (json::iterator it = axis.begin(); it != axis.end(); it++) {
          fouts[idx] << std::setw(30) << double(it.value());
        }
      }
      if (output_format.spin) {
        auto spin = body_json["spin"];
        for (json::iterator it = spin.begin(); it != spin.end(); it++) {
          fouts[idx] << std::setw(30) << double(it.value());
        }
      }
      if (output_format.obliquity)
        fouts[idx] << std::setw(30) << double(body_json["obliquity"]);
      if (output_format.axial_tilt)
        fouts[idx] << std::setw(30) << double(body_json["axial_tilt"]);
      if (output_format.Ftf) {
          auto NFtf = body_json["Ftf"];
          for (json::iterator it = NFtf.begin(); it != NFtf.end(); it++) {
            fouts[idx] << std::setw(30) << double(it.value());
          }
        }
      if (output_format.Eheat)
        fouts[idx] << std::setw(30) << double(body_json["Eheat"]);
    }
    fouts[idx] << std::endl;
  }
}

void System::OutputHamiltonian(const json &j, std::ofstream &fout) const {
  fout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
       << std::left << std::scientific << std::left << std::setw(30)
       << double(j["time"]) << std::setw(30) << j["Hamiltonian"] << std::endl;
}

void System::OutputMomentum(const json &j, std::ofstream &fout) const {
  fout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
       << std::left << std::scientific << std::left << std::setw(30)
       << double(j["time"]);
  auto linear_mom = j["linear_momentum"];
  auto ang_mom = j["angular_momentum"];
  for (int j = 0; j < 3; j++)
    fout << std::setw(30) << linear_mom[j].dump();
  for (int j = 0; j < 3; j++)
    fout << std::setw(30) << ang_mom[j].dump();
  fout << std::endl;
}

void System::OutputCustomizeData(const json &j, std::ofstream &fout) const {
  fout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
       << std::left << std::scientific << std::left << std::setw(30)
       << double(j["time"]);

  auto customize = j["customize"];

  for (json::iterator it = customize.begin(); it != customize.end(); it++) {
    fout << std::setw(30) << it.value();
  }
  fout << std::endl;
}

} // namespace rb_sim
