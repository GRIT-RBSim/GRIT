/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

template <typename T>
std::string ToStringWithPrecision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}

inline std::string Basename(const std::string &pathname) {
  return {std::find_if(pathname.rbegin(), pathname.rend(),
                       [](char c) { return c == '/'; })
              .base(),
          pathname.end()};
}

inline std::string GetFileContents(const char *filename) {
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in.bad())
    return "";

  std::ostringstream contents;
  contents << in.rdbuf();
  return (contents.str());
}

inline bool IsReadable(const std::string &file) {
  return std::ifstream(file.c_str()).good();
}

inline void ProgressBar(std::ostream &os, float percentage, int elapsed_in_s,
                        float current_time) {
  static unsigned int cc = 0;
  const char pp[4] = {'-', '\\', '|', '/'};
  os.flush();
  os << "\r [";
  int i = 0;
  for (; i < percentage * 50; i++)
    os << "=";
  os << ">";
  for (i++; i <= 50; i++)
    os << "-";
  os << "] " << std::setfill(' ') << std::setw(3) << int(percentage * 100)
     << "\%";
  int days = elapsed_in_s / 86400;
  elapsed_in_s %= 86400;
  int hours = elapsed_in_s / 3600;
  elapsed_in_s %= 3600;
  int minutes = elapsed_in_s / 60;
  elapsed_in_s %= 60;
  os << " " << pp[(cc++) % 4] << " " << std::right << std::setfill(' ')
     << std::setw(5) << days;
  os << ":" << std::setfill('0') << std::setw(2) << hours << ":" << std::setw(2)
     << minutes << ":" << std::setw(2) << elapsed_in_s
     << " elapsed. t=" << std::setfill(' ') << std::left << std::setw(10)
     << current_time;
}

#endif
