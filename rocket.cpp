#include "rocket.hpp"
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

static void error(int line, string msg) {
  stringstream msg_format;
  msg_format << msg << " at line " << line;
  throw runtime_error(msg_format.str());
}

Rocket::Rocket(std::istream &input) {
  string line;
  int lineNumber = 1;
  while (getline(input, line)) {
    istringstream iss(line);
    // Ignore comments lines starting with #
    skipws(iss);
    if (iss.peek() == '#') {
      lineNumber += 1;
      continue;
    }
    // Extract parameters
    double dryMass, wetMass, thrust, isp;
    Stage stage;
    iss >> dryMass;
    if (!iss) error(lineNumber, "First element is invalid");
    // Ignore whitespace in case of comment tag
    while (iss.peek() == ' ') iss.ignore(1);
    if (iss && iss.peek() != '#' && iss.peek() != '\r') {
      // 4 value case
      if (!(iss >> wetMass)) error(lineNumber, "Second element is invalid");
      if (!(iss >> thrust)) error(lineNumber, "Third element is invalid");
      if (!(iss >> isp)) error(lineNumber, "Fourth element is invalid");
      stage = Stage(dryMass, wetMass, thrust, isp);
    } else {
      stage = Stage(dryMass);
    }
    addStage(stage);
    lineNumber += 1;
  }
}

void Rocket::addStage(const Stage &stage) {
  stages.push_back(stage);
}

int Rocket::firstStage() const {
  return stages.size()-1;
}

const Rocket::Stage& Rocket::getStage(int id) const {
  return this->stages[id];
}

double Rocket::Stage::getDryMass() const {
  return dryMass;
}

double Rocket::Stage::getWetMass() const {
  return wetMass;
}

double Rocket::Stage::getThrust() const {
  return thrust;
}

double Rocket::Stage::getISP() const {
  return isp;
}

double Rocket::Stage::getPropellantMass() const {
  return wetMass-dryMass;
}

Rocket::Stage::Stage(double mass): Stage(mass, mass, 0, 0) {}

Rocket::Stage::Stage(double dryMass, double wetMass, double thrust, double isp) {
  if (dryMass < 0 || wetMass < 0 || thrust < 0 || isp < 0 || wetMass < dryMass) throw runtime_error("Incorrect stage parameters");
  this->dryMass = dryMass;
  this->wetMass = wetMass;
  this->thrust = thrust;
  this->isp = isp;
}