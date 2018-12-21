#pragma once

#include <vector>
#include <istream>
#include <sstream>

class Rocket
{
public:
  /**
   * @brief Creates a rocket with no stages
   * @details
   */
  Rocket() = default;

  /**
   * @brief Creates a rocket from a formatted string
   * @details Each line is a stage, with either 1 (mass) or 4 values (dry mass, wet mass, thrust and isp).
   * Stages lowermost to topmost, so payload first and descending stages afterwards (Example : payload -> S-IVB -> S-II -> S-IC) 
   * 
   * @param input input string
   */
  Rocket(std::istream &input);

  /**
   * @brief Rocket stage parameters
   * @details [long description]
   * 
   */
  class Stage {
  public:
    Stage() = default;
    Stage(double mass);
    Stage(double dryMass, double wetMass, double thrust, double isp, double ispSL);
    double getDryMass() const { return dryMass; }
    double getWetMass() const { return wetMass; }
    double getThrust() const { return thrust; }
    double getISP() const { return isp; }
    double getISPSeaLevel() const { return ispSL; }
    double getPropellantMass() const { return wetMass-dryMass; }

  private:
    double dryMass;
    double wetMass;
    double thrust;
    double isp;
    double ispSL;
  };
  
  /**
   * @brief Adds a stage under the rocket (will burn before the ones that were already added)
   * @details [long description]
   * 
   * @param stage [description]
   */
  void addStage(const Stage &stage) {
    stages.push_back(stage);
  }

  int firstStage() const {
    return stages.size()-1;
  }
  const Stage& getStage(int stageNumber) const {
    return stages[stageNumber];
  }


private:

  std::vector<Stage> stages;
};

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
    double dryMass, wetMass, thrust, isp, ispSL;
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
      if (!(iss >> ispSL)) error(lineNumber, "Fifth element is invalid");
      stage = Stage(dryMass, wetMass, thrust, isp, ispSL);
    } else {
      stage = Stage(dryMass);
    }
    addStage(stage);
    lineNumber += 1;
  }
}

Rocket::Stage::Stage(double mass): Stage(mass, mass, 0, 0, 0) {}

Rocket::Stage::Stage(double dryMass, double wetMass, double thrust, double isp, double ispSL) {
  if (dryMass < 0 || wetMass < 0 || thrust < 0 || isp < 0 || ispSL < 0 || wetMass < dryMass) throw runtime_error("Incorrect stage parameters");
  this->dryMass = dryMass;
  this->wetMass = wetMass;
  this->thrust = thrust;
  this->isp = isp;
  this->ispSL = ispSL;
}