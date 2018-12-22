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
    Stage(double dryMass, double wetMass, double thrust, double isp, double ispSL, double drag);
    double getDryMass() const { return dryMass; }
    double getWetMass() const { return wetMass; }
    double getThrust() const { return thrust; }
    double getISP() const { return isp; }
    double getISPSeaLevel() const { return ispSL; }
    double getDrag() const { return drag; }
    double getPropellantMass() const { return wetMass-dryMass; }

  private:
    double dryMass;
    double wetMass;
    double thrust;
    double isp;
    double ispSL;
    double drag;
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
    // Remove comment
    string line2 = line.substr(0, line.find_first_of("#"));
    // Remove leading and trailing spaces
    auto start = line2.find_first_not_of(" \t\r");
    if (start == string::npos) {
      lineNumber += 1;
      continue;
    }
    auto endSpace = line2.find_last_not_of(" \t\r");
    string line3 = line2.substr(start, endSpace-start+1);
    istringstream iss(line3);

    // Extract parameters
    double dryMass, wetMass;
    double drag = 0;
    double thrust = 0;
    double isp = 1;
    double ispSL = 1;
    Stage stage;
    if (!(iss >> dryMass) || dryMass < 0) error(lineNumber, "dry mass invalid");
    wetMass = dryMass;
    if (!iss.eof()) {
      // 4 value case
      if (!(iss >> wetMass) || wetMass < dryMass) error(lineNumber, "wet mass invalid");
      if (!(iss >> thrust) || thrust < 0) error(lineNumber, "thrust invalid");
      if (!(iss >> isp) || isp < 0) error(lineNumber, "isp invalid");
      if (!(iss >> ispSL) || ispSL < 0) error(lineNumber, "sea level isp invalid");
      if (!iss.eof()) {
        if (!(iss >> drag)) error(lineNumber, "drag invalid");
      }
    }
    stage = Stage(dryMass, wetMass, thrust, isp, ispSL, drag);
    addStage(stage);
    lineNumber += 1;
  }
}

Rocket::Stage::Stage(double dryMass, double wetMass, double thrust, double isp, double ispSL, double drag) {
  if (dryMass < 0 || wetMass < 0 || thrust < 0 || isp < 0 || ispSL < 0 || wetMass < dryMass || drag < 0) throw runtime_error("Incorrect stage parameters");
  this->dryMass = dryMass;
  this->wetMass = wetMass;
  this->thrust = thrust;
  this->isp = isp;
  this->ispSL = ispSL;
  this->drag = drag;
}