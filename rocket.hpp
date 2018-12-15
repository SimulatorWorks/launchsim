#pragma once

#include <vector>
#include <istream>

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
    Stage(double dryMass, double wetMass, double thrust, double isp);
    double getDryMass() const;
    double getWetMass() const;
    double getThrust() const;
    double getISP() const;
    double getPropellantMass() const;

  private:
    double dryMass;
    double wetMass;
    double thrust;
    double isp;
  };
  
  /**
   * @brief Adds a stage under the rocket (will burn before the ones that were already added)
   * @details [long description]
   * 
   * @param stage [description]
   */
  void addStage(const Stage &stage);

  int firstStage() const;
  const Stage& getStage(int stageNumber) const;


private:

  std::vector<Stage> stages;
};