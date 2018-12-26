#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>
#include <iomanip>

#include "rocket.hpp"
#include "guidance.hpp"

using namespace std;

/**
 * @brief Converts radians to degrees
 * @details [long description]
 * 
 * @param a [description]
 * @return [description]
 */
double deg(double a) {
  return a*180/M_PI;
}

/**
 * @brief Computes length of vector
 * @details [long description]
 * 
 * @param a x component of vector
 * @param b y component of vector
 * 
 * @return [description]
 */
double length(double a, double b) {
  return sqrt(a*a+b*b);
}

class Sim {

private:
  /// Earth GM constant
  const double mu = 3.9860044188E14;
  /// Earth sealevel radius
  const double seaLevel = 6357000;
  /// G constant (for isp and ve calculations)
  const double g0 = 9.80665;

  /// Rocket model used for simulation
  const Rocket& rocket;
  /// Whether engine apply thrust or not
  bool enginesOn = true;
  /// Currently burning stage
  int currentStageId = 0;
  /// Current mass of propellant left in current stage
  double currentPropellantMass = 0;
  /// Position on x axis (inertial, origin is center of planet, meters)
  double posx = 0.0;
  /// Position on y axis (inertial, origin is center of planet, meters)
  double posy = seaLevel;
  /// Velocity on x axis (inertial, m/s)
  double velx = 0.0;
  /// Velocity on y axis (inertial, m/s)
  double vely = 0.0;

  /// Time since launch (seconds)
  double time = 0.0;

  /// Guidance algorithms
  Guidance guidance;

  /// Pitch angle of rocket (radians, pi/2 is straight up, 0 is parallel to ground)
  double pitch = M_PI/2;

  enum class FailReason {
    NO_FAIL,
    GUIDANCE,
    DEVIATION,
    CRASH,
    AERODYNAMIC,
  };

  double getCurrentMass();
  double step(double dt);
  void displayHeader() {
    cout << "TIME(s) DOWNRANGE(km) ALT(km) HVEL(m/s) VVEL(m/s) PITCH(deg) PERI(km) APO(km) e" << endl;
  }
  auto computeOrbital() {
    double r = length(posx, posy);
    double v = length(velx, vely);
    double E = (v*v)/2 - mu/r;
    double a = -mu/(2*E);

    double h = posx*vely-posy*velx;
    double e = sqrt(1+(2*E*h*h)/(mu*mu));
    double apo  = ((1+e)*a)-seaLevel;
    double peri = ((1-e)*a)-seaLevel;
    return make_tuple(apo, peri, e);
  }

  void displayInfo(double hvel, double vvel) {
    double downrange = atan2(posx, posy)*seaLevel;
    double alt = length(posx, posy) - seaLevel;

    double apo, peri, e;
    tie(apo, peri, e) = computeOrbital();

    cout 
      << (int)round(time) << " "
      << downrange/1000 << " "
      << alt/1000 << " "
      << hvel << " " 
      << vvel << " " 
      << deg(pitch) << " "
      << peri/1000 << " " 
      << apo/1000 << " " 
      << e << endl;
  }

  void integrate(double impulse, double dragAcc, double dt);
  double computeExpandedDeltaV();

public:
  Sim(const Rocket &rocket) : rocket(rocket) {
    currentStageId = rocket.firstStage();
    currentPropellantMass = rocket.getStage(currentStageId).getPropellantMass();
  }
  
  void launch(double start, double end, double pitchRate, double targetAltitude);
};

double earthPressure(double alt) {
  return 101325*exp(-0.0289644*9.80665*alt/(288.15*8.31447));
}

double earthDensity(double alt) {
  return 1.2250*exp(-0.0289644*9.80665*alt/(288.15*8.31447));
}

double Sim::getCurrentMass() {
  const Rocket::Stage& currentStage = rocket.getStage(currentStageId);
  double payloadMass = currentStage.getDryMass();
  for (int i=currentStageId-1;i>=0;i--) {
    payloadMass += rocket.getStage(i).getWetMass();
  }
  return payloadMass + currentPropellantMass;
}

double Sim::step(double dt) {
  // No impulse if engines off, duh
  if (!enginesOn) return 0;

  // Staging (instantly stage empty stages)
  while (currentPropellantMass <= 0.0 && currentStageId > 0) {
    cout << "STAGED #" << currentStageId << " at " << time << " seconds" << endl;
    currentStageId -= 1;
    currentPropellantMass = rocket.getStage(currentStageId).getPropellantMass();
  }

  const Rocket::Stage& currentStage = rocket.getStage(currentStageId);
  // Get mass flow rate
  double ve = g0*currentStage.getISP();
  double massflowrate = currentStage.getThrust()/ve;
  // Accumulate wet mass of upper stages
  double payloadMass = currentStage.getDryMass();
  for (int i=currentStageId-1;i>=0;i--) {
    payloadMass += rocket.getStage(i).getWetMass();
  }
  // Get mass at start and end of impulse
  double mass0 = currentPropellantMass;
  double mass1 = max(currentPropellantMass-massflowrate*dt, 0.0);
  double totalMass0 = payloadMass+mass0;
  double totalMass1 = payloadMass+mass1;
  // Get actual exhaust velocity
  double pressure = earthPressure(length(posx, posy)-seaLevel)/earthPressure(0);
  double atmoVe = g0*(currentStage.getISPSeaLevel()*pressure + currentStage.getISP()*(1-pressure));
  // Impulse (m/s)
  double impulse = atmoVe*log(totalMass0/totalMass1);
  // Remove propellant mass
  currentPropellantMass = mass1;
  return impulse;
}

void Sim::integrate(double impulse, double dragAcc, double dt) {
  double r2 = posx*posx+posy*posy;
  double r = sqrt(r2);
  double dirx = posx/r;
  double diry = posy/r;

  double cosp = cos(pitch);
  double sinp = sin(pitch);

  // Get thrust direction
  double thrustDirx = diry*cosp+dirx*sinp;
  double thrustDiry = diry*sinp-dirx*cosp;
  // Gravity
  double gravity = dt*mu/r2;
  // Subtract gravity and drag
  double accxdt = thrustDirx*impulse - dirx*gravity - velx*dragAcc;
  double accydt = thrustDiry*impulse - diry*gravity - vely*dragAcc;

  velx += accxdt;
  vely += accydt;

  time += dt;
  posx += dt*velx + 0.5*accxdt*dt;
  posy += dt*vely + 0.5*accydt*dt;
}

void Sim::launch(double start, double end, double pitchRate, double targetAltitude) {
  double targetRadius = seaLevel + targetAltitude;
  double targetHoriVel = sqrt(mu/targetRadius);

  guidance.setPitchProgram(start, end, pitchRate);
  guidance.setPEGParameters(targetRadius, targetHoriVel, 1.0);

  const double dt = 0.02;
  const double displaydt = 1.0;

  displayHeader();
  double displayTime = 0.0;
  double acc = 0.0;
  FailReason reason = FailReason::NO_FAIL;
  while (reason == FailReason::NO_FAIL && enginesOn) {
    // Get data for guidance
    double r = length(posx, posy);
    if (r - seaLevel < 0) reason = FailReason::CRASH;
    double ve = g0*rocket.getStage(currentStageId).getISP();
    double hvel = (posy*velx-posx*vely)/r;
    double vvel = (posx*velx+posy*vely)/r;
    // Check AoA
    double aoa = abs(atan2(vvel, hvel)-pitch);
    // Dynamic pressure
    double airDensity = earthDensity(r-seaLevel);
    double dynPressure = 0.5*airDensity*(velx*velx+vely*vely);
    if (sin(aoa)*dynPressure >= 1000) reason = FailReason::AERODYNAMIC;
    // Run guidance
    double dtEngines = guidance.run(time, dt, mu, r, ve, acc, hvel, vvel);
    // Peg fail
    if (guidance.failed()) reason = FailReason::GUIDANCE;
    // Set pitch
    const double maxPitchRate = rad(2.0);
    pitch += min(maxPitchRate*dt, max(-maxPitchRate*dt, guidance.getPitch()-pitch));
    // Shutdown engines if guidance requested 
    if (dtEngines <= 0.0) enginesOn = false;
    // Drag calculation
    double dragAcc = dt*0.5*airDensity*sqrt(velx*velx+vely*vely)*rocket.getStage(currentStageId).getDrag()/getCurrentMass();
    // Engine function generating thrust
    double impulse = step(dtEngines);
    // Integrate rocket position and velocity
    integrate(impulse, dragAcc, dtEngines);
    acc = impulse/dtEngines;
    // Shutdown engines if last stage expanded
    if (currentStageId == 0 && currentPropellantMass <= 0.0) enginesOn = false;
    // Display log
    displayTime += dt;
    if (displayTime >= displaydt) {
      displayInfo(hvel, vvel);
      displayTime = 0;
    }
  }
  if (reason == FailReason::NO_FAIL) {
    double apo, peri, e;
    tie(apo, peri, e) = computeOrbital();
    cout << "SHUTDOWN" << endl;
    cout << "Insertion in a " << peri/1000 << "kmx" << apo/1000 << "km (e=" << e << ") orbit in " << time << "s" << endl;
    cout << "Remaining propellant mass : " << currentPropellantMass << endl;
    cout << "Delta-V expanded at cutoff : " << computeExpandedDeltaV() << "m/s" << endl;
  } else if (reason == FailReason::GUIDANCE) {
    cout << "Flight terminated because guidance can't find a trajectory" << endl;
  } else if (reason == FailReason::CRASH) {
    cout << "Rocket crashed" << endl;
  } else if (reason == FailReason::AERODYNAMIC) {
    cout << "Aerodynamic failure due to aoa too big in thick parts of the atmosphere" << endl;
  }
}

double Sim::computeExpandedDeltaV() {
  const Rocket::Stage& currentStage = rocket.getStage(currentStageId);
  double currentMass = 0.0;
  for (int i=0;i<currentStageId;i++) currentMass += rocket.getStage(i).getWetMass();
  double deltaV = currentStage.getISP()*g0*log((currentStage.getWetMass()+currentMass)/(currentPropellantMass+currentStage.getDryMass()+currentMass));
  currentMass += currentStage.getWetMass();
  for (int i=currentStageId+1;i<=rocket.firstStage();i++) {
    const Rocket::Stage& s = rocket.getStage(i);
    deltaV += s.getISP()*g0*log((s.getWetMass()+currentMass)/(s.getDryMass()+currentMass));
    currentMass += s.getWetMass();
  }
  return deltaV;
}

int main(int argc, char **argv) {
  // Saturn V for skylab configuration
  double start = 30;
  double end = 167;
  double pitchRate = 0.42;
  double targetAltitude = 438000;

  if (argc < 2) {
    cout << "No file given" << endl;
    exit(0);
  } else if (argc == 6) {
    start = atof(argv[2]);
    end = atof(argv[3]);
    pitchRate = atof(argv[4]);
    targetAltitude = atof(argv[5]);
  }

  ifstream file(argv[1]);
  Rocket rocket(file);
  file.close();

  Sim sim(rocket);
  sim.launch(start, end, pitchRate, targetAltitude);
  return 0;
}