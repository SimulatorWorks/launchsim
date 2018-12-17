#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>

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

  double time = 0.0;

  Guidance guidance;
  double pitch = M_PI/2;

public:
  Sim(const Rocket &rocket) : rocket(rocket) {
    currentStageId = rocket.firstStage();
    currentPropellantMass = rocket.getStage(currentStageId).getPropellantMass();
  }
  double step(double dt) {
    // No impulse if engines off, duh
    if (!enginesOn) return 0;

    // Staging (instantly stage empty stages)
    while (currentPropellantMass <= 0.0 && currentStageId >= 0) {
      cout << "STAGED #" << currentStageId << " at " << time << " seconds" << endl;
      currentStageId -= 1;
      currentPropellantMass = rocket.getStage(currentStageId).getPropellantMass();
    }

    // If last stage expanded, its over
    if (currentStageId == -1) return 0;

    const Rocket::Stage& currentStage = rocket.getStage(currentStageId);
    // Get change in velocity
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
    // Impulse (m/s)
    double impulse = ve*log(totalMass0/totalMass1);
    // Remove propellant mass
    currentPropellantMass = mass1;
    return impulse;
  }

  double length(double a, double b) {
    return sqrt(a*a+b*b);
  }

  void displayHeader() {
    cout << "TIME DOWNRANGE(km) ALT(km) HVEL(m/s) VVEL(m/s) PITCH(deg) PERI(km) APO(km) e" << endl;
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

    cout << (int)round(time) << " " << downrange/1000 << " " << alt/1000 << " " << hvel << " " << vvel << " " << deg(pitch) << " " << peri/1000 << " " << apo/1000 << " " << e << endl;
  }

  void integrate(double impulse, double dt) {
    double dirx = posx/length(posx, posy);
    double diry = posy/length(posx, posy);

    double cosp = cos(pitch);
    double sinp = sin(pitch);

    // gravity
    double radius = length(posx, posy);
    double r2 = radius*radius;
    double thrustDirx = diry*cosp+dirx*sinp;
    double thrustDiry = diry*sinp-dirx*cosp;
    double accxdt = thrustDirx*impulse - dt*dirx*mu/r2;
    double accydt = thrustDiry*impulse - dt*diry*mu/r2;
    velx += accxdt;
    vely += accydt;

    time += dt;
    posx += dt*velx + 0.5*accxdt*dt;
    posy += dt*vely + 0.5*accydt*dt;
  }

  void launch() {
    double targetRadius = seaLevel + 185000;
    double targetHoriVel = sqrt(mu/targetRadius);

    guidance.setPitchProgram(30, 170, 0.53);
    guidance.setPEGParameters(targetRadius, targetHoriVel, 1.0, 5.0);

    const double dt = 0.02;
    const double displaydt = 1.0;

    displayHeader();
    double displayTime = 0.0;
    double acc = 0.0;
    while (enginesOn) {
      // Get data for guidance
      double r = length(posx, posy);
      double ve = g0*rocket.getStage(currentStageId).getISP();
      double hvel = (posy*velx-posx*vely)/r;
      double vvel = (posx*velx+posy*vely)/r;
      // Run guidance
      double dtEngines = guidance.run(time, dt, mu, r, ve, acc, hvel, vvel);
      pitch = guidance.getPitch();
      // Shutdown engines if guidance requested
      if (dtEngines <= 0.0) enginesOn = false;
      double impulse = step(dtEngines);
      integrate(impulse, dtEngines);
      acc = impulse/dtEngines;
      if (displayTime >= displaydt) {
        displayInfo(hvel, vvel);
        displayTime = 0;
      } else displayTime += dt;
    }
    cout << "SHUTDOWN" << endl;
    double apo, peri, e;
    tie(apo, peri, e) = computeOrbital();
    cout << "Insertion in a " << peri/1000 << "kmx" << apo/1000 << "km (e=" << e << ") orbit in " << time << "s" << endl;
    cout << "Remaining propellant mass : " << currentPropellantMass << endl;
  }
};

int main(int argc, char **argv) {
  // Saturn V for skylab configuration
  if (argc < 2) {
    cout << "No file given" << endl;
    exit(0);
  }
  ifstream file(argv[1]);
  Rocket rocket(file);
  file.close();

  Sim sim(rocket);
  sim.launch();
  return 0;
}