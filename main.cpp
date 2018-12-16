#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#include "rocket.hpp"
#include "guidance.hpp"

using namespace std;

double rad(double a) {
  return a*M_PI/180;
}

class Sim {

private:
  double acc = 0.0;
  bool enginesOn = true;
  const Rocket& rocket;
  int currentStageId = 0;
  float currentPropellantMass = 0;

  const double mu = 3.9860044188E14;
  const double seaLevel = 6357000;
  const double g0 = 9.80665;

  double posx = 0.0;
  double posy = seaLevel;
  double velx = 0.0;
  double vely = 0.0;
  double thrustDirx = 0.0;
  double thrustDiry = 0.0;
  double accx = 0.0;
  double accy = 0.0;
  double hvel = 0.0;
  double vvel = 0.0;

  double time = 0.0;
  const double dt = 0.02;
  const double displaydt = 1.0;

  Guidance guidance;
  double pitch = 90.0;

  float e;
  float peri;
  float apo;

public:
  Sim(const Rocket &rocket) : rocket(rocket) {
    currentStageId = rocket.firstStage();
    currentPropellantMass = rocket.getStage(currentStageId).getPropellantMass();
  }
  void step() {
    acc = 0;
    if (!enginesOn) return;

    // Staging
    while (currentPropellantMass <= 0.0 && currentStageId >= 0) {
      cout << "STAGED #" << currentStageId << " at " << time << " seconds" << endl;
      currentStageId -= 1;
      currentPropellantMass = rocket.getStage(currentStageId).getPropellantMass();
    }

    // If last stage expanded, its over
    if (currentStageId == -1) {
      return;
    }

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
    // between 0 and dt if burn is cut short
    double timeBurn = (mass0-mass1)/(massflowrate);
    double totalMass0 = payloadMass+mass0;
    double totalMass1 = payloadMass+mass1;
    // Impulse (m/s)
    double impulse = ve*log(totalMass0/totalMass1);

    acc = impulse/dt;
    // Remove propellant mass
    currentPropellantMass = mass1;
  }

  double length(double a, double b) {
    return sqrt(a*a+b*b);
  }

  void displayHeader() {
    cout << "TIME DOWNRANGE(km) ALT(km) HVEL(m/s) VVEL(m/s) PITCH(deg) PERI(km) APO(km) e" << endl;
  }

  void displayInfo() {
    double downrange = atan2(posx, posy)*seaLevel;
    double alt = length(posx, posy) - seaLevel;

    cout << (int)round(time) << " " << downrange/1000 << " " << alt/1000 << " " << hvel << " " << vvel << " " << pitch << " " << peri/1000 << " " << apo/1000 << " " << e << endl;
  }

  void integrate() {
    double dirx = posx/length(posx, posy);
    double diry = posy/length(posx, posy);

    double cosp = cos(rad(pitch));
    double sinp = sin(rad(pitch));

    // gravity
    double radius = length(posx, posy);
    double r2 = radius*radius;
    thrustDirx = diry*cosp+dirx*sinp;
    thrustDiry = diry*sinp-dirx*cosp;
    accx = thrustDirx*acc - dirx*mu/r2;
    accy = thrustDiry*acc - diry*mu/r2;
    velx += dt*accx;
    vely += dt*accy;

    time += dt;
    posx += dt*velx + 0.5*accx*dt*dt;
    posy += dt*vely + 0.5*accy*dt*dt;

    hvel = diry*velx-dirx*vely;
    vvel = dirx*velx+diry*vely;
  }

  void computeOrbital() {
    double r = length(posx, posy);
    double v = length(velx, vely);
    double E = (v*v)/2 - mu/r;
    double a = -mu/(2*E);

    double h = posx*vely-posy*velx;
    e = sqrt(1+(2*E*h*h)/(mu*mu));
    apo  = ((1+e)*a)-seaLevel;
    peri = ((1-e)*a)-seaLevel;
  }

  void launch() {
    double targetRadius = seaLevel + 185000;
    double targetHoriVel = sqrt(mu/targetRadius);

    guidance.setPitchProgram(30, 170, 0.53);
    guidance.setPEGParameters(targetRadius, targetHoriVel, 1.0, 5.0);

    displayHeader();
    double displayTime = 0.0;
    bool shutdown  = false;
    while (!shutdown) {
      if (!shutdown) {
        double r = length(posx, posy);
        const Rocket::Stage& currentStage = rocket.getStage(currentStageId);
        double ve = g0*currentStage.getISP();
        shutdown = guidance.run(time, dt, mu, r, ve, acc, hvel, vvel);
        pitch = guidance.getPitch();
        if (shutdown) {
          enginesOn = false;
        }
      }
      computeOrbital();
      step();
      integrate();
      if (displayTime >= displaydt) {
        displayInfo();
        displayTime = 0;
      } else displayTime += dt;
    }
    cout << "SHUTDOWN" << endl;
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