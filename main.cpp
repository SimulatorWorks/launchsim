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

  const double gravity = 3.9860044188E14;
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

  double downrange = 0;
  double alt = 0;

  double time = 0.0;
  const double dt = 0.1;
  const int displaydt = 10;

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

    double currentMass = 0;

    // Staging
    while (currentPropellantMass <= 0.0 && currentStageId >= 0) {
      cout << "STAGED #" << currentStageId << endl;
      currentStageId -= 1;
      currentPropellantMass = rocket.getStage(currentStageId).getPropellantMass();
    }

    // If last stage expanded, its over
    if (currentStageId == -1) {
      currentMass = rocket.getStage(0).getDryMass();
      return;
    }

    const Rocket::Stage& currentStage = rocket.getStage(currentStageId);

    // Get current mass
    currentMass = currentStage.getDryMass()+currentPropellantMass;
    // Accumulate wet mass of upper stages
    for (int i=currentStageId-1;i>=0;i--) {
      currentMass += rocket.getStage(i).getWetMass();
    }
    // Deduce acc from current mass
    acc = currentStage.getThrust()/currentMass;
    // Remove propellant mass
    double massflowrate = currentStage.getThrust()/(g0*currentStage.getISP());
    currentPropellantMass -= dt*massflowrate;
  }

  double length(double a, double b) {
    return sqrt(a*a+b*b);
  }

  void displayHeader() {
    cout << "TIME DOWNRANGE(km) ALT(km) HVEL(m/s) VVEL(m/s) PITCH(deg) PERI(km) APO(km) e" << endl;
  }

  void displayInfo() {
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
    accx = thrustDirx*acc - dirx*gravity/r2;
    accy = thrustDiry*acc - diry*gravity/r2;
    velx += dt*accx;
    vely += dt*accy;

    time += dt;
    posx += dt*velx + 0.5*accx*dt*dt;
    posy += dt*vely + 0.5*accy*dt*dt;

    downrange = atan2(posx, posy)*seaLevel;
    alt = length(posx, posy) - seaLevel;

    hvel = diry*velx-dirx*vely;
    vvel = dirx*velx+diry*vely;
  }

  void computeOrbital() {
    double r = length(posx, posy);
    double v = length(velx, vely);
    double E = (v*v)/2 - gravity/r;
    double a = -gravity/(2*E);

    double h = posx*vely-posy*velx;
    e = sqrt(1+(2*E*h*h)/(gravity*gravity));
    apo  = ((1+e)*a)-seaLevel;
    peri = ((1-e)*a)-seaLevel;
  }

  void launch() {
    double targetRadius = seaLevel + 185000;
    double targetHoriVel = sqrt(gravity/targetRadius);

    guidance.setPitchProgram(30, 170, 0.53);
    guidance.setPEGParameters(targetRadius, targetHoriVel, 1.0, 5.0);

    displayHeader();
    bool shutdown  = false;
    while (!shutdown) {
      for (int i=0;i<displaydt && !shutdown;i++) {
        if (!shutdown) {
          double r = length(posx, posy);
          const Rocket::Stage& currentStage = rocket.getStage(currentStageId);
          double ve = g0*currentStage.getISP();
          shutdown = guidance.run(time, dt, gravity, r, ve, acc, hvel, vvel);
          pitch = guidance.getPitch();
          if (shutdown) {
            enginesOn = false;
          }
        }
        computeOrbital();
        step();
        integrate();
      }
      displayInfo();
    }
    cout << "SHUTDOWN" << endl;
    cout << "Insertion in a " << peri/1000 << "kmx" << apo/1000 << "km (e=" << e << ") orbit in " << time << "s" << endl;
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