#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#include "rocket.hpp"

using namespace std;

double rad(double a) {
  return a*M_PI/180;
}

double deg(double a) {
  return a*180/M_PI;
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

  const double target_alt = 185000;

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

  double pitch = 90;
  double T = -10;
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double lastT;
  bool PEGinit = false;

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

  void guidance(double oldT, double r_t) {
    const Rocket::Stage& currentStage = rocket.getStage(currentStageId);

    double r = length(posx, posy);
    double ve = g0*currentStage.getISP();
    double tau = ve/acc;

    double b0 = -ve*log(1-oldT/tau);
    double b1 = b0*tau - ve*oldT;
    double c0 = b0*oldT - b1;
    double c1 = c0*tau - ve*oldT*oldT/2;

    double mbx = -vvel;
    double mby = r_t-r-vvel*oldT;
    float det = (b0*c1-c0*b1);
    A = (c1*mbx-b1*mby)/det;
    B = (b0*mby-c0*mbx)/det;
  }

  void estimation(double r_t, double deltaT, double v_theta_T) {
    const Rocket::Stage& currentStage = rocket.getStage(currentStageId);

    double oldT = T;

    double r = length(posx, posy);
    double ve = g0*currentStage.getISP();
    double tau = ve/acc;

    if (oldT < -5) {
      oldT = 0.9*tau;
      guidance(oldT, r_t);
    }

    // navigation, basis vectors, and additional target conditions
    double h = r*hvel;
    double ht = r_t*v_theta_T;
    double dh = ht - h;
    double rbar = (r + r_t)/2;

    // Vehicle performance
    C = (gravity/(r*r) - hvel*hvel/r) / acc;
    double fr = A + C;

    double CT = (gravity/(r_t*r_t) - v_theta_T*v_theta_T/r_t) / (acc / (1 - oldT/tau));
    double frT = A + B*oldT + CT;
    double frdot = (frT - fr)/oldT;
    double ftheta = 1 - fr*fr/2;
    double fthetadot = -fr*frdot;
    double fthetadotdot = -frdot*frdot/2;

    double currentT = oldT-deltaT;
    double dv = (dh/rbar + ve*currentT*(fthetadot + fthetadotdot*tau) + fthetadotdot*ve*currentT*currentT/2) / (ftheta + fthetadot*tau + fthetadotdot*tau*tau);

    T = tau * (1 - exp(-dv/ve));

    //cout << time << " A = " << A << " B = " << B << " T = " << oldT << " C = " << C << endl;

    // Update until terminal guidance
    if (T >= 7.5) {
      guidance(T, r_t);
    }
  }

  bool pitchProgram() {
    // clear tower
    if (time < 30) pitch = 90.0;
    // first stage pitch program
    else if (currentStageId == 2) {
      pitch -= 0.53*dt;
    }
    // peg
    else {
      double r_t = seaLevel + target_alt;
      double v_theta_T = sqrt(gravity/r_t);

      if (!PEGinit) {
        estimation(r_t, 0, v_theta_T);
        lastT = T;
        PEGinit = true;
      }

      // Major loop
      double cycleTime = 1.0;
      if ((lastT - T) >= cycleTime) {
        estimation(r_t, cycleTime, v_theta_T);
        lastT = T;
      } else {
        T = T - dt;
      }

      // Minor loop
      double sinPitch = A + (lastT - T)*B + C;
      double targetPitch = deg(asin(sinPitch));
      // Check for Nan
      if (targetPitch == targetPitch) {
        pitch = targetPitch;
        // Termination
        if (T <= 0.0) return true;
      }
    }
    return false;
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
    displayHeader();
    bool shutdown  = false;
    while (!shutdown) {
      for (int i=0;i<displaydt && !shutdown;i++) {
        if (!shutdown) {
          shutdown = pitchProgram();
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
    if (peri < target_alt - 5000) {
      cout << "Failure : ";
    } else {
      cout << "Success! ";
    }
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