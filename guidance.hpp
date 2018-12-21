#pragma once

/**
 * @brief Converts degrees to radians
 * @details [long description]
 * 
 * @param a [description]
 * @return [description]
 */
double rad(double a) {
  return a*M_PI/180;
}

class Guidance {
public:
  Guidance() = default;
  /**
   * @brief Set the parameters for the first phase of ascent (vertical ascent and pitching)
   * @details [long description]
   * 
   * @param start time at which vertical ascent ends and pitch program starts
   * @param end end of pitch program
   * @param pitchRate degrees/s pitching rate
   */
  void setPitchProgram(double start, double end, double pitchRate) {
    this->start = start;
    this->end = end;
    this->pitchRate = rad(pitchRate);
  }
  /**
   * @brief Set Powered Explicit Guidance parameters after pitch program ends
   * @details [long description]
   * 
   * @param targetRadius radius of target orbit (meters from center of body)
   * @param targetHoriVel horizontal velocity of target orbit (meters/seconds)
   * @param cycleTime time between major loops
   * @param maxPitchRate maximum change of pitch allowed
   */
  void setPEGParameters(double targetRadius, double targetHoriVel, double cycleTime) {
    this->targetRadius = targetRadius;
    this->targetHoriVel = targetHoriVel;
    this->cycleTime = cycleTime;
  }
  /**
   * @brief Runs guidance
   * @details [long description]
   * 
   * @param time time since launch (s)
   * @param dt delta-t for calculations
   * @param mu gravity constant for current body
   * @param r distance from center of body (m)
   * @param ve exhaust velocity (m/s)
   * @param acc current thrust (m/s²)
   * @param hvel horizontal velocity (m/s)
   * @param vvel vertical velocity (m/s)
   * @return time for which engines should run
   */
  double run(double time, double dt, double mu, double r, double ve, double acc, double hvel, double vvel);
  /**
   * @brief Returns target pitch of vehicle
   * @details [long description]
   * @return [description]
   */
  double getPitch() {
    return pitch;
  }

  /**
   * @brief Returns true if PEG failed to compute parameters
   * @details [long description]
   * @return [description]
   */
  bool failed() {
    return fail;
  }

private:
  void guidance(double oldT, double r, double ve, double acc, double vvel);
  void estimation(double dt, double mu, double r, double ve, double acc, double hvel, double vvel);

  double T = -1.0;
  double A,B,C,lastT,pitch;
  double pitchRate = -1.0;
  double cycleTime = 1.0;
  double start = -1.0;
  double end = -1.0;
  double targetRadius = -1.0;
  double targetHoriVel = -1.0;
  bool PEGinit = false;
  bool fail = false;
};

using namespace std;

void Guidance::guidance(double oldT, double r, double ve, double acc, double vvel) {
  double tau = ve/acc;

  double b0 = -ve*log(1-oldT/tau);
  double b1 = b0*tau - ve*oldT;
  double c0 = b0*oldT - b1;
  double c1 = c0*tau - ve*oldT*oldT/2;

  double mbx = -vvel;
  double mby = targetRadius-r-vvel*oldT;
  float det = (b0*c1-c0*b1);
  A = (c1*mbx-b1*mby)/det;
  B = (b0*mby-c0*mbx)/det;
}

void Guidance::estimation(double dt, double mu, double r, double ve, double acc, double hvel, double vvel) {
  double oldT = T;
  double tau = ve/acc;

  if (oldT <= -1.0) {
    oldT = 0.9*tau;
    guidance(oldT, r, ve, acc, vvel);
  }

  // navigation, basis vectors, and additional target conditions
  double h = r*hvel;
  double ht = targetRadius*targetHoriVel;
  double dh = ht - h;
  double rbar = (r + targetRadius)/2;

  // Vehicle performance
  C = (mu/(r*r) - hvel*hvel/r) / acc;
  double fr = A + C;

  double CT = (mu/(targetRadius*targetRadius) - targetHoriVel*targetHoriVel/targetRadius) / (acc / (1 - oldT/tau));
  double frT = A + B*oldT + CT;
  double frdot = (frT - fr)/oldT;
  double ftheta = 1 - fr*fr/2;
  double fthetadot = -fr*frdot;
  double fthetadotdot = -frdot*frdot/2;

  double currentT = oldT-dt;
  double dv = (dh/rbar + ve*currentT*(fthetadot + fthetadotdot*tau) + fthetadotdot*ve*currentT*currentT/2) / (ftheta + fthetadot*tau + fthetadotdot*tau*tau);

  T = tau * (1 - exp(-dv/ve));
  //cout << time << " A = " << A << " B = " << B << " T = " << oldT << " C = " << C << endl;
  // Update until terminal guidance
  if (T >= 7.5) {
    guidance(T, r, ve, acc, vvel);
  }
}

double Guidance::run(double time, double dt, double mu, double r, double ve, double acc, double hvel, double vvel) {
  // Vertical ascent
  if (time < start) pitch = rad(90.0);
  // Pitch program
  else if (time < end) {
    pitch -= pitchRate*dt;
  }
  // peg
  else {
    if (!PEGinit) {
      estimation(0, mu, r, ve, acc, hvel, vvel);
      lastT = T;
      PEGinit = true;
    }

    // Major loop
    if ((lastT - T) >= cycleTime) {
      estimation(cycleTime, mu, r, ve, acc, hvel, vvel);
      lastT = T;
    } else {
      T = T - dt;
    }

    // Minor loop
    double sinPitch = A + (lastT - T)*B + C;
    if (sinPitch <= 1.0 && sinPitch >= -1.0) {
      fail = false;
      pitch = asin(sinPitch);
      return max(0.0, min(dt, T));
    } else {
      fail = true;
      return dt;
    }
  }
  fail = false;
  return dt;
}