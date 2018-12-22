# launchsim

Rocket launch simulator (2D only). Plug in your guidance parameters, and watch whether your rocket either inserts itself into a stable orbit or fails miserably.

## What does it do!?

It's a simulation of how a rocket launch would behave, coupled with guidance algorithms.
The guidance algorithms starts with a standard straight-up pitch hold to build up vertical speed, then a pitch program runs to steer the rocket at a constant rate (so the angle of attack of the rocket doesnt endanger its structural integrity) and once the rocket can pitch freely, a Powered Explicit Guidance program (PEG) takes care of inserting the rocket in a circular orbit at the target altitude. 
The PEG algorithm is explained [here](https://www.orbiterwiki.org/wiki/Powered_Explicit_Guidance) and was made by smart people at NASA.

## Usage

For now there is only one rocket file : `skylab_saturnv.txt` which describes the Skylab space station mounted atop a two-stage Saturn V rocket.

You can add your own rockets using this sample file, its structure should be as following :
* A `#` symbol means the rest of the line is ignored
* Each line describes a stage, the topmost line is the topmost stage (or payload), the next is the stage below...
* A stage can either have no engines : in this case the only value on the line is its mass
* A stage with engines needs 5 or 6 values on the line : dry mass, wet mass, thrust ISP of the engine(s) in vacuum and at sea level and optionaly an atmospheric drag coefficient
* The ISP is in seconds
* The dry mass and the wet mass need to be the same unit (but can be whatever, kg, ton, slugs, lb) for ALL the stages
* Thrust needs to be consistent with the previous mass unit (UNIT.m/s², so N if mass is in kg, kN if mass is in tons....)
* The drag coefficient has to be tweaked as it is a completely arbitrary value
* No drag coefficient means no atmospheric drag at all.

Let's analyze our sample file : 
```
# SATURN V / SKYLAB (tons, kN)
88.6                       #SKYLAB
39 491 5116 421 200        #S-II
135.2 2286 38703 304 265 0.08   #S-IC
```

The mass unit is the ton, and the unit of force is the kilonewton. So we can deduce that :
* The topmost stage (SKYLAB) doesn't have an engine, and weighs 88.6 tons
* The stage below it (S-II) has engines, weighs 39 tons dry, 491 tons wet, produces 5116kN of thrust and its ISP is 421s in vacuum and 200s at sea level.
* The lowermost stage (S-IC) has engines, weighs 135.2 tons dry, 2286 tons wet, produces 38703kN of thrust, its ISP is 304s in vacuum and 265 at sea level, and it has a drag coefficient of 0.08.

You can find all those numbers for existing rockets on the internet.

`launchsim` has default guidance parameters for inserting the sample rocket in its historical orbit, but you can pass parameters in other cases : 
```./launchsim filename [pitchProgramStart pitchProgramEnd pitchRate targetAltitude]```
* filename is the name of the rocket file obviously, duh
* pitchProgramStart is the time for which the rocket stands straight until the start of the pitch program
* pitchProgramEnd is the time at which the pitch program ends and Powered Explicit Guidance starts
* pitchRate is the degrees the rocket steers per second during the pitch program
* targetAltitude is the altitude in meters from sea level at which the rocket will try to insert itself

Here's an excerpt of the output for successful parameters :
```
TIME(s) DOWNRANGE(km) ALT(km) HVEL(m/s) VVEL(m/s) PITCH(deg) PERI(km) APO(km) e
1 4.47784e-19 0.00198063 8.28887e-16 3.67321 90 -6357 0.00269246 1
2 1.76006e-18 0.00781752 1.67825e-15 7.48351 90 -6357 0.0107139 1
3 3.94086e-18 0.0175765 2.53157e-15 11.3587 90 -6357 0.0242054 1
4 6.99425e-18 0.0313242 3.38891e-15 15.2993 90 -6357 0.0433109 1
5 1.09243e-17 0.0491276 4.25028e-15 19.3061 90 -6357 0.0681774 1
...
519 1443.07 437.969 7550.36 19.6241 -10.0101 67.1274 438.753 0.0281112
520 1450.3 437.985 7586.79 12.9236 -10.1436 189.441 438.494 0.0186669
521 1457.56 437.994 7623.54 6.16223 -10.2765 315.752 438.226 0.00909373
522 1464.61 437.997 7660.13 -0.571114 -10.4089 437.932 441.945 0.000295239
SHUTDOWN
Insertion in a 437.932kmx441.945km (e=0.000295239) orbit in 522.207s
Remaining propellant mass : 10.2061
```

Use your preferred software to make pretty graphs and stuff.
* TIME is the time since launch
* DOWNRANGE is the distance on the ground from the pad to the rocket
* ALT is the altitude of the rocket from sea level
* HVEL is the tangential velocity (parallel to ground)
* VVEL is the radial velocity (perpendicular to ground)
* PITCH is the pitch of the rocket relative to the ground (90° = straight up, 0 = parallel to ground)
* PERI is the periapsis of the orbit relative to sea level (lowest altitude)
* APO is the apoapsis of the orbit relative to sea level (highest altitude)
* e is the eccentricity of the orbit (closer to 0 = more circular)

If the resulting orbit has a periapsis or apoapsis very much below the target altitude, it's a fail. The 'game' is to find the parameters that result in the least fuel wasted
