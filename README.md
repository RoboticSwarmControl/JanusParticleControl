# Janus Particle Control Simulation
Code and algorithms for Steering Catalytic Janus Particles

Authors: Li Huang and Aaron T. Becker

Email: lhuang21@uh.edu

All rights reserved.

#### Abstract

A catalytic Janus particle is a two-faced particle with one face that reacts with the surrounding medium to produce thrust. By using a permanent magnet core, the particle can be steered. Unlike many current microrobots that are steered and propelled by an external magnetic field, these particles have independent steering and propulsion mechanisms. Janus particles can be manufactured in large numbers. An offset angle between their thrust and magnetization vectors can provide kinematic heterogeneity in a uniform magnetic field, which is the key for controlling multiple microrobots. In the 2D case, only two degrees-of-freedom (DOF) are controllable. We review controllability results in 2D, and then show that interesting things happen in 3D. We provide control laws for steering up to nine DOF, which can be mapped in various ways, including
to control the x; y; z position of three particles, make four particles meet, or reduce the spread of n particles. A hardware implementation is described.

#### MATLAB Simulation Code

- [`Open-loop control Janus spheres using linear programming`](./src/OpenloopControlwithLP.m) 
- [`Closed-loop control using random rotation matrices`](./src/RandomRotation)
- [`Closed-loop control using linear programming`](./src/ClosedloopLP.m)
- [`Closed-loop control using greedy optimal control`](./src/GreedyOptimalControl.m)
- [`Control four spheres to collide`](./src/FourSpheresControl.m)
- [`Control 9 spheres up to 9 DOF `](./src/NineDoF.m)
- [`Steer 10 spheres closer to each other`](./src/TenSpheresControl.m)


#### Demo

Janus sphere simulations of open-loop control using linear programming. 
The goal locations are indicated by green orbits. 
In each figure, all spheres move the same total distance and reach the goal location at the same time. 
A colored line describes the trajectory of each Janus sphere. Black arrows indicate the magnetic moment orientation for the subsequent move.
<p align="center">
<img src="https://github.com/RoboticSwarmControl/JanusParticleControl/blob/master/Media/OpenloopLP.gif" width="600">
</p>

Janus sphere simulations of (closed-loop) linear programming. The demo shows three spheres’ trajectories, and the corresponding sum of squared distance error.
<p align="center">
<img src="https://github.com/RoboticSwarmControl/JanusParticleControl/blob/master/Media/ClosedloopLP.gif" width="600">
</p>

Janus sphere simulations of (closed-loop) optimal greedy control. The demo shows three spheres’ trajectories, and the corresponding sum of squared distance error.
<p align="center">
<img src="https://github.com/RoboticSwarmControl/JanusParticleControl/blob/master/Media/ClosedloopGreedy.gif" width="600">
</p>

Trajectories of four spheres moving towards their mean position till they collide, 
and trajectories of three spheres chasing the fourth sphere till they collide. 
These two cases have the same initial conditions with different target position, but the four spheres collide at the same spot.
<p align="center">
<img src="https://github.com/RoboticSwarmControl/JanusParticleControl/blob/master/Media/4SpheresControl.gif" width="600">
</p>

Closed-loop control simulation of nine Janus spheres with 3D view, xz plane view, and xy plane view. Three groups
of spheres (red, blue, and green) are delivered to x = 0, z = 0, and y = 0 planes respectively
<p align="center">
<img src="https://github.com/RoboticSwarmControl/JanusParticleControl/blob/master/Media/9SpheresControl.gif" width="600">
</p>
