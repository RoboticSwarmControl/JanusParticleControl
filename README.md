# Janus Particle Control Simulation
Code and algorithms for Steering Catalytic Janus Particles

### Demo

Janus sphere simulations of open-loop control using linear programming. 
The top two figures each have two Janus spheres, and the bottom two figures each have three spheres.
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
