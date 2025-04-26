# AAE590ACA-Stochastic-SCP-Rocket-Landing

Stochastic Sequential Convex Programming using the Penaliezd Trust Region (PTR) method for robust control of a rocket landing.

## Level 1 Roadmap
- [x] [4/6] Define the dynamics of the planar system and the constraints. 
- [x] [4/6] Convexify nonconvex elements of deterministic problem where appropriate
- [x] [4/6] Solve deterministic optimal control problem with SCP using [SCPToolbox.jl](https://github.com/UW-ACL/SCPToolbox.jl/tree/master) to get a baseline solution to help us implement SCP ourselves
- [Milestone:] Problem formulation
- [x] [4/12] Implement the LTV Kalman filter from equations in Oguri paper for accurate state estimation throughout the trajectory.
- [x] [4/12] Convexify stochastic elements of problem
- [Milestone:] Convex sub-problem formulation
- [x] [4/12] Implement SCP and solve the baseline deterministic optimal control problem
- [ ] [4/16] Solve stochastic optimal control problem (convexified and turned into an approximate deterministic one) with stochastic SCP
- [Milestone:] Basic stochastic implementation
- [x] [4/16] Create a Monte Carlo simulation to help analyze the constraint satisfaction
## Level 2 Roadmap
- [ ] [4/19] Implement Gates control error model and unmodeled external disturbances
- [Milestone:] Full stochastic implementation
- [ ] [4/20] Implement FOH
- [ ] [4/20] Add final time as parameter to be optimized
- [ ] [4/23] Sensitivity analysis to see what sources of uncertainty are most influential
- [Milestone:] Full performance analysis
## Level 3 Roadmap
- [X] [4/26] 6DoF dynamics and optimal control problem defined (roll controlled by cold gas thrusters) 
- [X] [4/28] 6DoF baseline solution of deterministic optimal control problem with SCP
- [ ] [4/28] 6DoF stochastic optimal control problem solution with stochastic SCP
- [ ] [4/28] Assess accuracy of uncertainty propagation
- [Milestone:] Full rocket landing problem
- [ ] [4/28] Higher order uncertainty propagation methods if needed
