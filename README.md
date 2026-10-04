## Overview
Simulation of articulated rigid bodies in Unity (C#), built on a physics engine where rigid bodies are connected through point constraints. The project implements three different ways of integrating and solving the dynamics, making it possible to compare how each one handles constraints in terms of stability and behavior.

## Features
**Rigid bodies**
* Full 6-DoF dynamics: position and orientation, with linear and angular velocity
* Gravity, plus linear and angular damping
* Gyroscopic term in the angular dynamics, since the inertia tensor changes with the orientation

**Point constraints**
* Joins two points, each belonging to a rigid body, or one body and a fixed point in space
* Allows building chains and other articulated structures out of rigid bodies

**Integration and constraint solvers**
* **Symplectic Euler with weak constraints:** constraints behave as penalty springs that add forces to the system. Simple and fast, but the joints stretch and a high stiffness is needed to keep them tight
<p align = "center">
  <img width="538" height="300" alt="Symplectic" src="https://github.com/user-attachments/assets/3b3193ab-f6d4-4846-b13d-4e55cde3c2ef" />
</p>

* **Implicit Euler with weak constraints:** a linear system is solved every step, using the Jacobians of the forces with respect to position and velocity (damping, gyroscopic term and constraint force). Much more stable with stiff constraints and large time steps
<p align = "center">
  <img width="536" height="300" alt="Implicit" src="https://github.com/user-attachments/assets/c85219dd-6875-4e81-a0c9-c17f41552b86" />
</p>

* **Symplectic Euler with strong constraints:** constraints are solved exactly through the constraint vector and its Jacobian, using Lagrange multipliers, so the joints stay together without relying on stiffness
<p align = "center">
  <img width="538" height="300" alt="Symplectic_strong" src="https://github.com/user-attachments/assets/eb584b5a-ae54-4364-8274-3524a55133e0" />
</p>

## Authorship
I only developed the simulation code in `PhysicsManager.cs`, `RigidBody.cs` and `PointConstraint.cs`. The rest of the project was provided by the professor and is included only so the project can run.

## Technologies
* Unity
* C#
* Math.NET Numerics
