# Stokes-FEM-Jetflow-Solver
This project implements a finite element model (FEM) to simulate incompressible viscous flow through an idealised stenosed aortic valve geometry at peak systole.

The governing incompressible Navier–Stokes equations are solved using a mixed velocity–pressure finite element formulation with Newton–Raphson iteration. The model resolves full velocity and pressure fields and computes the mean transvalvular pressure drop as a physics-based analogue of the clinically measured pressure gradient.

The project was developed in the context of investigating pressure loss mechanisms in aortic valve stenosis and evaluating the suitability of rigid-wall FEM models for pressure-focused cardiovascular flow analysis.

Model Overview
- Governing equations: Incompressible Navier–Stokes
- Flow assumption: Steady-state peak systolic flow
- Fluid model: Newtonian, incompressible
- Geometry: Idealised 2D stenosed valve domain
- Discretisation: Unstructured triangular mesh
- Mixed formulation: Quadratic basis functions for velocity + Linear basis functions for pressure
- Nonlinear solver: Newton–Raphson method
- Convergence criterion: residual norm

This project demonstrates 
- PDE modelling of cardiovascular flow
- Weak formulation derivation and finite element discretisation
- Mixed velocity-pressure FEM implementation
- Assembly of residual vectors and block Jacobians
- Newton–Raphson nonlinear solution procedure
- Sparse linear system solving
- Mesh refinement and convergence analysis
- Post-processing of derived hemodynamic quantities
- Application of computational modelling to biomedical problems
