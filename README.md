# Immersed Boundary Methods for Rapid Prototyping (IBM-RP)

A MATLAB-based computational fluid dynamics (CFD) solver for investigating immersed boundary methods. This solver was developed as part of my Master's project which investigated how different boundary reconstruction methods affect the resulting solution.

When I'm able to, I intend to make changes to make the code more readable/clear as well as add clarification in the README of how to use the GUI. Currently it is presented here "as-is", i.e. with minimal modification, to allow for anyone who may be interested to look through.

## Example Results

![Example flow visualisation gif showing temperature around a 2D immersed cylinder](/RW-05-IBM-for-Rapid-Prototyping/Icons/IBM_FLOW_VIS.gif)
_Flow visualisation of temperature around an immersed (2D) cylinder._

![Image of IBM RP graphical user interface](/RW-05-IBM-for-Rapid-Prototyping/Icons/CD_GUI.png)
_IBM-RP GUI during execution of a simulation._

## Solver Description

The solver makes use of Immersed Boundary Methods. This means that the geometry is represented by forcing functions which add terms to the governing equations which aim to replicate the effect of the geometry on the flow.

The solver is capable of estimating 2D flow fields for various shapes on a rectilinear Cartesian grid. It uses the following numerical methods:

* Governing Equations - Compressible Navier-Stokes
* Spatial Discretisation - Cell-centred finite volume
* Reconstruction Scheme - MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws)
  * with Venkatakrishnan limiter
* Convective Flux Scheme - HLLC (Harten–Lax–van Leer–Contact)
* Time Integration - Explicit second-order Runge-Kutta

The overall algorithm is described in the flow chart below.

![Algorithmic flowchart describing IBM RP solver execution](/RW-05-IBM-for-Rapid-Prototyping/Icons/Algorithm_Flow_Chart.png)
_Algorithmic flowchart describing IBM RP solver execution._

### Functions

The main_loop function executes the RK2 time looping, within which the flux computations occur for each timestep.

The various computations are performed by several different helper functions. For example:

* IBMForcing
  * Computes the required forcing to enforce the boundary condition within the timestep.
* HLLC
  * Computes the fluxes at cell boundaries (for a given direction) using the HLLC approximate Riemann solver.

## Using the Solver

### Requirements

Requires MATLAB 2023b or later with the following toolboxes:

* App Designer Toolbox
* PDE Toolbox
* Image Processing Toolbox
* Partial Differential Equation Toolbox

### How to run the solver:

Simulations can be set up via the GUI by launching the CFD_Controller.mlapp in App Designer.

## To-Dos

* Add further clarification explaining how to set-up and run simulations via the GUI.
* Create a standalone executable to remove the need for specific toolboxes.
* Improve code clarity and remove unnecessary code.

---
<font size=2>Note on use of GenAI:

_An initial simple 2D compressible solver was generated using AI to assist in the early stages of this project. It has since been modified significantly from the AI generated code. The original AI generated code is contained in the directory "GPTOriginal". 2 further files can be found here which were also AI generated but not directly used/modified but used as a guide to develop certain function._
