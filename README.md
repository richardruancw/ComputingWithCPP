# CME212 Submission Repository

### Parallel v.s Sequantial Test for HW 4

The following are results of average runing time of each call of symp_euler_step:

#### Set up: 

* Hardware: virtual machine with 4cpu assigned

* mass_spring file: 

Use ``GravityForce()``, ``MassSpringForce()``, ``DampingForce()``, and ``PlainConstraint()``, ``SphereConstraint()``, ``SphereRemoveConstraint``.

Simulation parameters follow from HW4 suggestions

#### Results
| Grids    | Poly        | Average running time each step  |
| ---------|:-----------:| -------------------------------:|
| grid 1   | sequantial  | 0.000249 |
| grid 1   | parallel    | 0.001705 |
| grid 2   | sequantial  | 0.000732 |
| grid 2   | parallel    | 0.002851 |
| grid 3   | sequantial  | 0.003767 |
| grid 3   | parallel    | 0.011301 |
| grid 4   | sequantial  | 0.241682 |
| grid 4   | parallel    | 0.091106 |
