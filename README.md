# mujoco-parallel-example
how can a parallel mechanism be controlled in mujoco?

# General Problem Description
## Given
A parallel mechanism, where the kinematic tree between "ground" and the end effector includes loops. MuJoCo uses a serial kinematic tree, so loops are formed using the `equality/connect` constraint.

## Find
A constrained jacobian which maps from actuator (joint) velocity to end effector (cartesian) velocity, taking into account the equality/connect constraints (which define the parallel mechanism).

## Approach
0. Use dense jacobian representation and elliptic solver.
1. Get the unconstrained jacobian for the body with `mj_jac()`.
2. Get the constraint jacobian `d->efc_J`, limit to rows with at least one nonzero entry (constraint is active in some way).
3. Follow the formula `x = J * (I - A' * inv(A * A') * A) * v`.
   - Calculate the transpose manually since it is apparently not done automatically.
   - Use LAPACKE to calculate the matrix pseudo inverse.
4. Extract the portion of the jacobian that we are actually interested in, perform coordinate transformations.
5. Apply the jacobian.

# Example Problem 1
A 3DOF planar mechanism. The end effector is connected to ground via three legs. Each 3R leg has two passive joints and one actuated joint.

# Compiling
## Setup
- Dependent upon LAPACKE (`liblapacke-dev` on my Linux Mint system) package to perform matrix pseudo-inversion.
- system? gcc, make, gl, glew..
- Includes a copy of MuJoCo 2.10 that it references to compile. 
## Compile
Type `make`. 

## Run
Run the executable `./main.bin`

# Notes
- Reference [Jacobian for Contact](https://roboti.us/forum/index.php?threads/jacobian-for-contact.3482/). Use "dense" jacobian representation and elliptic solver to make manual jacobian manipulations easier.
  - This is checked by assertions in `main.cpp`.
- Reference [Programming, data layout and buffer allocation](https://mujoco.readthedocs.io/en/latest/programming.html). The transpose `mjData.efc_JT` is **only** computed if the sparse jacobian is used. 
- Reference [How to get full submatrix of sparse jacobian?](https://roboti.us/forum/index.php?threads/how-to-get-full-submatrix-of-sparse-jacobian.4068). Different constraints have different dimensions at runtime, scan the vector `efc_type` and look for `mjCNSTR_EQUALITY` in it.
  - Use if I need to extract just a portion of the constrained jacobian?
- Reference [How are constraints and efc_J initialized after initializing mjData and mjModel?](https://roboti.us/forum/index.php?threads/how-are-constraints-and-efc_j-initialized-after-initializing-mjdata-and-mjmodel.4100). `mjData.nefc` is the current number of constraints. Run `mj_forward()` to calculate constrained Jacobian `efc_J` (supported by API reference).
- Reference [End effector jacobian for model with equality constraints](https://roboti.us/forum/index.php?threads/end-effector-jacobian-for-model-with-equality-constraints.3478/#post-4049). Jacobian from `mj_jac()` is computed for the serial-topology kinematic tree. This can be combined with `mjData.efc_J` to compute a constrained Jacobian. Let `x` be the cartesian velocity, `v` be the joint velocities, such that the unconstrained map is `x=Jv`. Then `x = J * (I - A' * inv(A * A') * A) * v`.
  - `A` is the linearized equality constraint `Av=0`.. but <mark> how to get `A` from `mjData.efc_J`? </mark>
  - `mjData.efc_J` is `njmax` (maximum allowable constraints) by `nv` (number of DOF). <mark> So expect I need use the submatrix with `mjData.nefc` rows (only active constraints; ignore currently un-used areas of the array). Is `mjData.nefc` the correct number of rows? Or do I need to search `mjData.efc_J_rownnz` to determine that? Is `mjData.efc_J` zeroed between steps? </mark>

- Reference [Programming, jacobians](https://mujoco.readthedocs.io/en/latest/programming.html#jacobians). `mjData.efc_J` is a jacobian matrix of all scalar constraint violations. 
  - So values of joint velocity `v` which satisfy `mjData.efc_J * v = 0` will maintain the current amount of constraint violations (ideally zero). So conclude that `A` from above is simply the active portion of `mjData.efc_J`.
- Attempted to use the subset of `mjData.efc_J` which is simply `mjData.nefc` by `m->nv`. This was unsuccessful. Resulted in a singular (non invertible) matrix. Looking at the matrix `A*A'` to be inverted, it has several columns and rows of all zeros in the example case, likely due to the planar scenario being described. <mark>Is `A*A'` expected to be perfectly invertable, or is the expectation I need to use a pseudo inverse? </mark>
- Attempted to truncate `mjData.efc_J` by only using rows which `mjData.efc_J_rownnz` reports as having some nonzeros. Not successful: `mjData.efc_J_rownnz=0` for all rows. 

