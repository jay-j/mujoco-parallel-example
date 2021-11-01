/*  Parallel Mechanism Control in MuJoCo
 */

#include "mujoco.h"
#include "glfw3.h"
#include "mjgraphics.hpp"

#include "stdio.h"
#include "assert.h"

#include "lapacke.h"

// MuJoCo data structures
mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data
GUI gui;                            // window and UI elements
int simulation_exit;                // 0 to continue, 1 to call to exit loop
int simulation_mode;                // various per the defines
int simulation_step_frames;         // counter to track/limit progress in step mode
#define MODE_SIM 0
#define MODE_INSPECT 1
#define MODE_STEP 2

void mode_switch_inspect(){
    simulation_mode = MODE_INSPECT;
}
void mode_switch_step(){
    simulation_mode = MODE_STEP;
    simulation_step_frames = 0;
}
void mode_switch_resume(){
    simulation_mode = MODE_SIM;
}

// keyboard callback
void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods){
    // backspace: reset simulation
    if( act==GLFW_PRESS && key==GLFW_KEY_BACKSPACE )
    {
        mj_resetData(m, d);
        mj_forward(m, d);
    }

    // time management; quit, pause, step, go
    if (key == GLFW_KEY_ESCAPE){
        simulation_exit = 1;
    }
    if (key == GLFW_KEY_F1){
        mode_switch_inspect();
    }
    if (key == GLFW_KEY_F2){
        mode_switch_step();
    }
    if (key == GLFW_KEY_F3){
        mode_switch_resume();
    }

}


void identity(mjtNum* result, size_t size){ // TODO more efficient or not every loop
    for (size_t i=0; i<size; i++){
        for(size_t j=0; j<size; j++){
            if (i == j){
                result[i*size + j] = 1;
            }
            else{
                result[i*size + j] = 0;
            }
        }
    }
}


// compute matrix pseudo-inverse using the LAPACK function dgelss
// http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_gaa6ed601d0622edcecb90de08d7a218ec.html#gaa6ed601d0622edcecb90de08d7a218ec
// the resulting inverse is stored in the B matrix
void pseudo_inverse(double* result, double* matrix, int nr, int nc){
    //int nr = 3; // (m) number of rows of matrix A
    //int nc = 3; // (n) number of columns of matrix A
    int nrhs = nc; // number of columns of matrices B and X
    int lda = nr; // TODO??!
    int ldb = nr > nc ? nr : nc ;

    // matrix TODO size.. and fill.. and type...]
    double A[lda*nc];
    mju_copy(A, matrix, nr*nc);

    double B[ldb*nrhs];
    // B is the identity matrix of size N
    for (size_t i=0; i<ldb*nrhs; ++i){
        B[i] = 0;
    }
    for (size_t i=0; i<nc; ++i){
        B[i*nc + i] = 1;
    }

    double S[12]; // auto TODO minimum size

    double rcond = -1;
    int rank;

    int matrix_layout = LAPACK_ROW_MAJOR;
    int info = LAPACKE_dgelss(matrix_layout, nr, nc, nrhs, A, lda, B, ldb, S, rcond, &rank);

    if (info < 0){
        printf("If INFO = -i, the i-th argument had an illegal value. info = %d\n", info);
        assert(0);
    }
    if (info > 0){
        printf("The algorithm for computing the SVD failed to converge; if INFO = i, i off-diagonal elements of an intermediate bidiagonal form did not converge to zero. info = %d\n", info);
        assert(0);
    }

    mju_copy(result, B, nr*nc);
}


// THIS IS THE FUNCTION I NEED HELP WITH
void controller(const mjModel* m, mjData* d){
    printf("\n================== in the controller ==================\n");

    int ee_body_id = mj_name2id(m, mjOBJ_BODY, "endeffector");
    mjtNum* ee_pos = &(d->xpos[ee_body_id*3]);

    // get the unconstrained jacobian
    // argument "point" is given in global frame, acting as if welded to the body
    mjtNum jacp[3*m->nv];
    mjtNum jacr[3*m->nv];
    mj_jac(m, d, jacp, jacr, ee_pos, ee_body_id);

    // printf("Unconstrained Jacobian JacP:\n");
    // mju_printMat(jacp, 3, m->nv);
    // printf("Unconstrained Jacobian JacR:\n");
    // mju_printMat(jacr, 3, m->nv);

    // get the raw constraint jacobian A
    // printf("Raw constraint jacobian d->efc_J:\n");
    // mju_printMat(d->efc_J, m->njmax, m->nv);

    // get the active part of d->efc_J, grabbing all rows that have nonzero entries
    mjtNum efc_J_active[m->nv * m->njmax]; // make it bigger than it needs to be
    int njused = 0; // rows used; from 0 to m->njmax
    int index;
    int row_in_use;
    for (size_t i=0; i<m->njmax; ++i){
        row_in_use = 0;

        // copy row every time, overwriting if previous row copied was all zeros
        for (size_t j=0; j<m->nv; ++j ){
            index = i*m->nv + j; // index to check in d->efc_J
            if ((d->efc_J[index] > 1e-6) || (d->efc_J[index] < -1e-6)){
                row_in_use = 1;
            }
            efc_J_active[njused*m->nv + j] = d->efc_J[index];
        }

        if(row_in_use){
            ++njused;
        }
    }

    // printf("truncated efc_J:\n");
    // mju_printMat(efc_J_active, njused, m->nv); // (njused x m->nv)

    // calculate A' (A transpose)
    mjtNum efc_JT[njused*m->nv]; // create variable for JT since in dense mode it is not automatically created
    mju_transpose(efc_JT, efc_J_active, njused, m->nv);

    // variable to store A * A'
    mjtNum AAT[njused*njused];
    mju_mulMatMat(AAT, d->efc_J, efc_JT, njused, m->nv, njused);
    // printf("Need to invert this matrix: A * A'\n");
    // mju_printMat(AAT, njused, njused);

    mjtNum inv_AAT[njused*njused];
    pseudo_inverse(inv_AAT, AAT, njused, njused);
    // printf("Inverted matrix inv_AAT:\n");
    // mju_printMat(inv_AAT, njused, njused);

    // finish calculations of constraint correction matrix
    // invAAT * A : (njused x njused) * (njused x m->nv) = (njused x m->nv)
    mjtNum tmp1[njused*m->nv];
    mju_mulMatMat(tmp1, inv_AAT, efc_J_active, njused, njused, m->nv);

    // A' * tmp1 : (m->nv x njused) * (njused x m->nv) = (m->nv x m->nv)
    mjtNum tmp2[m->nv*m->nv];
    mju_mulMatMat(tmp2, efc_JT, tmp1, m->nv, njused, m->nv);

    mjtNum eye[m->nv*m->nv];
    identity(eye, m->nv);
    mjtNum constraint[m->nv*m->nv];
    mju_sub(constraint, eye, tmp2, m->nv*m->nv);

    // printf("Jacobian Adjustment Matrix:\n");
    // mju_printMat(constraint, m->nv, m->nv);

    // combine with unconstrain jacobian - get the constrained jacobian!
    // J * constraint : (3 x m->nv) * (m->nv x m->nv) = (3 x m->nv)
    mjtNum jacobian_mod[6*m->nv];
    mju_mulMatMat(jacobian_mod, jacp, constraint, 3, m->nv, m->nv);
    mju_mulMatMat(jacobian_mod+(3*m->nv), jacr, constraint, 3, m->nv, m->nv);

    // printf("modified jacoban:\n");
    // mju_printMat(jacobian_mod, 6, m->nv);

    // given actuators, compute which DOF they correspond to
    int actuator_qty = 3; // TODO non-generic
    int actuator_id[actuator_qty]; 
    int dof_actuated[3]; // TODO non-generic
    actuator_id[0] = mj_name2id(m, mjOBJ_ACTUATOR, "act_leg0_joint0");
    dof_actuated[0] = m->jnt_dofadr[m->actuator_trnid[2*actuator_id[0]]];

    actuator_id[1] = mj_name2id(m, mjOBJ_ACTUATOR, "act_leg1_joint0");
    dof_actuated[1] = m->jnt_dofadr[m->actuator_trnid[2*actuator_id[1]]];

    actuator_id[2] = mj_name2id(m, mjOBJ_ACTUATOR, "act_leg2_joint0");
    dof_actuated[2] = m->jnt_dofadr[m->actuator_trnid[2*actuator_id[2]]];

    printf("dof actuated: %d, %d, %d\n", dof_actuated[0], dof_actuated[1], dof_actuated[2]);

    int axis[6];
    axis[0] = 1; // want to control x position
    axis[1] = 1; // want to control y position
    axis[2] = 0;
    axis[3] = 0;
    axis[4] = 0;
    axis[5] = 1; // want to control z rotation

    mjtNum rot[9]; // rotation of the base - for commands relative to a mobile robot base
    identity(rot, 3);
    //mju_copy(rot, &(d->xmat[9*ee_body_id]), 9);

    // truncate constrained jacobian according to the directions we want to control and the actuators we have available to do it
    mjtNum jacobian_active[3*3]; // TODO generic
    for (size_t i=0; i<actuator_qty; ++i){  // loop over actuators
        mjtNum jacobian_column_global[6]; 

        for (size_t j=0; j<6; ++j){
            jacobian_column_global[j] = jacobian_mod[j*m->nv + dof_actuated[i]];
        }

        mjtNum jacobian_column_local[6];
        mju_transformSpatial(jacobian_column_local, jacobian_column_global, 0, ee_pos, ee_pos, rot); // TODO point old..

        // copy data in jacobian_active - for only the axis we care about
        int k = 0;
        for(size_t j=0; j<6; j++){
            if (axis[j] == 1){
                jacobian_active[k*actuator_qty + i] = jacobian_column_local[j];
                ++k;
            }
        }
    }

    printf("Now reduced to the active jacobian:\n");
    mju_printMat(jacobian_active, 3, 3);

    mjtNum jacobian_active_T[3*3]; // TODO generic size
    mju_transpose(jacobian_active_T, jacobian_active, 3, 3); // TODO generic size

    // TODO a bunch of relative position error calculating stuff!!!!
    mjtNum state_target[3];
    state_target[0] = 0.1*cos(0.05*d->time);
    state_target[1] = 0.1*sin(0.05*d->time);
    state_target[2] = 0;

    mjtNum state_current[3];
    state_current[0] = d->xpos[3*ee_body_id];
    state_current[1] = d->xpos[3*ee_body_id+1];

    mjtNum state_target_quat[4];
    state_target_quat[0] = 1;
    state_target_quat[1] = 0;
    state_target_quat[2] = 0;
    state_target_quat[3] = 0;

    mjtNum state_error_quat[3];
    mju_subQuat(state_error_quat, state_target_quat, &(d->xquat[4*ee_body_id]));
    printf("quat error: %lf, %lf, %lf\n", state_error_quat[0], state_error_quat[1], state_error_quat[2]);

    state_current[2] = 2*state_error_quat[2]; // TODO is this what I want?

    mjtNum desired_velocity_cart[3];
    mju_sub3(desired_velocity_cart, state_target, state_current); // TODO better error for angles, generic size
    desired_velocity_cart[2] *= 0.05; // TODO is there a more sensible way to control angular axis proportions? auto based on type
    printf("Desired Velocity Cart: %lf, %lf, %lf\n", desired_velocity_cart[0], desired_velocity_cart[1], desired_velocity_cart[2]);

    mjtNum desired_velocity_joint[actuator_qty];
    mju_mulMatVec(desired_velocity_joint, jacobian_active_T, desired_velocity_cart, actuator_qty, 3); // todo generic axis qty

    // use the jacobian, very simple proportional controller
    for (size_t i=0; i<actuator_qty; ++i){
        // get current location of actuator
        mjtNum loc = d->actuator_length[actuator_id[i]];

        // increment using desired joint velocity calc
        d->ctrl[actuator_id[i]] = loc + 5*desired_velocity_joint[i];
        printf("cmd: %lf  vel: %lf\n", d->ctrl[actuator_id[i]],  desired_velocity_joint[i]);
    }
}


int main(int argc, const char** argv){

    // load model, setup data
    char error[1000] = "Could not load binary model";
    m = mj_loadXML("scene3.xml", 0, error, 1000);
    if (!m){
        mju_error_s("Load model error: %s", error);
    }
    d = mj_makeData(m);

    // check settings we think are good for these kind of jacobian manipulations
    assert(m->opt.jacobian == mjJAC_DENSE);
    assert(m->opt.cone == mjCONE_ELLIPTIC);

    simulation_exit = 0;

    // initialize graphics, send the keyboard callback function to the window
    graphics_init(&gui, keyboard); 

    // run main loop, target real-time simulation and 60 fps rendering
    while( (!glfwWindowShouldClose(gui.window)) && (simulation_exit == 0) ){

        if ((simulation_mode == MODE_SIM) || (simulation_mode == MODE_STEP)){

            mjtNum simstart = d->time;
            while( d->time - simstart < 1.0/60.0 ){
                // calculate as much as possible
                mj_kinematics(m, d);
                mj_forward(m, d);

                // run the controller on the parallel mechanism
                controller(m, d);

                // finally, advance time
                mj_step(m, d);
            }

            // if in step mode, pause after a couple frames
            simulation_step_frames++;
            if (simulation_mode == MODE_STEP){
                if (simulation_step_frames > 1){
                    simulation_mode = MODE_INSPECT;
                }
            }
        }
        graphics_draw(m, d);
    }
    
    cleanup();

    return 1;
}
