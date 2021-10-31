/*  Parallel Mechanism Control in MuJoCo
 */

#include "mujoco.h"
#include "glfw3.h"
#include "mjgraphics.hpp"

#include "stdio.h"
//#include "stdlib.h"
//#include "string.h"
//#include "math.h"
#include "assert.h"

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

// MuJoCo data structures
mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data
GUI gui;                            // window and UI elements
int simulation_exit;                // 0 to continue, 1 to call to exit loop


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

    printf("Unconstrained Jacobian JacP:\n");
    mju_printMat(jacp, 3, m->nv);
    printf("Unconstrained Jacobian JacR:\n");
    mju_printMat(jacr, 3, m->nv);

    // get the raw constraint jacobian A
    printf("Raw constraint jacobian d->efc_J:\n");
    mju_printMat(d->efc_J, m->njmax, m->nv);

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


    printf("truncated efc_J:\n");
    mju_printMat(efc_J_active, njused, m->nv); // (njused x m->nv)

    // calculate A' (A transpose)
    mjtNum efc_JT[njused*m->nv]; // create variable for JT since in dense mode it is not automatically created
    mju_transpose(efc_JT, efc_J_active, njused, m->nv);

    // variable to store A * A'
    mjtNum AAT[njused*njused];
    mju_mulMatMat(AAT, d->efc_J, efc_JT, njused, m->nv, njused);
    printf("Need to invert this matrix: A * A'\n");
    mju_printMat(AAT, njused, njused);

    // use GSL to invert the matrix TODO eventually need to move memory allocations out of every-timestep code
    gsl_matrix* matrix1 = gsl_matrix_alloc(njused, njused); // want to invert this
    gsl_permutation* matrix2 = gsl_permutation_alloc(njused); // a working matrix
    gsl_matrix* matrix3 = gsl_matrix_alloc(njused, njused); // will store resulting inverted matrix

    // copy matrix to GSL
    for (size_t i=0; i<njused; ++i){ // rows
        for (size_t j=0; j<njused; ++j){ // columns
            gsl_matrix_set(matrix1, i, j, AAT[i*njused + j]);
        }
    }
    
    // actually do the inversion
    int s;
    gsl_linalg_LU_decomp(matrix1, matrix2, &s);
    gsl_linalg_LU_invert(matrix1, matrix2, matrix3);

    // bring data back into mujoco types
    mjtNum inv_AAT[njused*njused];
    mju_copy(inv_AAT, matrix3->data, njused*njused);
    printf("Inverted matrix inv_AAT:\n");
    mju_printMat(inv_AAT, njused, njused);

    // cleanup GSLdata
    gsl_matrix_free(matrix1);
    gsl_permutation_free(matrix2);
    gsl_matrix_free(matrix3);


    // finish calculations of constraint correction matrix


    // combine with unconstrain jacobian - get the constrained jacobian!


    // truncate constrained jacobian according to the directions we want to control and the actuators we have available to do it


    // use the jacobian
    // d->ctrl[]


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
            
        graphics_draw(m, d);
    }
    
    cleanup();

    return 1;
}
