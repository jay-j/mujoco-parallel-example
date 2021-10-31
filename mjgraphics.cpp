#include "mjgraphics.hpp"

extern mjModel* m;
extern mjData* d;


void graphics_init(GUI* gui, void (*keyboard)( GLFWwindow* window, int key, int scancode, int act, int mods)){

    // init GLFW
    if( !glfwInit() ){
        mju_error("Could not initialize GLFW");
    }

    // create window, make OpenGL context current, request v-sync
    gui->window = glfwCreateWindow(1200, 900, "MuJoCo Simulation", NULL, NULL);
    glfwMakeContextCurrent(gui->window);
    glfwSwapInterval(1);

    // initialize visualization data structures
    mjv_defaultCamera(&gui->cam);
    mjv_defaultOption(&gui->opt);
    mjv_defaultScene(&gui->scn);
    mjr_defaultContext(&gui->con);

    // create scene and context
    mjv_makeScene(m, &gui->scn, 2000);
    mjr_makeContext(m, &gui->con, mjFONTSCALE_150);

    // install GLFW mouse and keyboard callbacks
    glfwSetKeyCallback(gui->window, keyboard);
    glfwSetCursorPosCallback(gui->window, mouse_move);
    glfwSetMouseButtonCallback(gui->window, mouse_button);
    glfwSetScrollCallback(gui->window, scroll);
}


// draw one frame
void graphics_draw(const mjModel* m, mjData* d){
    extern GUI gui;

    // get framebuffer viewport
    mjrRect viewport = {0, 0, 0, 0};
    glfwGetFramebufferSize(gui.window, &viewport.width, &viewport.height);

    // update scene and render
    mjv_updateScene(m, d, &gui.opt, NULL, &gui.cam, mjCAT_ALL, &gui.scn);
    mjr_render(viewport, &gui.scn, &gui.con);

    // swap OpenGL buffers (blocking call due to v-sync)
    glfwSwapBuffers(gui.window);

    // process pending GUI events, call GLFW callbacks
    glfwPollEvents();
}


// cleanup of graphics and mujoco model variables
void cleanup(){
    extern GUI gui;
    mjv_freeScene(&gui.scn);
    mjr_freeContext(&gui.con);

    mj_deleteData(d);
    mj_deleteModel(m);

    // terminate GLFW (crashes with Linux NVidia drivers)
    #if defined(__APPLE__) || defined(_WIN32)
        glfwTerminate();
    #endif
}


// mouse button callback
void mouse_button(GLFWwindow* window, int button, int act, int mods){
    extern GUI gui;

    // update button state
    gui.button_left =   (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS);
    gui.button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS);
    gui.button_right =  (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS);

    // update mouse position
    glfwGetCursorPos(window, &gui.lastx, &gui.lasty);
}


// mouse move callback
void mouse_move(GLFWwindow* window, double xpos, double ypos){
    extern GUI gui;

    // no buttons down: nothing to do
    if( !gui.button_left && !gui.button_middle && !gui.button_right )
        return;

    // compute mouse displacement, save
    double dx = xpos - gui.lastx;
    double dy = ypos - gui.lasty;
    gui.lastx = xpos;
    gui.lasty = ypos;

    // get current window size
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // get shift key state
    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT)==GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS);

    // determine action based on mouse button
    mjtMouse action;
    if( gui.button_right )
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    else if( gui.button_left )
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    else
        action = mjMOUSE_ZOOM;

    // move camera
    mjv_moveCamera(m, action, dx/height, dy/height, &gui.scn, &gui.cam);
}


// scroll callback
void scroll(GLFWwindow* window, double xoffset, double yoffset){
    extern GUI gui;

    // emulate vertical mouse motion = 5% of window height
    mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05*yoffset, &gui.scn, &gui.cam);
}
