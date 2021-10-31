#ifndef MJGRAPHICS_H
#define MJGRAPHICS_H

#include "mujoco.h"
#include "glfw3.h"

struct GUI {
    mjvCamera cam;                      // abstract camera
    mjvOption opt;                      // visualization options
    mjvScene scn;                       // abstract scene
    mjrContext con;                     // custom GPU context

    // mouse interaction
    bool button_left = false;
    bool button_middle = false;
    bool button_right =  false;
    double lastx = 0;
    double lasty = 0;

    GLFWwindow* window;
};

// initialize the gui window
void graphics_init(GUI* gui, void (*keyboard)( GLFWwindow* window, int key, int scancode, int act, int mods));

// draw one frame
void graphics_draw(const mjModel* m, mjData* d);

// cleanup of graphics and mujoco model variables
void cleanup();

// mouse stuff
void mouse_button(GLFWwindow* window, int button, int act, int mods);
void mouse_move(GLFWwindow* window, double xpos, double ypos);
void scroll(GLFWwindow* window, double xoffset, double yoffset);


#endif