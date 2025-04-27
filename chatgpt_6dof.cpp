#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept> // For std::runtime_error

// Constants
const double mass = 1.0;       // Mass (kg)
const double dt = 0.01;         // Time step (s)
const int steps = 1000;         // Number of simulation steps
const double gravity = 9.81;    // Gravity (m/s^2)

// State structure
struct State {
    double u, v, w;
    double p, q, r;
    double phi, theta, psi;
    double x, y, z;
};

// Separate Force and Moment structures
struct Force {
    double x, y, z;
};

struct Moment {
    double x, y, z;
};

struct Inertia_Tensor
{
    double xx, xy, xz;
    double yx, yy, yz;
    double zx, zy, zz;
};


State state = {10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
Force force = {0.0, 0.0, 0.0};
Moment moment = {0.0, 0.0, 0.0};
Inertia_Tensor inertia = {
    0.1, 0.0, 0.0,
    0.0, 0.1, 0.0,
    0.0, 0.0, 0.2
};


// Invert 3x3 inertia tensor once
Inertia_Tensor invert_inertia(const Inertia_Tensor& I){
    // Calculate determinant
    double det = I.xx * (I.yy * I.zz - I.yz * I.zy)
               - I.xy * (I.yx * I.zz - I.yz * I.zx)
               + I.xz * (I.yx * I.zy - I.yy * I.zx);

    if (std::abs(det) < 1e-12) {
        throw std::runtime_error("Inertia matrix is singular, cannot invert.");
    }

    Inertia_Tensor inv;

    inv.xx =  (I.yy * I.zz - I.yz * I.zy) / det;
    inv.xy = -(I.xy * I.zz - I.xz * I.zy) / det;
    inv.xz =  (I.xy * I.yz - I.xz * I.yy) / det;

    inv.yx = -(I.yx * I.zz - I.yz * I.zx) / det;
    inv.yy =  (I.xx * I.zz - I.xz * I.zx) / det;
    inv.yz = -(I.xx * I.yz - I.xy * I.zx) / det;

    inv.zx =  (I.yx * I.zy - I.yy * I.zx) / det;
    inv.zy = -(I.xx * I.zy - I.xy * I.zx) / det;
    inv.zz =  (I.xx * I.yy - I.xy * I.yx) / det;

    return inv;
}

Inertia_Tensor inv_inertia = invert_inertia(inertia);

// External Forces and Moments (time-varying based on user command)
void compute_forces_and_moments(Force& f, Moment& m, double current_time) {
    f.x = 0.0;
    f.y = 0.0;
    f.z = 0.0 * mass * gravity;

    m.x = 0.0;
    m.y = 0.0;
    m.z = 0.0;

    if (current_time > 2.0 && current_time <= 4.0) {
        f.y = 50.0; // Apply lateral force between 2s and 4s
    }

    if (current_time > 5.0 && current_time <= 6.0) {
        m.z = 10.0; // Apply yaw moment between 5s and 6s
    }
}

// Derivative function
State compute_derivatives(const State& state, const Force& f, const Moment& m) {
    State dstate;

    double cphi = cos(state.phi), sphi = sin(state.phi);
    double ctheta = cos(state.theta), stheta = sin(state.theta);
    double cpsi = cos(state.psi), spsi = sin(state.psi);

    // Translational accelerations
    dstate.u = f.x / mass + state.r * state.v - state.q * state.w;
    dstate.v = f.y / mass + state.p * state.w - state.r * state.u;
    dstate.w = f.z / mass + state.q * state.u - state.p * state.v;

    // Compute Rotational accelerations

    // Compute I * omega
    double Ix_omega = inertia.xx * state.p + inertia.xy * state.q + inertia.xz * state.r;
    double Iy_omega = inertia.yx * state.p + inertia.yy * state.q + inertia.yz * state.r;
    double Iz_omega = inertia.zx * state.p + inertia.zy * state.q + inertia.zz * state.r;
    
    // Compute omega x (I * omega)
    double cross_x = state.q * Iz_omega - state.r * Iy_omega;
    double cross_y = state.r * Ix_omega - state.p * Iz_omega;
    double cross_z = state.p * Iy_omega - state.q * Ix_omega;
    
    // Compute right-hand side (M - omega x (I*omega))
    double rhs_x = moment.x - cross_x;
    double rhs_y = moment.y - cross_y;
    double rhs_z = moment.z - cross_z;
    
    // Multiply Ineverse Inertia Tensor across right-hand side to obtain omega_dot
    dstate.p = inv_inertia.xx * rhs_x + inv_inertia.xy * rhs_y + inv_inertia.xz * rhs_z;
    dstate.q = inv_inertia.yx * rhs_x + inv_inertia.yy * rhs_y + inv_inertia.yz * rhs_z;
    dstate.r = inv_inertia.zx * rhs_x + inv_inertia.zy * rhs_y + inv_inertia.zz * rhs_z;

    // Euler angle rates
    dstate.phi = state.p + state.q * sin(state.phi) * tan(state.theta) + state.r * cos(state.phi) * tan(state.theta);
    dstate.theta = state.q * cos(state.phi) - state.r * sin(state.phi);
    dstate.psi = state.q * sin(state.phi) / cos(state.theta) + state.r * cos(state.phi) / cos(state.theta);

    // Position derivatives (in inertial frame)
    dstate.x = ctheta * cpsi * state.u + (sphi * stheta * cpsi - cphi * spsi) * state.v + (cphi * stheta * cpsi + sphi * spsi) * state.w;
    dstate.y = ctheta * spsi * state.u + (sphi * stheta * spsi + cphi * cpsi) * state.v + (cphi * stheta * spsi - sphi * cpsi) * state.w;
    dstate.z = -stheta * state.u + sphi * ctheta * state.v + cphi * ctheta * state.w;

    return dstate;
}

// Runge-Kutta 4th order integration step
void rk4_step(State& state, const Force& f, const Moment& m) {
    State k1 = compute_derivatives(state, f, m);

    State temp_state = state;
    temp_state.u += 0.5 * dt * k1.u; temp_state.v += 0.5 * dt * k1.v; temp_state.w += 0.5 * dt * k1.w;
    temp_state.p += 0.5 * dt * k1.p; temp_state.q += 0.5 * dt * k1.q; temp_state.r += 0.5 * dt * k1.r;
    temp_state.phi += 0.5 * dt * k1.phi; temp_state.theta += 0.5 * dt * k1.theta; temp_state.psi += 0.5 * dt * k1.psi;
    temp_state.x += 0.5 * dt * k1.x; temp_state.y += 0.5 * dt * k1.y; temp_state.z += 0.5 * dt * k1.z;

    State k2 = compute_derivatives(temp_state, f, m);

    temp_state = state;
    temp_state.u += 0.5 * dt * k2.u; temp_state.v += 0.5 * dt * k2.v; temp_state.w += 0.5 * dt * k2.w;
    temp_state.p += 0.5 * dt * k2.p; temp_state.q += 0.5 * dt * k2.q; temp_state.r += 0.5 * dt * k2.r;
    temp_state.phi += 0.5 * dt * k2.phi; temp_state.theta += 0.5 * dt * k2.theta; temp_state.psi += 0.5 * dt * k2.psi;
    temp_state.x += 0.5 * dt * k2.x; temp_state.y += 0.5 * dt * k2.y; temp_state.z += 0.5 * dt * k2.z;

    State k3 = compute_derivatives(temp_state, f, m);

    temp_state = state;
    temp_state.u += dt * k3.u; temp_state.v += dt * k3.v; temp_state.w += dt * k3.w;
    temp_state.p += dt * k3.p; temp_state.q += dt * k3.q; temp_state.r += dt * k3.r;
    temp_state.phi += dt * k3.phi; temp_state.theta += dt * k3.theta; temp_state.psi += dt * k3.psi;
    temp_state.x += dt * k3.x; temp_state.y += dt * k3.y; temp_state.z += dt * k3.z;

    State k4 = compute_derivatives(temp_state, f, m);

    state.u += (dt / 6.0) * (k1.u + 2 * k2.u + 2 * k3.u + k4.u);
    state.v += (dt / 6.0) * (k1.v + 2 * k2.v + 2 * k3.v + k4.v);
    state.w += (dt / 6.0) * (k1.w + 2 * k2.w + 2 * k3.w + k4.w);
    state.p += (dt / 6.0) * (k1.p + 2 * k2.p + 2 * k3.p + k4.p);
    state.q += (dt / 6.0) * (k1.q + 2 * k2.q + 2 * k3.q + k4.q);
    state.r += (dt / 6.0) * (k1.r + 2 * k2.r + 2 * k3.r + k4.r);
    state.phi += (dt / 6.0) * (k1.phi + 2 * k2.phi + 2 * k3.phi + k4.phi);
    state.theta += (dt / 6.0) * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta);
    state.psi += (dt / 6.0) * (k1.psi + 2 * k2.psi + 2 * k3.psi + k4.psi);
    state.x += (dt / 6.0) * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
    state.y += (dt / 6.0) * (k1.y + 2 * k2.y + 2 * k3.y + k4.y);
    state.z += (dt / 6.0) * (k1.z + 2 * k2.z + 2 * k3.z + k4.z);
}

int main() {
    std::ofstream file("rigidbody6dof_output.csv");
    file << "time,u,v,w,p,q,r,phi,theta,psi,x,y,z,fx,fy,fz,mx,my,mz\n";

    for (int i = 0; i < steps; ++i) {
        double current_time = i * dt;
        compute_forces_and_moments(force, moment, current_time);
        rk4_step(state, force, moment);

        file << current_time << ","
             << state.u << "," << state.v << "," << state.w << ","
             << state.p << "," << state.q << "," << state.r << ","
             << state.phi << "," << state.theta << "," << state.psi << ","
             << state.x << "," << state.y << "," << state.z << ","
             << force.x << "," << force.y << "," << force.z << ","
             << moment.x << "," << moment.y << "," << moment.z << "\n";
    }

    file.close();
    std::cout << "6-DOF rigid body simulation with RK4 and time-varying force/moment inputs complete. Output written to rigidbody6dof_output.csv\n";
    return 0;
}
