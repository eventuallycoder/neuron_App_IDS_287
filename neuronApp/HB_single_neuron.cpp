/*
 Huber-Braun neuron model - termosensitive neuron
 Changing gsr parameter generates bursting

 compile with: g++ HB_single_neuron.cpp -o neuron
 run with: ./neuron

 output: voltage.txt

 visualize with: gnuplot 
    set xlabel "time (ms)"
    set ylabel "voltage (mV)"
    plot "voltage.txt" w l

 @Illinois State University -  Dr. Follmann & Dr. Rosa
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

using namespace std;

//------------------------------------------------------------------
//------------------------------------------------------------------
// Definitions and Classes
//------------------------------------------------------------------
//------------------------------------------------------------------

// Parameters
#define max_steps 200000
#define transient 200000 // 2
#define step_size 0.1

class Neurons
{
public:
    double X[4];

    Neurons();
    void set_parameters(double neuronParameters[]);
    void RK4(double Isyn);
    double equations(double Z[], int var, double Isyn);

private:
    // HB model parameters
    double gLeak, gNa, gK, gSD, gSR;
    double tauK, tauSD, tauSR;
    double sNa, sK, sSD;
    double V0Na, V0K, V0SD;
    double VLeak, VNa, VK, VSD, VSR;
    double ro, phi, theta, nu;
    double Cm;
};

//------------------------------------------------------------------
//------------------------------------------------------------------
// Main program
//------------------------------------------------------------------
//------------------------------------------------------------------

int main()
{
    FILE *fp = fopen("voltage.txt", "w");
    double gSR_var = 0.46;
    double time;
    double Iinj = 1.0;

    Neurons nrn; // Creates different neurons
    /*Parameters in the model only need to be changed here */
    double nrn_params[] = {
        0.1,     // gLeak
        1.5,     // gNa
        2.0,     // gK
        0.25,    // gSD
        gSR_var, // gSR
        2.0,     // tauK
        10.0,    // tauSD
        20.0,    // tauSR
        0.25,    // sNa
        0.25,    // sK
        0.09,    // sSD
        -25.0,   // V0Na
        -25.0,   // V0K
        -40.0,   // V0SD
        -60.0,   // Vleak
        50.0,    // VNa
        -90.0,   // VK
        50.0,    // VSD
        -90.0,   // VSR
        0.607,   // ro
        0.124,   // phi
        0.17,    // theta
        0.012,   // nu
        1.0      // Cm
    };

    // Initializes neuron parameters
    nrn.set_parameters(nrn_params);

    /** Transient run*/
    for (int n = 0; n < transient; n++)
    {
        nrn.RK4(Iinj);
    }

    /** --- Recording run-------
     *--------------------------
     */
    for (int n = 0; n < max_steps; n++)
    {
        time = n * step_size;
        nrn.RK4(Iinj); // solves neurons equation
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", time, nrn.X[0], nrn.X[1], nrn.X[2], nrn.X[3]);

    } // end time loop

    fclose(fp);

    cout << "Program Finished" << endl;
    return 0;
}

//------------------------------------------------------------------
//------------------------------------------------------------------
// Neurons Class
//------------------------------------------------------------------
//------------------------------------------------------------------

// Constructor
Neurons::Neurons()
{
    // Initialized states are based off Int1 Parameters
    gLeak = 0.1;
    gNa = 1.5;
    gK = 2.0;
    gSD = 0.25;
    gSR = 0.52;
    tauK = 2.0;
    tauSD = 10.0;
    tauSR = 20.0;
    sNa = 0.25;
    sK = 0.25;
    sSD = 0.09;
    V0Na = -25.0;
    V0K = -25.0;
    V0SD = -40.0;
    VLeak = -60.0;
    VNa = 50.0;
    VK = -90.0;
    VSD = 50.0;
    VSR = -90.0;
    ro = 0.607;
    phi = 0.124;
    theta = 0.17;
    nu = 0.012;
    Cm = 1.0;

    // Initial conditions
    X[0] = -60.0; // Soma
    X[1] = 0.1;   // K activation
    X[2] = 0.3;   // SD activation
    X[3] = 0.5;   // SR activation
}

// Parameter declaration
void Neurons::set_parameters(double neuronParameters[])
{
    gLeak = neuronParameters[0];
    gNa = neuronParameters[1];
    gK = neuronParameters[2];
    gSD = neuronParameters[3];
    gSR = neuronParameters[4];
    tauK = neuronParameters[5];
    tauSD = neuronParameters[6];
    tauSR = neuronParameters[7];
    sNa = neuronParameters[8];
    sK = neuronParameters[9];
    sSD = neuronParameters[10];
    V0Na = neuronParameters[11];
    V0K = neuronParameters[12];
    V0SD = neuronParameters[13];
    VLeak = neuronParameters[14];
    VNa = neuronParameters[15];
    VK = neuronParameters[16];
    VSD = neuronParameters[17];
    VSR = neuronParameters[18];
    ro = neuronParameters[19];
    phi = neuronParameters[20];
    theta = neuronParameters[21];
    nu = neuronParameters[22];
    Cm = neuronParameters[23];
}

// Standard RK4 for solving numerically differential equations
void Neurons::RK4(double Isyn)
{
    double k1[4], k2[4], k3[4], k4[4];
    double Xn[4];
    // Computes k1
    for (int i = 0; i < 4; i++)
        k1[i] = step_size * equations(X, i, Isyn);

    // Update temporary array
    for (int i = 0; i < 4; i++)
        Xn[i] = X[i] + k1[i] / 2.;

    // Computes k2
    for (int i = 0; i < 4; i++)
        k2[i] = step_size * equations(X, i, Isyn);

    // updates temporary array
    for (int i = 0; i < 4; i++)
        Xn[i] = X[i] + k2[i] / 2.;

    // Computes k3
    for (int i = 0; i < 4; i++)
        k3[i] = step_size * equations(X, i, Isyn);

    // updates temporary array
    for (int i = 0; i < 4; i++)
        Xn[i] = X[i] + k3[i];

    // computes k4
    for (int i = 0; i < 4; i++)
        k4[i] = step_size * equations(X, i, Isyn);

    // updates main array
    for (int i = 0; i < 4; i++)
        X[i] += (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]) / 6.0;
}
/** Neuron differetial equations */
double Neurons::equations(double Z[], int var, double Isyn)
{
    double val;
    double Ileak, INa, IK, ISD, ISR;
    double aKinf, aNainf, aSDinf;

    switch (var)
    {
    case 0:
        aNainf = 1.0 / (1.0 + exp(-sNa * (Z[0] - V0Na)));
        Ileak = gLeak * (Z[0] - VLeak);
        INa = ro * gNa * aNainf * (Z[0] - VNa);
        IK = ro * gK * Z[1] * (Z[0] - VK);
        ISD = ro * gSD * Z[2] * (Z[0] - VSD);
        ISR = ro * gSR * Z[3] * (Z[0] - VSR);

        val = -Ileak - INa - IK - ISD - ISR - Isyn;
        val /= Cm;
        break;
    case 1:
        aKinf = 1.0 / (1.0 + exp(-sK * (Z[0] - V0K)));
        val = phi * (aKinf - Z[1]) / tauK;
        break;

    case 2:
        aSDinf = 1.0 / (1.0 + exp(-sSD * (Z[0] - V0SD)));
        val = phi * (aSDinf - Z[2]) / tauSD;
        break;

    case 3:
        ISD = ro * gSD * Z[2] * (Z[0] - VSD);
        val = -phi * (nu * ISD + theta * Z[3]) / tauSR;
        break;
    }
    return val;
}
