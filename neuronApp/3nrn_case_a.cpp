/*
Written by: Quinton Skilling
Under: Dr. Epaminondas Rosa
Institute: Illinois State University
Date: July 2013

This code is intended to simulate a two-neuron coupled system. Each
HB neuron has different dynamics relating to the slow repolarization and slow
depolarization conductances, gSR and gSD, respectively, connected via chemical
synapse only OR chemical synapse AND electrical gap-junction.

The goal is to investigate the effect of inhibitory synaptic connections on the
membrane potential of neurons. Briefly, presynaptic inhibitory connections cause repolarizations
of the postsynaptic neuron which can slow down the interspike interval (ISI) considerably, even
to the point of "shutting down," i.e. no generation of action potentials, the postsynaptic neuron.
By creating a bifurcation diagram of the ISI vs the inhibitory synaptic conduction,
the effect of inhibitory connections can be visualized.

In case (a), one of the neurons is set in the fast tonic regime while the other is set
in the regular bursting regime (see Figure 2). When only the inhibitory synapse is active and sufficiently
strong (both neurons experience the same weight of inhibition) the spiking of the fast tonic
neuron eventually shuts down the bursting neuron. Similarly, when the electrical gap-junction is
turned on, causing the two neurons to synchronize, along with the inhibitory synapse, the fast neuron again shuts
down the bursting neuron, though the synaptic weight required is less than in the case where only the synapse is active.

For all cases dealing with the electrical gap-junction, the corresponding conductance was fixed at 0.06 mS/cm^2.
For a larger frequency difference in neurons, however, it should be noted that a stronger electrical
conductance may be needed in order to synchronize the neurons.

For more information about specific algorithms used in this code, please consult the comments listed with the
corresponding algorithms below.
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>

using namespace std;

//------------------------------------------------------------------
//------------------------------------------------------------------
// Definitions and Classes
//------------------------------------------------------------------
//------------------------------------------------------------------

// Parameters
#define max_steps 1000000
#define transient 700000 // 2
#define step_size 0.05
#define nrnNumber 3

// creates a class to handle numerical methods (RK4)
class Neurons
{
public:
    double X[8];

    Neurons();
    void set_parameters(int number);
    void set_tempParameters(double T);
    void RK4(int iterator, int number, double V_ps[], bool SC[], double V_elec[]);

private:
    double equations(double Z[], int var);

    // HB model parameters
    double gLeak, gNa, gK, gSD, gSR;
    double tauK, tauSD, tauSR;
    double sNa, sK, sSD;
    double V0Na, V0K, V0SD;
    double VLeak, VNa, VK, VSD, VSR;
    double ro, phi, theta, nu;
    double Cm;

    double V_preSynaptic[2], V_spike[2], V_gapJunction[2];
    double C_release[2], C_spike[2];
    double t_accum[2], t_elim[2], t_act[2], t_inact[2];
    double s_syn[2], p_mt[2], r_mt[2], C_mt[2];
    double g_ch[2], E_mt[2];
    double g_elec[2], g_elec_base[2];

    bool spikeCondition[2];
    bool ChemConnection[nrnNumber - 1], ElecConnection[nrnNumber - 1];
};

//------------------------------------------------------------------
//------------------------------------------------------------------
// Main program
//------------------------------------------------------------------
//------------------------------------------------------------------

double g_chemical = 0.0;  // initial value for g_chemical
double g_electric = 1.50; // static value for g_electric
bool is_chemical = true;
bool is_electric = false;

// g_elec_base is eta
int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "Usage: " << argv[0] << " <neuron temp> " << endl;
        return 1;
    }
    double neuronTemp = atof(argv[1]);

    stringstream ss;
    ss << fixed << setprecision(0) << neuronTemp;
    string tempString = ss.str();

    FILE *fp;
    FILE *fp_ISI0;
    string outputName = "isi1_" + tempString + ".txt";
    fp_ISI0 = fopen(outputName.c_str(), "w");

    FILE *fp_ISI1;
    string outputName2 = "isi2_" + tempString + ".txt";
    fp_ISI1 = fopen(outputName2.c_str(), "w");

    // Added file printer for neuron 2
    FILE *fp_ISI2;
    string outputName3 = "isi3_" + tempString + ".txt";
    fp_ISI2 = fopen(outputName3.c_str(), "w");

    // opening the synchronized file
    FILE *SynchronizedFile;
    SynchronizedFile = fopen("synchronized.txt", "w");

    cout << "Running with Temp:" << neuronTemp << endl;

    for (g_chemical = 0.0; g_chemical < 0.05; g_chemical += 0.01)
    {

        // keeps track of the number of spikes for each neuron;
        int num_spikes_nrn0 = 0;
        int num_spikes_nrn1 = 0;
        int num_spikes_nrn2 = 0;

        cout << "gch Value: " << g_chemical << endl;
        double V_currentVal[nrnNumber], V_previousVal[nrnNumber];
        double V_electric[nrnNumber - 1];
        double thresh = -20.0;
        double spikeTime[nrnNumber];

        int spikeNumber[nrnNumber];

        bool spike[nrnNumber];
        bool V_print = false;

        if (V_print)
        {
            fp = fopen("data.txt", "w");
        }

        // Creates different neurons
        Neurons nrn[nrnNumber];

        // Initializes neuron parameters and variables
        for (int num = 0; num < nrnNumber; num++)
        {
            nrn[num].set_parameters(num);
            nrn[num].set_tempParameters(neuronTemp);
            V_currentVal[num] = 0.0;
            V_previousVal[num] = 0.0;
            spikeTime[num] = 0.0;
            spike[num] = false;
            if (num < nrnNumber - 1)
            {
                V_electric[num] = 0.0;
            }
        }

        // Initial transient
        for (int n = 0; n < transient; n++)
        {
            for (int num = 0; num < nrnNumber; num++)
            {
                // If it is the first neuron, the potetial difference is calcualted with respect
                // to its right-most nearest neighbor and the last neuron
                if (num == 0)
                {
                    // Calculates the potential difference
                    V_electric[0] = nrn[num].X[0] - nrn[num + 1].X[0];
                    V_electric[1] = nrn[num].X[0] - nrn[nrnNumber - 1].X[0];
                }
                // If it is the last neuron, the potetial difference is calcualted with respect
                // to its left-most nearest neighbor and the first neuron
                else if (num == nrnNumber - 1)
                {
                    V_electric[0] = nrn[num].X[0] - nrn[num - 1].X[0];
                    V_electric[1] = nrn[num].X[0] - nrn[0].X[0];
                }
                // if it is a middle neuron, the potetial difference is calcualted with respect
                // to its nearest neighbors (1 on each side)
                else
                {
                    V_electric[0] = nrn[num].X[0] - nrn[num - 1].X[0];
                    V_electric[1] = nrn[num].X[0] - nrn[num + 1].X[0];
                }

                // calls the RK4
                nrn[num].RK4(n, num, V_currentVal, spike, V_electric);

                // checks if the neuron is spiking (depolarized)
                if (n != 0 && nrn[num].X[0] > thresh)
                {
                    spike[num] = true;
                    V_currentVal[num] = nrn[num].X[0];
                }
                else if (n != 0 && nrn[num].X[0] < thresh)
                {
                    spike[num] = false;
                }
            }

            // updates the previous value and resets the spiking condition boolean
            for (int num = 0; num < nrnNumber; num++)
            {
                V_previousVal[num] = nrn[num].X[0];
            }
        }

        int max_steps_refined = max_steps;
        // if(g_chemical>0.02){max_steps_refined /= 10;}

        // Recordable data
        for (int n = 0; n < max_steps_refined; n++)
        {
            // calls the RK4
            for (int num = 0; num < nrnNumber; num++)
            {
                // If it is the first neuron, the potetial difference is calcualted with respect
                // to its right-most nearest neighbor and the last neuron
                if (num == 0)
                {
                    // Calculates the potential difference
                    V_electric[0] = nrn[num].X[0] - nrn[num + 1].X[0];
                    V_electric[1] = nrn[num].X[0] - nrn[nrnNumber - 1].X[0];
                }
                // If it is the last neuron, the potetial difference is calcualted with respect
                // to its left-most nearest neighbor and the first neuron
                else if (num == nrnNumber - 1)
                {
                    V_electric[0] = nrn[num].X[0] - nrn[num - 1].X[0];
                    V_electric[1] = nrn[num].X[0] - nrn[0].X[0];
                }
                // if it is a middle neuron, the potetial difference is calcualted with respect
                // to its nearest neighbors (1 on each side)
                else
                {
                    V_electric[0] = nrn[num].X[0] - nrn[num - 1].X[0];
                    V_electric[1] = nrn[num].X[0] - nrn[num + 1].X[0];
                }

                nrn[num].RK4(n, num, V_currentVal, spike, V_electric);

                // checks for a spike
                if (n != 0 && nrn[num].X[0] > thresh)
                {
                    spike[num] = true;
                    V_currentVal[num] = nrn[num].X[0];
                }
                else if (n != 0 && nrn[num].X[0] < thresh)
                {
                    spike[num] = false;
                }

                if (nrn[num].X[0] > thresh && V_previousVal[num] < thresh)
                {
                    switch (num)
                    {
                    case 0:
                        if (spikeTime[num] > 0.0)
                        {

                            num_spikes_nrn0++; // increment the number of spikes by 1
                            fprintf(fp_ISI0, "%5.3lf\t%7.3lf\n", g_chemical, n * step_size - spikeTime[num]);
                            fflush(fp_ISI0);
                        }
                        break;

                    case 1:
                        if (spikeTime[num] > 0.0)
                        {
                            num_spikes_nrn1++; // increment the number of spikes by 1
                            fprintf(fp_ISI1, "%5.3lf\t%7.3lf\n", g_chemical, n * step_size - spikeTime[num]);
                            fflush(fp_ISI1);
                        }
                        break;

                    case 2:
                        if (spikeTime[num] > 0.0)
                        {
                            num_spikes_nrn2++; // increment the number of spikes by 1
                            fprintf(fp_ISI2, "%5.3lf\t%7.3lf\n", g_chemical, n * step_size - spikeTime[num]);
                            fflush(fp_ISI2);
                        }
                        break;
                    }
                    spikeTime[num] = n * step_size;
                }
            }

            if (V_print)
            {
                if (n % 1 == 0)
                {
                    // Prints the time, potential, K-activation, SD-Activation, and SR-activation
                    fprintf(fp, "%8.2lf", n * step_size);

                    for (int num = 0; num < nrnNumber; num++)
                    {
                        fprintf(fp, "\t%7.3lf", nrn[num].X[0]);
                    }

                    fprintf(fp, "\n");
                }
            }

            for (int num = 0; num < nrnNumber; num++)
            {
                V_previousVal[num] = nrn[num].X[0];
            }
        }

        std::cout << num_spikes_nrn0 << " " << num_spikes_nrn1 << " " << " " << num_spikes_nrn2 << std::endl;

        //check if neurons are synced or not
        if(num_spikes_nrn0 == num_spikes_nrn1 == num_spikes_nrn2)
        {
            fprintf(SynchronizedFile, "Neurons are synchronized for g_chem = %ld", g_chemical);
        }

        if (V_print)
        {
            fclose(fp);
        }
    }
    fclose(fp_ISI0);
    fclose(fp_ISI1);
    fclose(fp_ISI2);

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

    for (int num2 = 0; num2 < 2; num2++)
    {
        spikeCondition[num2] = false;
        V_preSynaptic[num2] = false;
    }

    X[0] = -60.0; // Soma
    X[1] = 0.1;   // K activation
    X[2] = 0.3;   // SD activation
    X[3] = 0.5;   // SR activation
    X[4] = 4.0;   // cleft concentration
    X[5] = 0.6;   // synaptic current activation
    X[6] = 0.4;   // cleft concentration for second synapse
    X[7] = 0.6;   // synaptic current activation for second synapse
}

void Neurons::set_tempParameters(double T)
{
    phi = pow(3.0, (T - 60.0) / 10.0);
    ro = pow(1.3, (T - 60.0) / 10.0);

    double tempFactor = pow(4.0, (T - 60.0) / 10.0);
    g_elec[0] = g_elec_base[0] * tempFactor;
    g_elec[1] = g_elec_base[1] * tempFactor;
}

// Parameter decleration
void Neurons::set_parameters(int number)
{
    gLeak = 0.1;
    gNa = 1.5;
    gK = 2.0;
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

    // These are the parameters you want to change
    // If you have more than two neurons, add a case statement with parameters (as below)
    switch (number)
    {
    case 0:
        gSD = 0.25;
        gSR = 0.20;
        g_ch[0] = -g_chemical;       // Chemical conductance of synapse from nrn1
        g_ch[1] = 0.02;              // Chemical conducatce of synapse from nrn2
        g_elec_base[0] = g_electric; // Electrical conductace from nrn1 to nrn0
        g_elec_base[1] = g_electric; // Electrical conductance form nrn2 to nrn0

        ChemConnection[0] = is_chemical; // Chemical Connection from nrn0 to nrn1
        ChemConnection[1] = is_chemical; // connection from nrn0 to nrn2

        ElecConnection[0] = is_electric; // Electircal Connection from nrn0 to nrn1
        ElecConnection[1] = is_electric; // Electircal Connection from nrn0 to nrn2
        break;
    case 1:
        gSD = 0.25;
        gSR = 0.38;
        g_ch[0] = -g_chemical;       // Chemical conductance of synapse from nrn0
        g_ch[1] = 0.02;              // Chemical conductance of synapse from nrn2
        g_elec_base[0] = g_electric; // Electrical conductace from nrn0
        g_elec_base[1] = g_electric; // Electrical conductance form nrn2

        ChemConnection[0] = is_chemical; // Connection from nrn1 to nrn0
        ChemConnection[1] = is_chemical; // connection from nrn1 to nrn2

        ElecConnection[0] = is_electric; // Electircal Connection from nrn0 to nrn1
        ElecConnection[1] = is_electric; // Electircal Connection from nrn1 to nrn1
        break;
    case 2:
        gSD = 0.25;
        gSR = 0.46;
        g_ch[0] = -g_chemical;       // Chemical conductance of synapse from nrn1
        g_ch[1] = 0.02;              // Chemical conductance of synapse from nrn2
        g_elec_base[0] = g_electric; // Electrical conductace from nrn0
        g_elec_base[1] = g_electric; // Electrical conductace from nrn1

        ChemConnection[0] = is_chemical; // connection from nrn2 to nrn0
        ChemConnection[1] = is_chemical; // connection from nrn2 to nrn1

        ElecConnection[0] = is_electric; // Electircal Connection from nrn0 to nrn2
        ElecConnection[1] = is_electric; // Electircal Connection from nrn1 to nrn2
        break;
    }

    for (int con = 0; con < 2; con++)
    {
        C_spike[con] = 1.0;
        s_syn[con] = 1.0;
        V_spike[con] = -20.0;
        t_accum[con] = 10.0;
        t_elim[con] = 20.0;
        r_mt[con] = 1.0;
        C_mt[con] = 1.0;
        p_mt[con] = 0.2;
        t_act[con] = 10.0;
        t_inact[con] = 15.0;
        E_mt[con] = 50.0;
    }
}

// Standard RK4 for 8 differential equations
void Neurons::RK4(int iterator, int number, double V_ps[], bool SC[], double V_elec[])
{
    double k1[8], k2[8], k3[8], k4[8];
    double Xn[8];

    int arrPos = 0;

    switch (number)
    {
    case 0:
        // gives the option to turn on/off chemical connections
        // if(iterator==int((max_steps+transient)/2))
        if (iterator == 0)
        {
            ChemConnection[0] = is_chemical;
            ChemConnection[1] = is_chemical;
        }
        break;

    case 1:
        // if(iterator==int((max_steps+transient)/2))
        if (iterator == 0)
        {
            ChemConnection[0] = is_chemical;
            ChemConnection[1] = is_chemical;
        }
        break;

    case 2:
        // if(iterator==int((max_steps+transient)/2))
        if (iterator == 0)
        {
            ChemConnection[0] = is_chemical;
            ChemConnection[1] = is_chemical;
        }
        break;
    }

    // save presynaptic voltage to be used in synapse. Checks if the connecting neurons are spiking
    for (int num2 = 0; num2 < nrnNumber; num2++)
    {
        if (number == num2)
        {
            continue;
        }

        V_preSynaptic[arrPos] = V_ps[num2];
        V_gapJunction[arrPos] = V_elec[arrPos];
        spikeCondition[arrPos] = SC[num2];

        arrPos++;
    }

    // Computes k1
    for (int i = 0; i < 8; i++)
    {
        // skips the synaptic equation if the connecting neuron isn't spiking
        if (i > 3 && i < 6 && !spikeCondition[0])
        {
            continue;
        }
        if (i > 5 && !spikeCondition[1])
        {
            continue;
        }

        k1[i] = step_size * equations(X, i);
    }

    // Update temporary array
    for (int i = 0; i < 8; i++)
    {
        if (i > 3 && i < 6 && !spikeCondition[0])
        {
            continue;
        }
        if (i > 5 && !spikeCondition[1])
        {
            continue;
        }

        Xn[i] = X[i] + k1[i] / 2.;
    }

    // Computes k2 with temp array
    for (int i = 0; i < 8; i++)
    {
        if (i > 3 && i < 6 && !spikeCondition[0])
        {
            continue;
        }
        if (i > 5 && !spikeCondition[1])
        {
            continue;
        }

        k2[i] = step_size * equations(Xn, i);
    }

    // updates temporary array
    for (int i = 0; i < 8; i++)
    {
        if (i > 3 && i < 6 && !spikeCondition[0])
        {
            continue;
        }
        if (i > 5 && !spikeCondition[1])
        {
            continue;
        }

        Xn[i] = X[i] + k2[i] / 2.;
    }

    // Computes k3 with temp array
    for (int i = 0; i < 8; i++)
    {
        if (i > 3 && i < 6 && !spikeCondition[0])
        {
            continue;
        }
        if (i > 5 && !spikeCondition[1])
        {
            continue;
        }

        k3[i] = step_size * equations(Xn, i);
    }

    // updates temporary array
    for (int i = 0; i < 8; i++)
    {
        if (i > 3 && i < 6 && !spikeCondition[0])
        {
            continue;
        }
        if (i > 5 && !spikeCondition[1])
        {
            continue;
        }

        Xn[i] = X[i] + k3[i];
    }

    // computes k4 with temp array
    for (int i = 0; i < 8; i++)
    {
        if (i > 3 && i < 6 && !spikeCondition[0])
        {
            continue;
        }
        if (i > 5 && !spikeCondition[1])
        {
            continue;
        }

        k4[i] = step_size * equations(Xn, i);
    }

    // updates main array
    for (int i = 0; i < 8; i++)
    {
        if (i > 3 && i < 6 && !spikeCondition[0])
        {
            continue;
        }
        if (i > 5 && !spikeCondition[1])
        {
            continue;
        }

        X[i] += (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]) / 6.0;
    }
}

// Huber-Braun differential equations
double Neurons::equations(double Z[], int var)
{
    double val;

    double Ileak, INa, IK, ISD, ISR;
    double aKinf, aNainf, aSDinf;
    double Isynaptic[2], I_gapJunction[2];

    double alpha = 1.0;

    switch (var)
    {
    // Soma -- "Houses" activation currents
    case 0:
        aNainf = 1.0 / (1.0 + exp(-sNa * (Z[0] - V0Na)));
        Ileak = gLeak * (Z[0] - VLeak);
        INa = ro * gNa * aNainf * (Z[0] - VNa);
        IK = ro * gK * Z[1] * (Z[0] - VK);
        ISD = ro * gSD * Z[2] * (Z[0] - VSD);
        ISR = ro * gSR * Z[3] * (Z[0] - VSR);

        // if the connecting neuron isn't spiking, the synaptic current is zero
        if (spikeCondition[0] && ChemConnection[0])
        {
            Isynaptic[0] = g_ch[0] * Z[5] * (Z[0] - E_mt[0]);
        }
        else
        {
            Isynaptic[0] = 0.0;
        }

        if (spikeCondition[1] && ChemConnection[1])
        {
            Isynaptic[1] = g_ch[1] * Z[5] * (Z[0] - E_mt[1]);
        }
        else
        {
            Isynaptic[1] = 0.0;
        }

        // If neurons are not electrically coupled, the gap junction current is zero
        if (ElecConnection[0])
        {
            I_gapJunction[0] = g_elec[0] * V_gapJunction[0];
        }
        else
        {
            I_gapJunction[0] = 0.0;
        }

        if (ElecConnection[1])
        {
            I_gapJunction[1] = g_elec[1] * V_gapJunction[1];
        }
        else
        {
            I_gapJunction[1] = 0.0;
        }

        val = -Ileak - INa - IK - ISD - ISR - I_gapJunction[0] - I_gapJunction[1] - Isynaptic[0] - Isynaptic[1] - 1.0;
        val /= Cm;
        break;
    // Activation for Potasium current
    case 1:
        aKinf = 1.0 / (1.0 + exp(-sK * (Z[0] - V0K)));
        val = phi * (aKinf - Z[1]) / tauK;
        break;
    // Activation for slow-depolarization (Sodium)
    case 2:
        aSDinf = 1.0 / (1.0 + exp(-sSD * (Z[0] - V0SD)));
        val = phi * (aSDinf - Z[2]) / tauSD;
        break;
    // Activation for slow-hyperpolarization (Potasium)
    case 3:
        ISD = ro * gSD * Z[2] * (Z[0] - VSD);
        val = -phi * (nu * ISD + theta * Z[3]) / tauSR;
        break;
    // Ionic concentration of synapse 1
    case 4:
        C_release[0] = C_spike[0] / (1. + exp(-s_syn[0] * (V_preSynaptic[0] - V_spike[0])));
        val = C_release[0] / t_accum[0] - Z[4] / t_elim[0];
        break;
    // synaptic activation of synapse 1
    case 5:
        p_mt[0] = r_mt[0] * Z[4] / (C_mt[0] + Z[4]);
        val = p_mt[0] / t_act[0] - Z[5] / t_inact[0];
        break;
    // Ionic concentration of synapse 1
    case 6:
        C_release[1] = C_spike[1] / (1. + exp(-s_syn[1] * (V_preSynaptic[1] - V_spike[1])));
        val = C_release[1] / t_accum[1] - Z[6] / t_elim[1];
        break;
    // synaptic activation of synapse 1
    case 7:
        p_mt[1] = r_mt[1] * Z[6] / (C_mt[1] + Z[6]);
        val = p_mt[1] / t_act[1] - Z[7] / t_inact[1];
        break;
    }
    return val;
}