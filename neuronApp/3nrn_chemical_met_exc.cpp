// This code is intended to simulate a network of neurons
// connected in a user-defined way via chemical synapse
// Illinois State Universtiy 2014, Dr. Rosa's Lab

//setting up standard library imports 
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

//------------------------------------------------------------------
//------------------------------------------------------------------
// Definitions and Classes
// understand code and where things are in the code
// play with nrnnum start at 2 and then 3 
//------------------------------------------------------------------
//------------------------------------------------------------------

// Parameters

//how long the simulation runs for
#define max_steps 500000

//A period of tie to run the simulation without recording data
#define transient 0 // 700000//2

//the time step  (dt) for the integreator
#define step_size 0.01

//The number of neurons in the model
#define nrnNumber 3

//the state variable number for each neuron (Voltage, ion channel states, synaptic concentrations, etc.)
#define nVar 8

// creates a class to handle numerical methods (RK4)
class Neurons
{
public:

/*this holds the current state of the neuron. X[0] is voltage, 
  and X[1] through X[7] represents internal variables like potassium 
  activations or neurotransmitter levels*/
  double X[nVar];

  Neurons();
  //customizes the neuron
  void set_parameters(int number);
  //the runge kutta 4 solver. This advances the simulation one step forward in time
  void RK4(int iterator, int number, double V_ps[], bool SC[], double V_elec[]);

private:
//the list of differential equations describing biological physics.
  double equations(double Z[], int var);

  // HB model parameters

  /*conductances. gK and gNa are the main sources of action potential
  whereas gSD and gSR are what causes bursting. gLeak is the leakage
  channel trying to keep the neuron at resting potential */
  double gLeak, gNa, gK, gSD, gSR;

  /*represent value speeds. A large tau means hte channel opens or closes
  very slowly, and a small tau means it reacts instantly to voltage change.
  Na's channel is assumed to open up instantly which is why it doesn't have a 
  time variable*/
  double tauK, tauSD, tauSR;

  /*These determine how sensitive the ion channels are to voltage changes. 
  A high value means the channel snaps open effectively all at once when a 
  threshold is hit; a low value means it opens gradually.*/
  double sNa, sK, sSD;

  /* The voltage threshold at which 50% of these specific ion channels are open*/
  double V0Na, V0K, V0SD;
  double VLeak, VNa, VK, VSD, VSR;
  double ro, phi, theta, nu;
  double Cm;

  double V_preSynaptic[2], V_spike[2], V_gapJunction[2];
  double C_release[2], C_spike[2];
  double t_accum[2], t_elim[2], t_act[2], t_inact[2];
  double s_syn[2], p_mt[2], r_mt[2], C_mt[2];
  double g_mt[2], E_mt[2];
  double g_elec[2];

  bool spikeCondition[2];
  bool ChemConnection[nrnNumber - 1], ElecConnection[nrnNumber - 1];
};

//------------------------------------------------------------------
//------------------------------------------------------------------
// Main program
//------------------------------------------------------------------
//------------------------------------------------------------------

double g_chemical = 0.0;

int main()
{
  FILE *fp;
  FILE *fp1;
  FILE *fp2;
  FILE *fp_ISI0 = fopen("ISI_nrn0_case1_MT_CE_dnmwtest1.txt", "w");
  FILE *fp_ISI1 = fopen("ISI_nrn1_case1_MT_CE_dnmwtest1.txt", "w");
  // FILE *fp_ISI2 = fopen("ISI_nrn2.txt","w");

  double Crelease1;
  double p_meta1;
  double Crelease2;
  double p_meta2;
  double C_spike1 = 1.0;
  double s_syn1 = 1.0;
  double V_spike1 = -20.0;
  double r_mt1 = 1.0;
  double C_mt1 = 1.0;
  double p_mt1 = 0.2;

  for (g_chemical = 0.; g_chemical <= 0.0001; g_chemical += 0.001)
  {

    cout << g_chemical << endl;

    double V_currentVal[nrnNumber], V_previousVal[nrnNumber];
    double V_electric[nrnNumber - 1];
    double thresh = -20.0;
    double spikeTime[nrnNumber];

    int spikeNumber[nrnNumber];

    bool spike[nrnNumber];

    //turn this off (false) when running with larger gChemical range.
    bool V_print = true;

    if (V_print)
    {
      fp = fopen("NeuronVoltages.txt", "w");
      fp1 = fopen("data1.txt", "w");
      fp2 = fopen("data2.txt", "w");
    }

    // Creates different neurons
    Neurons nrn[nrnNumber];

    // Initializes neuron parameters and variables
    for (int num = 0; num < nrnNumber; num++)
    {
      nrn[num].set_parameters(num);
      V_currentVal[num] = 0.0;
      V_previousVal[num] = 0.0;
      spikeTime[num] = 0.0;
      spike[num] = false;
      V_electric[num] = 0.0;
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
    // if(g_chemical>0.032){max_steps_refined /= 10;}

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

        Crelease1 = C_spike1 / (1 + exp(-s_syn1 * (V_currentVal[0] - V_spike1)));

        p_meta1 = r_mt1 * nrn[0].X[4] / (C_mt1 + nrn[0].X[4]);

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
              fprintf(fp_ISI0, "%5.4lf\t%7.2lf\n", g_chemical, n * step_size - spikeTime[num]);
              fflush(fp_ISI0);
            }
            break;

          case 1:
            if (spikeTime[num] > 0.0)
            {
              fprintf(fp_ISI1, "%5.4lf\t%7.2lf\n", g_chemical, n * step_size - spikeTime[num]);
              fflush(fp_ISI1);
            }
            break;

          case 2:
            if (spikeTime[num] > 0.0)
            {
              // fprintf(fp_ISI2,"%5.3lf\t%7.3lf\n",g_chemical,n*step_size-spikeTime[num]);
              // fflush(fp_ISI2);
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

          // Print neuron 0 data
          //  fprintf(fp1,"%7.3lf\t%7.3lf\t%7.3lf\n",n*step_size,nrn[0].X[0],nrn[0].X[1]);
          fprintf(fp1, "%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\n", n * step_size, nrn[0].X[0], nrn[0].X[1], nrn[0].X[2], nrn[0].X[3], nrn[0].X[4], nrn[0].X[5], nrn[0].X[6], nrn[0].X[7], Crelease1, p_meta1);
          // print neuron 1 data
          fprintf(fp2, "%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\n", n * step_size, nrn[1].X[0], nrn[1].X[1], nrn[1].X[2], nrn[1].X[3], nrn[1].X[4], nrn[1].X[5], nrn[1].X[6], nrn[1].X[7]);
        }
        // Print neuron 0 data
        //  fprintf(fp1,"%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\n",n*step_size,nrn[0].X[0],nrn[0].X[1],nrn[0].X[2],nrn[0].X[3],nrn[0].X[4],nrn[0].X[5],nrn[0].X[6],nrn[0].X[7]);

        // print neuron 1 data
        // fprintf(fp2,"%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\n",n*step_size,nrn[1].X[0],nrn[1].X[1],nrn[1].X[2],nrn[1].X[3],nrn[1].X[4],nrn[1].X[5],nrn[1].X[6],nrn[1].X[7]);
      }

      for (int num = 0; num < nrnNumber; num++)
      {
        V_previousVal[num] = nrn[num].X[0];
      }
    }

    if (V_print)
    {
      fclose(fp);
      fclose(fp1);
      fclose(fp2);
    }
  }
  fclose(fp_ISI0);
  fclose(fp_ISI1);

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
  X[4] = 1.;    // 4.0; //cleft concentration
  X[5] = .8;    // 0.6; //synaptic current activation
  X[6] = 0.4;   // cleft concentration for second synapse
  X[7] = 0.6;   // synaptic current activation for second synapse
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
    gSR = 0.26;
    g_mt[0] = -g_chemical; // Chemical conductance of synapse from nrn1
    g_mt[1] = 0.02;        // Chemical conducatce of synapse from nrn2
    g_elec[0] = 0.075;     // Electrical conductace from nrn1 to nrn0
    g_elec[1] = 0.05;      // Electrical conductance form nrn2 to nrn0

    ChemConnection[0] = true;  // Chemical Connection from nrn0 to nrn1
    ChemConnection[1] = false; // connection from nrn0 to nrn2

    ElecConnection[0] = false; // Electircal Connection from nrn0 to nrn1
    ElecConnection[1] = false; // Electircal Connection from nrn0 to nrn2
    break;
  case 1:
    gSD = 0.25;
    gSR = 0.38;
    g_mt[0] = -g_chemical; // Chemical conductance of synapse from nrn0
    g_mt[1] = 0.02;        // Chemical conductance of synapse from nrn2
    g_elec[0] = 0.075;     // Electrical conductace from nrn0
    g_elec[1] = 0.05;      // Electrical conductance form nrn2

    ChemConnection[0] = true;  // Connection from nrn1 to nrn0
    ChemConnection[1] = false; // connection from nrn1 to nrn2

    ElecConnection[0] = false; // Electircal Connection from nrn0 to nrn1
    ElecConnection[1] = false; // Electircal Connection from nrn1 to nrn1
    break;
  case 2:
    gSD = 0.25;
    gSR = 0.46;
    g_mt[0] = 0.02;   // Chemical conductance of synapse from nrn1
    g_mt[1] = 0.02;   // Chemical conductance of synapse from nrn2
    g_elec[0] = 0.75; // Electrical conductace from nrn0
    g_elec[1] = 0.05; // Electrical conductace from nrn1

    ChemConnection[0] = false; // connection from nrn2 to nrn0
    ChemConnection[1] = false; // connection from nrn2 to nrn1

    ElecConnection[0] = false; // Electircal Connection from nrn0 to nrn2
    ElecConnection[1] = false; // Electircal Connection from nrn1 to nrn2
    break;
  }

  for (int con = 0; con < 2; con++)
  {
    C_spike[con] = 1.0;
    s_syn[con] = -3.; // 1.0;
    V_spike[con] = -20.0;
    t_accum[con] = 10.; // 10.0;
    t_elim[con] = 20.;  // 20.0;
    r_mt[con] = 1.;     // 1.0;
    C_mt[con] = 1.0;    // 1.0;
    p_mt[con] = 1.0;
    t_act[con] = 1.;    // 10.0;
    t_inact[con] = 15.; // 15.0;
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
      ChemConnection[0] = true;
      ChemConnection[1] = false;
    }
    break;

  case 1:
    // if(iterator==int((max_steps+transient)/2))
    if (iterator == 0)
    {
      ChemConnection[0] = true; //   true
      ChemConnection[1] = false;
    }
    break;

  case 2:
    // if(iterator==int((max_steps+transient)/2))
    if (iterator == 0)
    {
      ChemConnection[0] = false;
      ChemConnection[1] = false;
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
      Isynaptic[0] = g_mt[0] * Z[5] * (Z[0] - E_mt[0]);
    }
    else
    {
      Isynaptic[0] = 0.0;
    }

    if (spikeCondition[1] && ChemConnection[1])
    {
      Isynaptic[1] = g_mt[1] * Z[5] * (Z[0] - E_mt[1]);
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

    val = -Ileak - INa - IK - ISD - ISR + Isynaptic[0] + Isynaptic[1] - I_gapJunction[0] - I_gapJunction[1] - 1.0;
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
  // Ionic/transmiter concentration of synapse 1
  case 4:
    C_release[0] = C_spike[0] / (1. + exp(-s_syn[0] * (V_preSynaptic[0] - V_spike[0])));
    val = C_release[0] / t_accum[0] - Z[4] / t_elim[0];
    break;
  // synaptic/receptor activation of synapse 1
  case 5:
    p_mt[0] = r_mt[0] * Z[4] / (C_mt[0] + Z[4]);
    val = p_mt[0] / t_act[0] - Z[5] / t_inact[0];
    break;
  // Ionic/transmiter concentration of synapse 1
  case 6:
    C_release[1] = C_spike[1] / (1. + exp(-s_syn[1] * (V_preSynaptic[1] - V_spike[1])));
    val = C_release[1] / t_accum[1] - Z[4] / t_elim[1];
    break;
  // synaptic activation of synapse 1
  case 7:
    p_mt[1] = r_mt[1] * Z[4] / (C_mt[1] + Z[4]);
    val = p_mt[1] / t_act[1] - Z[5] / t_inact[1];
    break;
    //  case 8:
    //    val = C_spike[0]/(1. + exp(-s_syn[0]*(V_preSynaptic[0]-V_spike[0])));
    //  break;
  }
  return val;
}
