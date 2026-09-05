/*
 Neuron Model - built up progressively alongside the Neuron Biology Study Plan

 Module 3: The Nernst Equation
 ------------------------------
 E_ion = (RT / zF) * ln(C_out / C_in)

 This gives the equilibrium potential for a single ion species - the
 voltage at which diffusive and electrical forces on that ion exactly
 cancel (net flux = 0).

 compile with: g++ neuron_model.cpp -o neuron_model
 run with: ./neuron_model
 */

#include <iostream>
#include <cmath>

using namespace std;

//------------------------------------------------------------------
// Physical constants
//------------------------------------------------------------------
#define R 8.314    // Gas constant, J/(mol*K)
#define F 96485.0  // Faraday's constant, C/mol
#define T_BODY 310.15 // Body temperature in Kelvin (37 C)

//------------------------------------------------------------------
// Nernst equation
//------------------------------------------------------------------
// c_out, c_in : concentrations (mM or any consistent unit)
// z           : valence (+1 for Na+/K+, -1 for Cl-, +2 for Ca2+)
// T           : temperature in Kelvin (defaults to body temp)
// returns     : E_ion in volts (multiply by 1000 for mV)
double nernst_potential(double c_out, double c_in, double z, double T = T_BODY)
{
    return (R * T) / (z * F) * log(c_out / c_in);
}

//------------------------------------------------------------------
// Module 4: The Goldman-Hodgkin-Katz (GHK) equation
//------------------------------------------------------------------
// Nernst (above) gives the equilibrium potential for a SINGLE ion in
// isolation. Real membranes have K+, Na+, and Cl- channels open at
// once, each pulling V_m toward its own E_ion. GHK computes the
// resulting compromise voltage, weighted by each ion's PERMEABILITY
// (not concentration alone) - an ion with more open channels pulls
// harder toward its own E_ion.
//
// V_m = (RT/F) * ln( (P_K*[K]_out + P_Na*[Na]_out + P_Cl*[Cl]_in)
//                   / (P_K*[K]_in  + P_Na*[Na]_in  + P_Cl*[Cl]_out) )
//
// Note Cl- is flipped (in on top, out on bottom) since it's the only
// negative ion here - same valence effect from Module 3, baked
// directly into the equation's structure this time.
//
// If P_Na == 0, this collapses exactly to the K+ Nernst equation -
// GHK is Nernst generalized to multiple simultaneous permeabilities,
// not a separate theory.
double ghk_voltage(double p_k, double p_na, double p_cl,
                    double k_out, double k_in,
                    double na_out, double na_in,
                    double cl_out, double cl_in,
                    double T = T_BODY)
{
    double numerator = p_k * k_out + p_na * na_out + p_cl * cl_in;
    double denominator = p_k * k_in + p_na * na_in + p_cl * cl_out;
    return (R * T) / F * log(numerator / denominator);
}

//------------------------------------------------------------------
// Main program
//------------------------------------------------------------------
int main()
{
    // Typical mammalian neuron concentrations (mM)
    struct Ion
    {
        const char *name;
        double c_out;
        double c_in;
        double z;
    };

    Ion ions[] = {
        {"K+", 5.0, 140.0, 1.0},
        {"Na+", 145.0, 15.0, 1.0},
        {"Cl-", 110.0, 10.0, -1.0},
    };

    cout << "Ion\tE_ion (mV)" << endl;
    for (const auto &ion : ions)
    {
        double e_volts = nernst_potential(ion.c_out, ion.c_in, ion.z);
        printf("%s\t%.1f\n", ion.name, e_volts * 1000.0);
    }

    // ------------------------------------------------------------
    // GHK: resting vs. spike-peak permeability regimes
    // ------------------------------------------------------------
    // Concentrations (mM), reused from the ions[] table above
    double k_out = 5.0, k_in = 140.0;
    double na_out = 145.0, na_in = 15.0;
    double cl_out = 110.0, cl_in = 10.0;

    // Resting regime: P_K dominant, small P_Na, small P_Cl
    double v_rest = ghk_voltage(1.0, 0.04, 0.45,
                                 k_out, k_in, na_out, na_in, cl_out, cl_in);
    cout << "\nGHK resting V_m (P_K dominant): "
         << v_rest * 1000.0 << " mV" << endl;

    // Spike-peak regime: P_Na now dominant (Na+ channels flung open)
    double v_peak = ghk_voltage(1.0, 20.0, 0.45,
                                 k_out, k_in, na_out, na_in, cl_out, cl_in);
    cout << "GHK spike-peak V_m (P_Na dominant): "
         << v_peak * 1000.0 << " mV" << endl;

    // ------------------------------------------------------------
    // Mapping to the Huber-Braun model's reversal potential constants
    // ------------------------------------------------------------
    // In HB_single_neuron.cpp, every current is written as g*(V - E),
    // exactly the driving-force form from Module 2. The E terms are
    // hardcoded there as fixed parameters instead of being computed
    // here, since the model assumes concentrations (and therefore
    // E_ion) stay constant on the timescale of a simulation run
    // (Module 5's non-equilibrium steady state assumption):
    //
    //   VNa  = +50 mV   <- should equal E_Na  (we computed +60.6 mV;
    //                       HB detunes it slightly to fit spike shape)
    //   VK   = -90 mV   <- should equal E_K   (matches our -89.1 mV
    //                       almost exactly)
    //   VSD  = +50 mV   <- biologically E_Ca (slow depolarizing current
    //                       is Ca2+ based), but HB reuses a Na-like
    //                       value since the gating dynamics matter more
    //                       than the exact reversal here
    //   VSR  = -90 mV   <- reuses E_K, since SR is a slow Ca2+-activated
    //                       K+ current (same ion family as IK)
    //
    // None of these are derived dynamically in the HB model - they are
    // constants baked into nrn_params[]. A more advanced extension
    // could replace VK with a live variable updated by its own
    // differential equation (e.g. tracking extracellular K+ buildup),
    // which would make E_K dynamic instead of fixed.

    return 0;
}
