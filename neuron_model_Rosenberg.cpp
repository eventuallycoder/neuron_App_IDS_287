// Personal implementation of the single neuron code -- Aaron Rosenberg

//compile with: g++ neuron_model_Rosenberg.cpp -o model_R
#include <math.h> //math.h library used to perform ln()
#include <iostream>
#include <string>
class Neuron
{
    public:
        static constexpr double gas_constant_R = 8.314; //represents the gas constant measured in (J/mol*K)
        static constexpr int faradays_constant_F = 96485; // represents Faraday's constant measured in C/mol
        static constexpr double body_temperature = 310.15; //represents the human body temperature, measured in Kelvin
        
        //nernst potential: computes a single ion's equillibrium potential
        double nernst_potential(double conductance_out, double conductance_in, double valence_z, double T = body_temperature)
        {
            return ((gas_constant_R * T) / (valence_z * faradays_constant_F)) * log(conductance_out/conductance_in); // E_ion in volts (not mV)
        }

        //GHK function, the summation of the nernst potentials from each ion
        double ghk_voltage(double p_k, double p_na, double p_cl, double k_out, double k_in, double na_out, double na_in, double cl_out, double cl_in, double T = body_temperature)
        {


            //std::cout << "Log number: " << (p_k * k_out +p_na * na_out + p_cl * cl_in)  /  (p_k * k_in + p_na * na_in  + p_cl * cl_out) << std::endl;

            return ((gas_constant_R * body_temperature ) / faradays_constant_F) * log( (p_k * k_out +p_na * na_out + p_cl * cl_in)  /  (p_k * k_in + p_na * na_in  + p_cl * cl_out));

            
        }

        /*leak channel: represents Na entering the cell and K leaving the cell moving at a constant slow rate
        using a 3:2 ratio (1.5:1.0 in the code) because we want to counteract the Na/K pump which pumps 3 Na 
        out for 2 K in*/
        void apply_leak(double &na_in, double &k_in, double leak_rate)
        {
            na_in += 1.5 * leak_rate;
            k_in -= 1.0 * leak_rate;
        }

        /*represents the Na/K pump which actively pumps out 3 Na for every 2 K in
        Using a 3:2 ratio (1.5:1.0 in the code) because we want to counteract the leak channel
        which has the opposite affect as the pump*/
        void apply_pump(double &na_in, double &k_in, double pump_rate)
        {
            na_in -= 3 * pump_rate;
            k_in += 2 * pump_rate;
        }



        

};

//represents one ion
struct Ion
{
    std::string name;
    double conductance_out;
    double conductance_in;
    double z;   
};

int main(int argc, char* argv[])
{
    Neuron n;

    Ion Potassium;
    Potassium.name = "K+   ";
    Potassium.conductance_out = 5.0;
    Potassium.conductance_in = 140.0;
    Potassium.z = 1;

    Ion Sodium;
    Sodium.name = "Na+   ";
    Sodium.conductance_out = 145;
    Sodium.conductance_in = 15.0;
    Sodium.z = 1;

    Ion Cloride;
    Cloride.name = "Cl-   ";
    Cloride.conductance_out = 110.0;
    Cloride.conductance_in = 10.0;
    Cloride.z = -1;

    std::cout << "---------------------" << std::endl;
    
    std::cout << "Printing out K, Na, and Cl nernst potentials individually" << std::endl;
    std::cout << Potassium.name << 1000 * n.nernst_potential(Potassium.conductance_out,Potassium.conductance_in,Potassium.z) << " mV" << std::endl;
    std::cout << Sodium.name << 1000 * n.nernst_potential(Sodium.conductance_out,Sodium.conductance_in,Sodium.z) << " mV" << std::endl;
    std::cout << Cloride.name  << 1000 * n.nernst_potential(Cloride.conductance_out,Cloride.conductance_in,Cloride.z) << " mV" << std::endl;

    std::cout << "---------------------" << std::endl;

    std::cout << "Using GHK equation to test voltage of the cell given each nernst potential of the ions in the cell" << std::endl;

    //Case A
    std::cout << "Case A: " << 1000 *  n.ghk_voltage(50, 2, 15, Potassium.conductance_out, Potassium.conductance_in, Sodium.conductance_out, Sodium.conductance_in, Cloride.conductance_out, Cloride.conductance_in) << " mV" << std::endl;

    //Case B
    std::cout << "Case B: " << 1000 *  n.ghk_voltage(5,100, 5, Potassium.conductance_out, Potassium.conductance_in, Sodium.conductance_out, Sodium.conductance_in, Cloride.conductance_out, Cloride.conductance_in) << " mV" << std::endl;

    std::cout << "---------------------" << std::endl;

    std::cout << "Simulating leak and pump channels" << std::endl;

    for(int loop_variable = 0; loop_variable < 500; loop_variable++)
    {
        n.apply_leak(Sodium.conductance_in, Potassium.conductance_in, 5);
        n.apply_pump(Sodium.conductance_in, Potassium.conductance_in, 3);

        std::cout << "New Voltage: " << 1000 *  n.ghk_voltage(50, 0, 15, Potassium.conductance_out, Potassium.conductance_in, Sodium.conductance_out, Sodium.conductance_in, Cloride.conductance_out, Cloride.conductance_in) << " mV" << std::endl;

    }

    std::cout << Potassium.name << 1000 * n.nernst_potential(Potassium.conductance_out,Potassium.conductance_in,Potassium.z) << " mV" << std::endl;
    std::cout << Sodium.name << 1000 * n.nernst_potential(Sodium.conductance_out,Sodium.conductance_in,Sodium.z) << " mV" << std::endl;
    std::cout << Cloride.name  << 1000 * n.nernst_potential(Cloride.conductance_out,Cloride.conductance_in,Cloride.z) << " mV" << std::endl;
    std::cout << "Potassium inside concentration: " << Potassium.conductance_in << std::endl;



    return 0;
}