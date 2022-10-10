#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


/////////////////////////////
/* Defining RNG parameters */
/////////////////////////////

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long int seed_start = 5;
long int * seed;


/////////////////////////
/* Defining path array */
/////////////////////////

double * y;


/////////////////////////////////////
/* Declaring simulation parameters */
/////////////////////////////////////

double T;
int N; // Discretization steps
double eta; // Lattice step in units of 1/omega, omega is the oscillator frequency
int init_flag; // Starting state of the matrix, either hot (1), cold (0) or previous final state
int measures; // Number of measures
int decorrel_len; // Decorrelation length
int thermal_len; // Number of measures to be discarded
double delta_metro; // Maximum variation during Metropolis function
int bins; // Number of bins in histogram
double y_max; // Defines range for gaussian fit

/* Pointers to the adjacent sites */
int * npp;
int * nmm;


///////////////
/* I/O files */
///////////////


ifstream input_Parameters;
ifstream input_lattice;
ofstream output_lattice;
ofstream output_means;
ofstream output_probs;
ofstream output_green2;
ofstream output_green4;

/* Simulation flag */

int sim_flag;

/////////////////////////
/* Declaring functions */
/////////////////////////

float Ran2(long *idum); // Random number generator
void Geometry(int * movePlus, int * moveMinus, int d); // Generates proximity arrays
void LatticeInit(double * array, int flag, long int * seed); // Initializes array in a state defined by iflag
void Metropolis(double * array, long int * seed, double delta, int deco); // Generates Markov chain
double MeanValue(double * array, int len_array);
double Mean2(double * array, int len_array);
double MeanDiff2(double * array, int len_array);
double Green2(double * array, int len_array, int dist);
double Green4(double * array, int len_array, int dist);
void Hist(double * array_in, double * array_out, int len_array, int bin);


////////////////////////////
/*Declare output variables*/
////////////////////////////

double mean;
double mean2;
double mean_delta2;
double * green2;
double * green4;
double * probs;


//////////////////
/* Main program */
//////////////////

int main() {


    /* Reading the simulation parameters from parameters.txt */
    input_Parameters.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/input/parameters.txt", ios::in); // Opening file
    if (input_Parameters.is_open()) { // Check if file is open, then read

        string line;
        string::size_type sz;

        getline(input_Parameters, line);
        T = stod(line,&sz); // T is actually the inverse temperature, we called it T for simplicity
        getline(input_Parameters, line);
        measures = stoi(line,&sz);
        getline(input_Parameters, line);
        N = stoi(line,&sz);
        getline(input_Parameters, line);
        decorrel_len = stoi(line,&sz);
        getline(input_Parameters, line);
        thermal_len = stoi(line,&sz);
        getline(input_Parameters, line);
        init_flag = stoi(line,&sz);
        getline(input_Parameters, line);
        getline(input_Parameters, line);
        bins = stoi(line,&sz);
        getline(input_Parameters, line);
        y_max = stod(line,&sz);

        input_Parameters.close();

    }
    else { // Error message
        cerr << "Unable to open parameters file.\n";
    }


    /* Derive eta and delta_metro from inverse temperature and N. eta is the lattice step in units of 1/omega*/
    eta = (double)(T/(double)(N));
    delta_metro = sqrt(eta);

    /* Ask the user what type of simulation they want to run */
    cout << "This program simulates a harmonic oscillator.\n";
    cout << "Since it is very resource intensive, you will need to select which type of quantities you need.\n";
    cout << "Enter (0) to run a simulation that calculates the energy of the system at a given temperature.\n";
    cout << "Enter (1) to run a simulation that calculates the wave function of the ground state (remember to run this at low temperature).\n";
    cout << "Enter (2) to run a simulation that calculates the value of the first two energy gaps.\n";
    cin >> sim_flag;


    /* Initialize seed for RNG */
    seed = &seed_start;


    /* Create proximity arrays */
    int npp_array[N];
    int nmm_array[N];
    npp = npp_array;
    nmm = nmm_array;


    /* Compute proximity conditions */
    Geometry(npp,nmm,1);


    /* Allocating space for the 1-dimensional lattice */
    y = new double [N];


    /* Allocating space for arrays of the connected Green functions */
    green2 = new double[N/2];
    green4 = new double[N/2];

    /* Allocating space for probability function */
    probs = new double[bins];


    /* Initialize lattice */
    input_lattice.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/input/lattice.txt", ios::in);
    if(input_lattice.is_open()){ // Check if file is open, then run
        LatticeInit(y, init_flag, seed);
        input_lattice.close();
    }
    else { // Error message
        cerr << "Unable to open Lattice file.\n";
    }
    input_lattice.close();


    /* Check simulation flag and whether output files are open then run algorithm */
 
    if(sim_flag==0){

        output_means.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/means.txt", ios::trunc);
    
        if(output_means.is_open()){

            /* Start Markov chain and take a measurement for each iteration */
            for(int l=0;l<measures;l++){


                /* Call Metropolis function decorrel_len times before taking a measurement*/

                Metropolis(y,seed,delta_metro,decorrel_len);


                /* Computing output variables */
                mean = MeanValue(y,N);
                mean2 = Mean2(y,N);
                mean_delta2 = MeanDiff2(y,N);


                /* Writing array onto output files */
                output_means << mean << "\t" << mean2 << "\t" << mean_delta2 << "\n";

            }

            output_means.close();

        }

        else{// print error message
            cerr << "Unable to open output files.\n";
        }
        
    }


    if(sim_flag==1){
        
        /* We store the number of times the oscillator visits an interval in this file, each row represents a measure and each column an interval */
        output_probs.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/probs.txt", ios::trunc);

        if(output_probs.is_open()){


            /* Start Markov chain and take a measurement for each iteration */
            for(int l=0;l<measures;l++){


                /* Call Metropolis function decorrel_len times before taking a measurement*/
                Metropolis(y,seed,delta_metro,decorrel_len);

                /* Call Hist funtion to generate probability data */
                Hist(y, probs, N, bins);

                /* Writing array onto output files */
                for(int s=0;s<bins;s++){
                    output_probs << probs[s] << "\t";
                }
                output_probs << "\n";

            }

            /* Closing output files */
            output_probs.close();
        
        }

        else{// print error message
            cerr << "Unable to open output file.\n";
        }
        
    }


    if(sim_flag==2){

        output_means.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/means.txt", ios::trunc);
        output_green2.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/green2.txt", ios::trunc);
        output_green4.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo3/_data/green4.txt", ios::trunc);

        if(output_means.is_open() && output_green2.is_open() && output_green4.is_open()){

            /* Start Markov chain and take a measurement for each iteration */
            for(int l=0;l<measures;l++){


                /* Call Metropolis function decorrel_len times before taking a measurement*/
                Metropolis(y,seed,delta_metro,decorrel_len);



                /* Computing output variables */
                mean = MeanValue(y,N);
                mean2 = Mean2(y,N);

                for(int k=0;k<N/2;k++){
                    green2[k] = Green2(y,N,k+1);
                    green4[k] = Green4(y,N,k+1);
                }


                /* Writing array onto output files */
                output_means << mean << "\t" << mean2 << "\n";


                for(int k=0;k<N/2;k++){
                    output_green2 << green2[k] << "\t";
                    output_green4 << green4[k] << "\t";
                }
                output_green2 << "\n";
                output_green4 << "\n";
            }

            /* Closing output files */
            output_means.close();
            output_green2.close();
            output_green4.close();

        }
        
        else{// print error message
            cerr << "Unable to open output files.\n";
        }

    }


    /* Deallocate memory */
    delete[] y;
    delete[] green2;
    delete[] green4;
    delete[] probs;

    return 0;
}


///////////////
/* Functions */
///////////////

/* The Ran2 function is shamelessly taken from "Numerical recipes in C" */
float Ran2(long *idum){
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

/* If d=1, this function creates two arrays that allow us to select the adjacent sites of a certain point, instead of computing them every time, which is more computationally demanding, while defining the "toroid" boundary conditions. */
void Geometry(int * movePlus, int * moveMinus, int d){
    for(int count=0; count<N; count++){
        *(movePlus + count) = (count+d)%(N);
        *(moveMinus + count) = (N+count-d)%(N);
    }
  return;
}

/* The following function creates a matrix whose configuration depends on the value of init_flag */
void LatticeInit(double * array, int flag, long int * seed){

    /* Initialize array whose elements are all 0.0 */
    if(flag==0){
        for(int k=0; k<N; k++){
            array[k] = 0.0;
        }
        array[0]=array[N-1];
    }

    /* Initialize array whose elements are random, uniformly distributed inside standard deviation*/
    else if(flag==1){
        for(int k=0; k<N; k++){
        array[k] =2*(Ran2(seed)-0.5)*sqrt(eta);
        }
        array[0]=array[N-1];
    }

    /* Initialize array whose elements are read from the last iteration of the algorithm*/
    else {
        for(int k=0; k<N; k++){
            string line;
            string::size_type sz;
            getline(input_lattice, line);
            array[k]= stod(line,&sz);
        }
    }

    return;
}

/* The Metropolis function defines the Markov chain. */
void Metropolis(double * array, long int * seed, double delta, int deco){

    for(int r=0;r<deco;r++){
        for(int i=0; i<(N); i++){
            // Using the results of Geometry to select the adjacent sites
            int ip = *(npp + i);
            int im = *(nmm + i);

            /* This is where "physics" kicks in, we need to compute the action of the system and update the array based on that value*/
            double yp = array[i]+2*(Ran2(seed)-0.5)*delta;

            double action_diff = ((double)(1.0/eta)+(double)(eta/2.0))*(yp*yp - array[i]*array[i]) + (array[ip]+array[im])*(double)((array[i]-yp)/eta);
            // Last thing we need is the Metropolis test
            double u = Ran2(seed);
            double p = (double)(exp(-action_diff));

            if(u<p){
                array[i]=yp;
            }
        }
    }
  return;
}

/* The Hist function is basically a binning function used to collect the data needed to fit the ground state of the oscillator. The range is selected by the user and read from parameters.txt */
void Hist(double * array_in, double * array_out, int len_array, int bin){

    double array_in_min = -y_max;
    double array_in_max = y_max;

    double range = array_in_max - array_in_min;
    double step = (double)(range/(double)(bin));

    for(int i = 0; i < bin; i++){
        array_out[i] = 0.0;
        for(int j = 0; j < len_array; j++){
            if(array_in[j] >= array_in_min + i*step && array_in[j] < array_in_min + (i+1)*step){
                array_out[i] = array_out[i] + 1.0; // Function normalized during data analysis
            }
        }
    }

    return;

}

/* Computes mean value of an array */
double MeanValue(double * array, int len_array){
    double sum = 0.0;
    for(int i=0;i<len_array;i++){
        sum = sum + array[i];
    }
    return (double)(sum/(double(len_array)));
}

/* Computes mean value of the squared elements of an array */
double Mean2(double * array, int len_array){
    double sum = 0.0;
    for(int i=0;i<len_array;i++){
        sum = sum + array[i]*array[i];
    }
    return (double)(sum/(double(len_array)));
}

/* Computes mean value of the squared difference of consecutive elements of an array */
double MeanDiff2(double * array, int len_array){
    double sum=0.0;
    for(int i=0;i<len_array;i++){
        int ip = *(npp + i);
        sum = sum + (array[ip]-array[i])*(array[ip]-array[i]);
    }
    return (double)(sum/(double(len_array)));
}

/* Computes mean value of the correlation between elements of an array which are separated by dist elements */
double Green2(double * array, int len_array, int dist){
    double sum=0.0;
    for(int i=0;i<len_array;i++){
        int i_dist = (i+dist)%(len_array);
        sum = sum + array[i_dist]*array[i];
    }
    return (double)(sum/(double(len_array)));
}

/* Computes mean value of the correlation between squared elements of an array which are separated by dist elements */
double Green4(double * array, int len_array, int dist){
    double sum=0.0;
    for(int i=0;i<len_array;i++){
        int i_dist = (i+dist)%(len_array);
        sum = sum + array[i_dist]*array[i]*array[i_dist]*array[i];
    }
    return (double)(sum/(double(len_array)));
}
