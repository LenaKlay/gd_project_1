// To compile:
// gcc stochwave.c -lm -o test

// To run:
// ./test 1000 0.575 1.25 0.1 10000 1 > simtext.csv
//         K     s    r   mig tfinal nrep

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> // to get system time

// Simulation parameters
#define DTIME 0.05 // Time interval at which population is updated (has to be <1)

#define MAXSTEPS 10000000000 // Maximum number of simulation steps (in case a loop goes wrong)
#define STEPPRINTINTERVAL 100 // Print result every ... steps 
#define NSITES 1000 // Number of sites, i.e. discrete spatial locations   // LK (100 before)


// Define variables
int nD[NSITES]; // number of DD individuals at each site
int nO[NSITES]; // number of OO individuals at each site
// Note: WT genotype written with letter O
int newnD[NSITES]; // temporary number of DD
int newnO[NSITES]; // temporary number of OO
int newFocal; // Whether the individual that we are considering survives
int newOffspring; // Number of offspring produced by the focal
double fecundity; // Expected number of offspring in the time interval
int nLocal; // Temporary number
int nLeft; // Temporary number (migration)
int nRight; // Temporary number (migration)
int nStay; // Temporary number (migration)

int totPopSize; // Total population size
int totDrive; // Total number of Drive alleles in the population
int totWT; // Total WT

int timestep; // Time Step
double currenttime; // Time value
double deltaTime; // Time interval, if death occurs within DTIME

int isite; // Counter
int iindiv; // Counter
double randnum; // Random number
double randDeathTime; // Random time

double pPois; // Temporary variable for Poisson sampling
double LPois; // Temporary variable for Poisson sampling

// Functions used in the script
void printPopulationState(void); // Print current state of the population
void oneStep(void); // Do one step of the simulation and update the population accordingly
void deathAndBirth(void); // substep: Death and Birth
void migrate(void); // substep: Migration

int g_argc; // Global variable for command line parameter entry
char **g_argv; // Global variable for command line parameter entry

// Introduce parameter names (values entered when running the script)
double K; // Carrying capacity
double s; // Zygote survival cost induced by the drive (0<= s <= 1)
double r; // Growth rate
double mig; // Emigration probability
double MAXTIME; // Maximum time until which simulation is run
int NREPS; // Number of simulation replicates

// Other parameters
double deathRate = 1.0; // Scaled parameter

/* Function to process the parameters, give them names */
void process_command_line(void)
{
  // Extract parameters from the command line
  /// Model parameters
  K = (double) strtod(g_argv[1], (char **) NULL); // K, carrying capacity
  s = (double) strtod(g_argv[2], (char **) NULL); // Zygote survival cost induced by the drive (0<= s <=1)
  r = (double) strtod(g_argv[3], (char **) NULL); // Growth rate
  mig = (double) strtod(g_argv[4], (char **) NULL); // Emigration probability
  /// Simulation parameters
  MAXTIME = (double) strtod(g_argv[5], (char **) NULL); // Maximal time value of the simulation
  NREPS = (int) strtod(g_argv[6], (char **) NULL); // Number of replicates
}

int main(int argc, char **argv)
{
  // Process command line arguments
  g_argc = argc;
  g_argv = argv;
  process_command_line();
  //printf("K = %.0f, s = %.4f, r = %.2f, mig = %.4f, MAXTIME = %.2f, NREPS = %d\n", K, s, r, mig, MAXTIME, NREPS); // LK

  // Random seed
  // There we set the random seed to ensure reproducibility
  unsigned long seed=98756411; // to be changed
  srand48(seed);

  int irep = 0; // REPLICATE number

  for(irep = 0; irep<((int)NREPS); irep++){

    // Initialize time
    currenttime = 0.0;

    // Initialize the population
      // WT on the left-hand side, at carrying capacity K
    for(isite=0; isite< floor(NSITES/2); isite++){               // LK j'inverse les CI !
      nD[isite] = (int) floor(K);
      nO[isite] = (int) 0;
    }
      // Drive on the right-hand side, at K too
      // NB: this is not the carrying capacity of an all-drive population when s>0
    for(isite=floor(NSITES/2); isite< NSITES; isite++){          // LK j'inverse les CI !
      nD[isite] = (int) 0;
      nO[isite] = (int) floor(K);
    }
      
    // Compute total population size and total drive
    totWT = 0;
    totDrive = 0;
    
    // If only one replicate, print population state                     //LK
    if(NREPS<2){ // If one sim, print output at regular time intervals 
      printPopulationState();                                          
    }                                                                
    
    for(isite=0; isite<NSITES; isite++){
      totWT += nO[isite];
      totDrive += nD[isite];
    }
    totPopSize = totWT + totDrive;
    
    for(timestep = 0; timestep<((int)MAXSTEPS); timestep++){
      // Check that the population is not extinct, that the drive is still present
      // and that time is not over
      if(totWT<1 || totDrive<1 || currenttime >= MAXTIME){
          break;
        }
      
      // Do one step of the simulation
      oneStep();

      // If only one replicate, print population state                            //LK
      if(NREPS<2){ // If one sim, print output at regular time intervals
        if(1.0*timestep - STEPPRINTINTERVAL * floor(timestep/STEPPRINTINTERVAL) < 1.0){
          // if round number of timestep in base STEPPRINTINTERVAL, print ouput
          printPopulationState();
        }
      }
    }
    
    // Simulation finished -> save final outcome
    printPopulationState();
  } //end rep
  
  return(0);

} // end main Functions



/*----------------------------------------------------------------------------*/
/* ONE SIMULATION STEP - ZYGOTE CONVERSION */
void oneStep(void){
  // Update time
  currenttime += (1.0)*DTIME;
  
  //-------------------------------------------
  // 1) Death and birth
  // Update each site one at a time
  for(isite=0; isite<NSITES; isite++){
    // Reset newpop
    newnD[isite] = 0;
    newnO[isite] = 0;

    // WT individuals
    for(iindiv=0; iindiv<nO[isite]; iindiv++){
      // Compute expected number of offspring
      fecundity = (1.0 + r * (1.0 - (1.0*(nD[isite] + nO[isite]))/(1.0*K))) * (nO[isite]) / (nO[isite] + nD[isite]);
      // Rescale it - can be negative because number of individuals can go over K (migration, independent reproduction)
      if(fecundity < 0.0){
        fecundity = 0.0;
      }
      deathAndBirth();
      newnO[isite] += newFocal + newOffspring;
    }
    
    // Drive individuals
    for(iindiv=0; iindiv<nD[isite]; iindiv++){
      // Compute expected number of offspring in the time interval
      fecundity = (1.0 + r * (1.0 - (1.0*(nD[isite] + nO[isite]))/(1.0*K))) * (1.0 - s) * (2.0 * nO[isite] + nD[isite]) / (nO[isite] + nD[isite]); // Expected number of offspring
      // Rescale it
      if(fecundity < 0.0){
        fecundity = 0.0;
      }
      deathAndBirth();
      newnD[isite] += newFocal + newOffspring;
    }
  }

  
  //--------------------------------------------
  // 2) Migration
  
  // Reset pop
  for(isite=0; isite<NSITES; isite++){
    nD[isite] = 0;
    nO[isite] = 0;
  }
  
  // Each site, excluding edges
  for(isite=1; isite<(NSITES-1); isite++){
    // Drive individuals
    nLocal = newnD[isite];
    migrate();
    nD[isite-1] += nLeft;
    nD[isite+1] += nRight;
    nD[isite] += nStay;
    // WT individuals
    nLocal = newnO[isite];
    migrate();
    nO[isite-1] += nLeft;
    nO[isite+1] += nRight;
    nO[isite] += nStay;
  }
  // Edges
  // Site 0 -> to the left in fact stay home
    // Drive individuals
    nLocal = newnD[0];
    migrate();
    nD[0] += nLeft + nStay;
    nD[1] += nRight;
    // WT individuals
    nLocal = newnO[0];
    migrate();
    nO[0] += nLeft + nStay;
    nO[1] += nRight;

  // Site (NSITES-1) -> to the right in fact stay home
    // Drive individuals
    nLocal = newnD[NSITES-1];
    migrate();
    nD[NSITES-1] += nRight + nStay;
    nD[NSITES-2] += nLeft;
    // WT individuals
    nLocal = newnO[NSITES-1];
    migrate();
    nO[NSITES-1] += nRight + nStay;
    nO[NSITES-2] += nLeft;

  //--------------------------------------------
  // Updates
  totWT = 0;
  totDrive = 0;
  
  for(isite=0; isite<NSITES; isite++){
    totDrive += nD[isite];
    totWT += nO[isite];
  }
  totPopSize = totDrive + totWT;
}

/*----------------------------------------------------------------------------*/
/* SUBSTEP: DEATH AND BIRTH */
void deathAndBirth(void){
  // !! Need to have specified fecundity before (depends on genotype)
  
  // Sample death time
  randnum = drand48(); // Generate uniform random number
  randDeathTime = -1.0/(deathRate) * log(1.0 - randnum); // Convert into exponential random number
  	
  if(randDeathTime < DTIME){ // If the individual dies within DTIME
    deltaTime = randDeathTime; // We are now considering the time until the focal dies
    newFocal = 0; // Individual does not survive
  }else{ // The individual has not died withing DTIME
    deltaTime = DTIME; // Do not change end time of the step
    newFocal = 1; // Individual survives
  };
  
  newOffspring = 0;
  if(fecundity > 0.0){
    // Sample number of offspring produced during the interval (Poisson)
    // Source: https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
    // (Lambda is below 1 so we are fine with the simple algorithm)
    LPois = exp(-fecundity*deltaTime);
    pPois = 1.0;
    while (pPois >= LPois) {
      newOffspring += 1;
      randnum = drand48();
      pPois *= randnum;
    }
    newOffspring += -1;
  }
  // output is (newFocal, newOffspring)
}

/*----------------------------------------------------------------------------*/
/* SUBSTEP: MIGRATION */
void migrate(void){
  // For a specified site and a specified type, nLocal
  
  // Reset output variables
  nLeft = 0;
  nRight = 0;
  nStay = 0;
  
  // Then for each individual... (binomial sampling)
  for(iindiv=0; iindiv<nLocal; iindiv++){
    randnum = drand48(); // Random uniform to decide emigration
    if(randnum < mig){
      // Emigrates
      randnum = drand48(); // Random uniform to pick direction
      if(randnum < 0.5){
        nLeft += 1;
      }else{
        nRight += 1;
      }
    }else{ // Does not emigrate
      nStay += 1;
    }
  }
}

/*----------------------------------------------------------------------------*/
/* IF ONE REP: PRINT POPULATION STATE */
/* PRINT STATE OF THE POPULATION */
void printPopulationState(void){
  // Print current time
  printf("%09.3f, ", currenttime);
  // Print number of individuals of each genotype
  for (isite=0; isite<(NSITES-1); isite++) {
    printf("%04d, ", nO[isite]);
    printf("%04d, ", nD[isite]);
  }
  printf("%04d, ", nO[NSITES-1]);
  printf("%04d\n", nD[NSITES-1]);
}
