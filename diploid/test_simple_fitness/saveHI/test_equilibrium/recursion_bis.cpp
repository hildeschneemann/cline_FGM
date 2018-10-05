#include "fisher.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
#include "H5Cpp.h"
using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;
extern MTRand rnd;


/*----------------------------------------------------------
Function recursion: simulates a 1D stepping-stone, with a barrier 
appearing after T1 generations, followed
by a secondary contact after T1 + ff1 generations.
 
Parameters:

dv: number of demes
Nv: size of the initial population
migv: migration rate (between adjacent demes)
bv: position of the barrier
propv: number of individuals in the first deme (after the split)
nv: number of phenotypic dimensions (complexity)
mv: number of dimensions affected by a mutation (pleiotropy)
sigv: size of mutational steps along each phenotypic axis
diffv: change in optimum after split
av: strength of selection (steepness of fitness function)
Qv: curvature of fitness function (cf Gros et al 2009 Genetics)
Uv: total mutation rate per genome (on selected loci)
nbSv: total number of selected loci per genome
Lv: genome map length (mean nb of cross-overs per meiosis)
	set to -1 for free recombination
HI: hybrid index of focal individual
p12: heterozygosity of focal individual 
-----------------------------------------------------------*/


void recursion(int dv, int Nv, double migv, int bv, int nv, int mv, double sigv, double av, double diffv, double Qv, double Uv, int nbSv, double Lv,
			   int T1v, int ff1v, int fpv, int rv, int pasv, int nov)
{
	// variables:

	int i, j, k, loc, gen, mut, par1, par2, ind, nb, nb1, nb2, nb3, nb4, nbMig, ns, part, nbCo, nbSign;
	double w, wbar, varw, rd, pp, d, x, sz2, HI, p12, p2, delta_p, p_old;
	vector<int> store;

	// various fixed quantities:
	
    int Nd = Nv * dv;
    int nd = nv * dv;
    int twoN = 2*Nv;
    int twoNd = 2*Nv*dv;
	int N_1 = Nv - 1;
	int n_1 = nv - 1;
	int nbS_1 = nbSv - 1;
	int nS = nbSv * nv;
	int Nn = Nv * nv;
	double hQ = Qv / 2.0;
    int NbGen = T1v + ff1v + fpv; // total number of generations
	int NbGen_1 = NbGen - 1;
    int Tcontact = T1v + ff1v; // time of secondary contact
	int Nd1 = twoN * bv;
	bool sign = true;
	bool equi = false;

	boost::dynamic_bitset<> tmp1;
	boost::dynamic_bitset<> tmp2;
	boost::dynamic_bitset<> tmp3;
	boost::dynamic_bitset<> tmp4;

	
		// HDF5 constants and variables

	string fileName;
	stringstream nameF;
	nameF << "result_d" << dv << "_N" << Nv << "_mig" << migv << "_b" << bv
		<< "_n" << nv << "_m" << mv << "_sig" << sigv << "_diff" << diffv
        << "_a" << av << "_Q" << Qv << "_U" << Uv << "_nbS" << nbSv
		<< "_L" << Lv << "_Ts" << T1v << "_ff1" << ff1v << "_fp" << fpv
		<< "_r" << rv << "_" << nov << ".h5";
	nameF >> fileName;

	const H5std_string FILE_NAME(fileName);
	const H5std_string DSET_FREQ_NAME("alleleFreq");
	const H5std_string DSET_W_NAME("w");
	const H5std_string DSET_GEN_NAME("savedGen");
	const H5std_string DSET_HI_NAME("finalHI");

	const int RANK_FREQ = 3; // dimensions of the dataset
	const int RANK_W = 2;
	const int RANK_HI = 2;
	const int DIM0_SUB = 1;	// subset dimensions
	const int DIM1_SUB = nbSv;
	const int DIM2_SUB = dv;
	const int DIM0 = NbGen / pasv + 1; // size of dataset frequencies
	const int DIM1 = nbSv;
	const int DIM2 = dv;
	const int DIM_W0 = DIM0; // size of dataset fitness
	const int DIM_W1 = dv;
	const int DIM_HI = Nd;
	hsize_t dim_sub_freq[RANK_FREQ];
	hsize_t dim_sub_w[RANK_W];
	int indexGen = 0; // where to write in HDF5 dataset
	double* sdata_freq = new double[DIM0_SUB*DIM1_SUB*DIM2_SUB]; // subset to write to HDF5 dataset
	double* sdata_w = new double[1*DIM_W1];
	double* sdata_hi = new double[DIM_HI];
	int* savedGen = new int[DIM0];
	
	double* freq_equi = new double[nbSv*3*10]; //check whether equilibrium is reached
	double* sdata_freq_old = new double[nbSv*3*10];
    // population: table of 2N*d chromosomes (two chromosomes per individual):

    chr * pop = new chr [twoNd];
	chr * temp = new chr [twoNd]; // used to generate next generation
    chr * cp;

	// "mutations" will hold the effect of the 1 allele at each locus
	// on each phenotypic axis:
	
	double * mutations = new double [nS];

	// "Wtot" will hold the fitness of each individual:

	double * Wtot = new double [Nd];
    
    // maximal fitness in each deme:
    
    double * Wmax = new double [dv];
	
	// tables for means and variances of phenotypic traits:

	double * m = new double [nv];
	double * v = new double [nv];

	// "HIend" will hold the hybrid index of each individual at the last generation
	
	double * HIend = new double [Nd];
	
	// create HDF5 file and open it

	createHDF5(FILE_NAME, DSET_FREQ_NAME, RANK_FREQ, DIM0, DIM1, DIM2,
		DSET_W_NAME, RANK_W, DIM_W0, DIM_W1, DSET_HI_NAME, DIM_HI, DSET_GEN_NAME);

	H5::H5File file;
	H5::DataSet dset_freq;
	H5::DataSet dset_w;
	H5::DataSet dset_hi;
	H5::DataSet dset_gen;
	file.openFile(FILE_NAME, H5F_ACC_RDWR);
	dset_freq = file.openDataSet(DSET_FREQ_NAME);
	dset_w = file.openDataSet(DSET_W_NAME);
	dset_hi = file.openDataSet(DSET_HI_NAME);
	dset_gen = file.openDataSet(DSET_GEN_NAME);

	dim_sub_freq[0] = DIM0_SUB;
	dim_sub_freq[1] = DIM1_SUB;
	dim_sub_freq[2] = DIM2_SUB;
	dim_sub_w[0] = 1;
	dim_sub_w[1] = DIM_W1;
	

	// write parameters as attributes of alleleFreq dataset
	//(int dv, int Nv, double migv, int bv, int nv, int mv, double sigv, double av, double diffv, double Qv, double Uv, int nbSv, double Lv, int T1v, int ff1v, int fpv, int pasv, int nov)
	writeAttributeHDF5(file, dset_freq, "d", dv);
	writeAttributeHDF5(file, dset_freq, "N", Nv);
	writeAttributeHDF5(file, dset_freq, "mig", migv);
	writeAttributeHDF5(file, dset_freq, "b", bv);
	writeAttributeHDF5(file, dset_freq, "n", nv);
	writeAttributeHDF5(file, dset_freq, "m", mv);
	writeAttributeHDF5(file, dset_freq, "sig", sigv);
	writeAttributeHDF5(file, dset_freq, "a", av);
	writeAttributeHDF5(file, dset_freq, "diff", diffv);
	writeAttributeHDF5(file, dset_freq, "Q", Qv);
	writeAttributeHDF5(file, dset_freq, "U", Uv);
	writeAttributeHDF5(file, dset_freq, "nbS", nbSv);
	writeAttributeHDF5(file, dset_freq, "L", Lv);
	writeAttributeHDF5(file, dset_freq, "T1", T1v);
	writeAttributeHDF5(file, dset_freq, "ff1", ff1v);
	writeAttributeHDF5(file, dset_freq, "fp", fpv);
	writeAttributeHDF5(file, dset_freq, "pas", pasv);
	writeAttributeHDF5(file, dset_freq, "no", nov);

	
	// initialization: allele 0 is fixed at all selected loci:

	for (i = 0; i < twoNd; i++)
    {
		pop[i].sel.resize(nbSv);
        temp[i].sel.resize(nbSv);
    }

		// initializes complete differentiation
	for (i = 0; i < Nd1; i++)
	{
		    for (j = 0; j < nbS_1; j++)
				pop[i].sel.flip(j);
                temp[i].sel.flip(j);
	}


	

    
    // generations:
gen=0;	
	while (equi==false & gen < NbGen)
	{ 
		// fitness of each individual, maximal fitnesses,
		// mean fitness and variance in fitness:
		
		wbar = 0;
		varw = 0;
        //for (k = 0; k < nv; k++)
        //{
        //    m[k] = 0;
        //    v[k] = 0;
        //}
        
		for (i = 0; i < dv; i++)  // for each deme
        {
            nb = i * Nv;
            Wmax[i] = 0;
            
            for (j = 0; j < Nv; j++)  // for each individual
            {
                nb2 = 2 * (nb + j);
                //sz2 = 0;
				// determine hybrid index and heterozygosity
				tmp1 = (pop[nb2].sel ^ pop[nb2+1].sel);
				p12 = tmp1.count(); 
				p12 /= nbSv; //creates a vector with a 0 for homozygote and 1 for heterozygote sites
				tmp2 = (pop[nb2].sel & pop[nb2+1].sel);
				p2 = tmp2.count();
				p2 /= nbSv;
				HI = p2+ p12/2; //so how do I determine parental alleles? do I simply say 0 vs 1 ? then I need to start with diverged pop and have no mutation I guess

				sz2 = fpv + (4 - 2*fpv)*4*HI*(1-HI)+(ff1v-1)*p12+(ff1v-1+rv)*p12*(1-p12);
                // fitness
                
                w = exp(-av * pow(sz2,hQ));
                Wtot[nb2/2] = w;
		//if (gen == NbGen_1)
			HIend[nb2/2] = HI;

                wbar += w;
                varw += w * w;
                if (Wmax[i] < w)
                    Wmax[i] = w;
            }
        }
		wbar /= Nd;
		varw /= Nd;
        //for (k = 0; k < nv; k++)
        //{
            //m /= Nd;
            //v /= Nd;
            //v -= (m * m);
        //}
        
        // sampling the next generation:
        
        for (i = 0; i < dv; i++) // for each deme i
        {
            nb3 = Nv*i;
            
            // number of immigrants in deme i:
            
            nbMig = int(binldev(migv, Nv));
            
            // migrant individuals:
            
            for (ind = 0; ind < nbMig; ind++)
            {
                nb4 = nb3 + ind;
                
                // selection of a deme of origin (j):
                

                    if ((i > 0) && (i < dv-1))
                    {
                        if (rnd.rand() < 0.5)
                            j = i - 1;
                        else
                            j = i + 1;
                    }
                    else
                    {
                        if (i == 0)
                            j = 1;
                        else if (i == dv - 1)
                            j = dv - 2;
                    }
                
                

                
                nb2 = Nv*j;
                
                // sampling first parent:
                
                do
                {
                    par1 = nb2 + int(rnd.randInt(N_1));
                    
                } while (rnd.rand()> (Wtot[par1]/Wmax[j]));
                
                // recombination
                
                rec(temp[2*nb4], pop[2*par1], pop[2*par1+1], Lv, nbSv);
                
                // sampling second parent:
                
                do
                {
                    par2 = nb2 + int(rnd.randInt(N_1));
                        
                } while (rnd.rand()> (Wtot[par2] / Wmax[j]));
                
                // recombination
                
                rec(temp[2*nb4+1], pop[2*par2], pop[2*par2+1], Lv, nbSv);
            }
            
            // philopatric individuals:
            
            for (ind = nbMig; ind < Nv; ind++)
            {
                nb4 = nb3 + ind;
                
                // sampling first parent:
                
                do
                {
                    par1 = nb3 + int(rnd.randInt(N_1)) ;
                    
                } while (rnd.rand()> (Wtot[par1] / Wmax[i]));
                
                // recombination
                
                rec(temp[2*nb4], pop[2*par1], pop[2*par1+1], Lv, nbSv);
                
                do
                {
                    par2 = nb3 + int(rnd.randInt(N_1));
                        
                } while (rnd.rand()> (Wtot[par2] / Wmax[i]));
                    
                // recombination
                
                rec(temp[2*nb4+1], pop[2*par2], pop[2*par2+1], Lv, nbSv);
            }
        }
        
        // mutation
        
        //for (i = 0; i < twoNd; i++)
        //{
        //    mut = int(poisdev(Uv)); // number of new mutations
        //    for (j = 0; j < mut; j++)
        //        temp[i].sel.flip(int(rnd.randInt(nbS_1)));
        //}

        // update population:

        cp = pop;
        pop = temp;
        temp = cp;
			

		// write result in HDF5 file
		if (gen % pasv == 0 || gen == NbGen_1)
		{

			// Allele frequency data --> THIS I HAD TO CHANGE FOR THE DIPLOID VERSION Nv --> twoN
			for (loc = 0; loc < nbSv; loc++)
			{
				for (i = 0; i < dv; i++)
				{
					nb = twoN * i;
					d = 0;
					for (j = 0; j < twoN; j++)
						if (pop[nb + j].sel[loc] == 1)
							d += 1;
							sdata_freq[loc * dv + i] = d / twoN;
							if (i == bv & loc == 1)
								delta_p = d / twoN - p_old;
								p_old = d / twoN;
						}					

					}
			//count how often the sign of the change in allele frequency changes
			if (delta_p ==0)
				nbSign += 1;
			if (delta_p < 0 & sign ==true)
				nbSign +=1;
				sign = false;
			if (delta_p > 0 & sign==false)
				nbSign+=1;
				sign=true;
			if (nbSign > 10)
				equi=true;			

			writeTimeStepHDF5(file, dset_freq, RANK_FREQ,
				dim_sub_freq, sdata_freq, indexGen);

			// mean fitness of demes --> THIS TOO!
			for (i = 0; i < dv; i++) // demes
			{
				sdata_w[i] = 0;
				for (j = 0; j < Nv; j++) // individuals
				{
					 sdata_w[i] += Wtot[i * Nv + j];
				}
				sdata_w[i] /= Nv;
			}
			writeTimeStepHDF5(file, dset_w, RANK_W,
				dim_sub_w, sdata_w, indexGen);

			// which generation was saved
			savedGen[indexGen] = gen + 1;

			indexGen +=	1;
		}
		
	gen+=1;
	} // end gen loop
	writeHISaved(file, dset_hi, HIend);

	// write savedGen to file
	writeGenSaved(file, dset_gen, savedGen);

	dset_freq.close();
	dset_w.close();
	dset_hi.close();
	dset_gen.close();
	file.close();


	// writes in output file:
	fprintf(fichierS, "\n\nResultats dans fichier ");
	fprintf(fichierS, "%s", fileName.c_str());
	fprintf(fichierS, "\n");
	         
		
	
	delete [] pop;
	delete [] temp;
	delete [] mutations;
	delete [] Wtot;
    delete [] Wmax;
	delete [] m;
	delete [] v;
}
