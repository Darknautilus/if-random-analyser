#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "von_neumann.h"
#include "aes.h"
#include "mersenne_twister.h"

#define ARRAY_MAX_SIZE 1000
#define OLDRAND_MAX 2147483647

#define VALS_VN "vals_von_newman"
#define VALS_MT "vals_mers_twist"
#define VALS_AES "vals_aes"
#define VALS_RAND_FAIBLE "vals_rand_faible"
#define VALS_RAND_FORT "vals_rand_fort"
#define NB_TIRAGES 1024

static unsigned int next;

int rdtsc()
{
	// cette fonction suivante cree un warning : c'est normal.
	__asm__ __volatile__("rdtsc");
}

void oldinit_rand(unsigned int seed)
{
	next = seed;
}

unsigned int oldrand()
{
	next = next * 1103515245 + 12345;
	return (unsigned int)(next % OLDRAND_MAX);
}

double Frequency(int nb_elem, word32 *tableau, int word_size)
{
    int S=0;
	double Sobs;
	double Pvaleur;
	int i;
	int n;
	
	
	n=nb_elem*word_size;
	
	for(i=0; i<nb_elem; i++)
	{
		int element=tableau[i];
		int j;
		for(j=0;j<word_size;j++)
		{
			S+=2*(element&0x1)-1;
			element=element>>1;
		}
	}
	Sobs=abs((double)S)/sqrt((double)n);
	Pvaleur=erfc(Sobs/sqrt(2));
	return Pvaleur;
}


int main()
{
	word16 x=1111; // nombre entre 1000 et 9999 pour Von Neumann
	struct mt19937p mt; // Pour Mersenne-Twister
	int tmp = rand(), seed; // Pour Mersenne-Twister
	u32 Kx[NK], Kex[NB*NR], Px[NB]; // pour l'AES

	unsigned int output_rand; // sortie du rand du C	
	word32 output_AES; // sortie pour l'AES
	word16 output_VN; // sortie pour pour Von Neumann
	word32 output_MT; // sortie pour Mersenne-Twister

               
	// initialisation des graines des generateurs

	srand(rdtsc());  // rand du C 
	seed = rand();
	oldinit_rand(seed);
	sgenrand(time(NULL)+(tmp), &mt); // Mersenne-Twister
	// Initialisation de la clé et du plaintext pour l'AES 
	// 45 est un paramètre qui doit changer à chaque initialisation
	init_rand(Kx, Px, NK, NB, 45);
	KeyExpansion(Kex,Kx); // AES : sous-clefs

	FILE *f_von_newman = fopen(VALS_VN,"w");
	FILE *f_mers_twist = fopen(VALS_MT,"w");
	FILE *f_aes = fopen(VALS_AES,"w");
	FILE *f_rand_faible = fopen(VALS_RAND_FAIBLE,"w");
	FILE *f_rand_fort = fopen(VALS_RAND_FORT,"w");
	
	word32 t_von_newman[NB_TIRAGES];
	word32 t_mers_twist[NB_TIRAGES];
	word32 t_aes[NB_TIRAGES];
	word32 t_rand_faible[NB_TIRAGES];
	word32 t_rand_fort[NB_TIRAGES];

	int i;
	for(i=0;i<NB_TIRAGES;i++)
	{
		// sorties des generateurs
		
		//output_rand = oldrand(); // rand du C
		output_rand = rand();
		output_VN = Von_Neumann(&x);
		output_MT = genrand(&mt);
		output_AES = AES(Px, Kex);
		
		fprintf(f_rand_faible,"%d\n",(output_rand&0xF)); // bits de poids faible
		fprintf(f_rand_fort,"%d\n",(output_rand>>27)); // bits de poids fort
		fprintf(f_von_newman,"%d\n",output_VN); // Von Neumann
		fprintf(f_mers_twist,"%d\n",output_MT); // Mersenne-Twister
		fprintf(f_aes,"%d\n",output_AES); // AES
		
		t_von_newman[i] = (word32)output_VN;
		t_mers_twist[i] = output_MT;
		t_aes[i] = output_AES;
		t_rand_faible[i] = (output_rand&0xF);
		t_rand_fort[i] = (output_rand>>27);
	}

	// affichage
	printf("- Generation de nombres aleatoires -\n");
	printf("Von Neumann written in : %s\n",VALS_VN);
	printf("Mersenne Twister written in : %s\n",VALS_MT);
	printf("AES written in : %s\n",VALS_AES);
	
	printf("Frequencies :\n");
	printf("Von Neumann : %f\n",Frequency(NB_TIRAGES,t_von_newman,16));
	printf("Mersenne Twister : %f\n",Frequency(NB_TIRAGES,t_mers_twist,32));
	printf("AES : %f\n",Frequency(NB_TIRAGES,t_aes,32));
	printf("Rand poids faible : %f\n",Frequency(NB_TIRAGES,t_rand_faible,4));
	printf("Rand poids fort : %f\n",Frequency(NB_TIRAGES,t_rand_fort,4));
	
	fclose(f_von_newman);
	fclose(f_mers_twist);
	fclose(f_aes);
	fclose(f_rand_faible);
	fclose(f_rand_fort);
	
	return 0;
}
