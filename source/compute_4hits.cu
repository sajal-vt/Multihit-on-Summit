/*
    <CombinationFinder identifies 4-hit combinations of carcinogenic genes.>
    Copyright (C) <2020>  <CombinationFinderDevs>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

/***********************************************************************
 * fComb.c
 *
 * Identify multi-hit combinations based on LR-
  ofstream writer;
           stringstream fname, gene1, gene2, gene3, gene4, fvalue;
           char buffer[100];
           snprintf(buffer, sizeof(buffer), "fvalues-rank-%d.txt", world_rank);
           writer.open(buffer);
           for(int i = 0; i < two_hit_genes_post; i++) {
               gene1 << comb[i].gene1;
               gene2 << comb[i].gene2;
               gene3 << comb[i].gene3;
               gene4 << comb[i].gene4;
               fvalue << comb[i].f_max;
               writer << gene1.str() + " " + gene2.str() + " " + gene3.str() + " " + gene4.str() + " : " + fvalue.str() + "\n";
           }
           writer.close();
           if(world_rank == 0) {
               cout << "device id: " << dev_i << " and maxF: " << comb[0].f_max << endl;
           }
 *
 * Calling Parameters: Tumor matrix file
 *                     Normal gene-sample list file
 * 
 * Output: list of 3-hit combinations
 ************************************************************************/

#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <bitset>
#include <mpi.h>
#include <fstream>
#include <sstream>

static const int NAME_LEN  = 20;
#define NUM_HITS 5
#define NUM_BITS 64
#define gpuErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }

int current_num_samples = 0;
int current_num_group_samples = 0;
unsigned long long int* tumor_matrix;
using namespace std;

inline void gpuAssert(cudaError_t code, const char * file, int line, bool abort=true) {
   if( code != cudaSuccess ) {
      fprintf(stderr, "GPU Assert: %s Error ID: %d %s %d \n", cudaGetErrorString(code), code,file, line);
      if( abort ) exit (code);
   }
}

typedef struct {
   float f_max;
   int   gene1;
   int   gene2;
   int   gene3;
   int   gene4;
} multi_hit_comb; 


int count_high_bits_64( unsigned long long int num )
{
   int count = 0;
   for ( count; num; count++ ) 
   {
      num &= num - 1;
   }
   return count;
}
/***********************************************************************
 *
 * get count of unique genes and samples from header of input file
 *
 ************************************************************************/
void getNumGenesSamplesTumor( FILE *fp_gene_sample_matrix, int *num_genes, int *num_samples )
{
   int     i, j, ret_value;
   char    *line = NULL;
   size_t  len = 0;
   ssize_t read;

   /* First line contains number of genes and samples */
   read = getline( &line, &len, fp_gene_sample_matrix );
   ret_value = sscanf( line, "%d %d", &i, &j );
   if (ret_value == 2 )
   {
      *num_genes   = i;
      *num_samples = j;
   }
   else
   {
      printf("ERROR: invalid input file header %d\n", ret_value);
      exit( 1 );
   }
}


/***********************************************************************
 * Load gene-sample matrix data from input file
 *
 * Calling Parameters: gene-sample matrix input file
 *                     number of genes
 *                     number of samples
 *                     gene_sample matrix (updated)
 *
 ************************************************************************/
void loadGeneSampleMatrixTumor( FILE *fp_gene_sample_matrix, int num_genes, 
      int num_samples, int num_group_samples, unsigned long long int *gene_sample_matrix, 
      unsigned long long int *gene_bit_matrix, char *gene_id, int *tumor_samples_per_gene )
{
   int     i, j, k, n, ret_value;
   char    *line = NULL;
   char    *gene, *sample;
   size_t  len = 0;
   ssize_t read;

   gene   = (char *)malloc( NAME_LEN * sizeof( char ) );
   sample = (char *)malloc( NAME_LEN * sizeof( char ) );

   /* initialize matrix */
   for ( n = 0; n < num_genes; n++ )
   {
      tumor_samples_per_gene[n] = 0;
      for ( j = 0; j < num_group_samples; j++ )
      {
	 gene_sample_matrix[n * num_group_samples + j] = 0;
      }
      for ( j = 0; j < num_samples; j++ )
      {
	 gene_bit_matrix[n * num_samples + j]    = 0;
      }
   }

   while ( !feof( fp_gene_sample_matrix ) )
   {
      read = getline( &line, &len, fp_gene_sample_matrix );
      ret_value = sscanf( line, "%d %d %d %s %s", &i, &j, &k, gene, sample );
      if ( ret_value == 5 )
      {
	 strcpy( gene_id+(i*NAME_LEN), gene );
	 gene_bit_matrix[i * num_samples + j] = ( k > 0  ? 1 : 0 );
	 //gene_sample_matrix[i * num_samples + j] = k;
	 if ( k > 0 )
	 {
	    tumor_samples_per_gene[i]++;
	 }
      }
      else
      {
	 printf("ERROR: reading data from input file %d\n", ret_value);
	 exit( 1 );
      }
   }
   for ( n = 0; n < num_genes; n++ )
   {
      for ( j = 0; j < num_group_samples; j++ )
      {
	 for ( k = 0; k < NUM_BITS; k++ )
	 {
	    if ( num_samples > j * NUM_BITS + k ) 
	    {
	       gene_sample_matrix[n * num_group_samples + j] |=
		  ( (gene_bit_matrix[n * num_samples + ( j * NUM_BITS ) + k]) << (NUM_BITS - 1 - k) );
	    }
	 }
      }
   }      

   return;

}


/***********************************************************************
 *
 * get count of samples in normal gene-sample file
 *
 ************************************************************************/
int getNumSamplesNormal( FILE *fp_gene_sample_list )
{
   int     num_samples, last_sample;
   int     sample, ret_value;
   char    gene[NAME_LEN];
   char    *line = NULL;
   size_t  len = 0;
   ssize_t read;
   num_samples = last_sample = 0;

   /* read to end of file */
   while ( !feof( fp_gene_sample_list ) )
   {
      read = getline( &line, &len, fp_gene_sample_list );
      ret_value = sscanf( line, "%s %d", gene, &sample );
      if (ret_value == 2 )
      {
	 if ( sample != last_sample )
	 {
	    num_samples++;
	    last_sample = sample;
	 }
      }
      else
      {
	 printf("ERROR: invalid line in normal gene-sample list %d\n", ret_value);
	 exit( 1 );
      }
   }
   rewind( fp_gene_sample_list );

   return( num_samples );
}


/***********************************************************************
 * Load gene-sample matrix data from input file
 *
 * Calling Parameters: gene-sample matrix input file
 *                     number of genes
 *                     number of samples
 *                     gene_sample matrix (updated)
 *
 ************************************************************************/
void loadGeneSampleMatrixNormal( FILE *fp_gene_sample_list, int num_genes, 
      int num_samples_normal, int num_group_samples_normal, unsigned long long int *normal_matrix,
      unsigned long long int *normal_bit_matrix, char *gene_id, int *normal_samples_per_gene )
{
   int     j, n, k, ret_value;
   char    *line = NULL;
   char    gene[NAME_LEN];
   int     sample, new_sample_id, matrix_sample_index;
   size_t  len = 0;
   ssize_t read;
   int     *sample_id; /* to translate from sample# in file to sequential # */
   int     max_samples = 1000; /* for allocation of sample_id list */

   /* list of old to new (sequential excluding missing numbers) sample ids */
   sample_id = (int *)malloc( max_samples * sizeof( int ) );
   if ( sample_id == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal gene_sample_matrix \n" );
      exit( 1 );
   }
   for ( j = 0; j < max_samples; j++ )
   {
      sample_id[j] = -1;
   }

   /* initialize matrix */
   for ( n = 0; n < num_genes; n++ )
   {
      normal_samples_per_gene[n] = 0;
      for ( j = 0; j < num_samples_normal; j++ )
      {
	 normal_bit_matrix[n * num_samples_normal + j] = 0;
      }
      for ( j = 0; j < num_group_samples_normal; j++ )
      {
	 normal_matrix[n * num_group_samples_normal + j] = 0;
      }
   }

   new_sample_id = 0;
   while ( !feof( fp_gene_sample_list ) )
   {
      read = getline( &line, &len, fp_gene_sample_list );
      ret_value = sscanf( line, "%s %d", gene, &sample );
      if ( ret_value == 2 )
      {
	 if ( sample_id[sample] < 0 )
	 {
	    sample_id[sample] = new_sample_id;
	    matrix_sample_index = new_sample_id;
	    new_sample_id++;
	 }
	 else
	 {
	    matrix_sample_index = sample_id[sample];
	 }
	 for ( n = 0; n < num_genes; n++ )
	 {
	    if ( strcmp( gene_id+(n*NAME_LEN), gene ) == 0 )
	    {
	       normal_bit_matrix[n * num_samples_normal + matrix_sample_index] = 1;
	       normal_samples_per_gene = 0;
	       break;
	    }
	 }
      }
      else
      {
	 printf("ERROR: reading data from normal gene-sample input file %d\n", ret_value);
	 exit( 1 );
      }
   }
   for ( n = 0; n < num_genes; n++ )
   {
      for ( j = 0; j < num_group_samples_normal; j++ ) 
      {
	 for ( k = 0; k < NUM_BITS; k++ )
	 {
	    if ( num_samples_normal > (j * NUM_BITS + k) )
	    {
	       normal_matrix[n * num_group_samples_normal + j] |= 
		  ( (normal_bit_matrix[n * num_samples_normal + ( j * NUM_BITS ) + k]) << (NUM_BITS - 1 - k) );
	    }
	 }
      }
   }
   free( sample_id );
}


/***********************************************************************
 * Find combination with smallest lr- 
 *
 * Calling Parameters: tumor gene-sample matrix
 *                     number of genes
 *                     number of samples
 *                     normal gene-sample matrix
 *                     tumor genes_per_sample for lr bound
 * 
 ************************************************************************/
__global__ void maxF( unsigned long long int *tumor_matrix, int num_genes, int two_hit_genes,
      int num_samples_tumor, int num_group_samples_tumor, 
      unsigned long long int *normal_matrix, int num_samples_normal,
      int num_group_samples_normal, int *tumor_samples_per_gene, 
      int *normal_samples_per_gene, unsigned long long int *excluded_samples,
      multi_hit_comb *comb, float beta, int dev_i, int numThreads, int start, int end)
{
   int   i1, i2, i3, i4, j, count;
   int k;
   int   true_pos, false_pos, true_neg, false_neg;
   float f;          /* f-measure */
   float f_bound;
   int   num_skipped;
   int   temp_tp;
   int localIdx = threadIdx.x;
   unsigned long long int bit_true_pos, bit_false_pos;
    
   num_skipped = 0;
   k  = blockDim.x * blockIdx.x + threadIdx.x;
   k  = k + start; // Adding offset: start is this GPU's start point, not the node's start point
   if (k >= two_hit_genes || k > end) return;
   i2 = (int) (sqrt(0.25 + (2 * k)) - 0.5);
   i1 = k - (i2 * (i2 + 1) * 0.5);
   
   if (i1 >= num_genes || i2 >= num_genes || i1 >= i2 ) return;
   comb[k].f_max = 0.0;
  

   for( i3 = i2 + 1; i3 < num_genes; i3++) {
      //true_pos = 0;
      for(i4 = i3 + 1; i4 < num_genes; i4++) {
              true_pos = 0;
	      for ( j = 0; j < num_group_samples_tumor; j++ )
	      {
	         bit_true_pos = tumor_matrix[i1 *num_group_samples_tumor + j]
		    & tumor_matrix[i2 * num_group_samples_tumor + j]
		    & tumor_matrix[i3 * num_group_samples_tumor + j]
                    & tumor_matrix[i4 * num_group_samples_tumor + j];

		 true_pos += __popcll(bit_true_pos);

	      }
	      false_pos = 0;
	      for ( j = 0; j < num_group_samples_normal; j++ )
	      {
		 bit_false_pos = normal_matrix[i1 * num_group_samples_normal + j]
		    & normal_matrix[i2 * num_group_samples_normal + j]
		    & normal_matrix[i3 * num_group_samples_normal + j]
                    & normal_matrix[i4 * num_group_samples_normal + j];
		 false_pos += __popcll(bit_false_pos);
	      }
	      if ( true_pos > 0 ) /* avoid divide by zero */
	      {
		 false_neg = num_samples_tumor  - true_pos;
		 true_neg  = num_samples_normal - false_pos; 
		 f         = (float)(beta * (float) true_pos + (float) true_neg) / (float)(num_samples_tumor + num_samples_normal);
		 if(f > comb[k].f_max) {
		    comb[k].f_max = f;
		    comb[k].gene1 = i1;
		    comb[k].gene2 = i2;
		    comb[k].gene3 = i3;
                    comb[k].gene4 = i4;
		 }
	      }
	   }
      }
}
/***********************************************************************
 * exclude samples contatining gene combinations found
 *
 * Calling Parameters: gene index
 *                     excluded samples list (for update)
 *                     tumor gene-sample matrix
 *                     number of samples
 * 
 ************************************************************************/
int excludeSamples( multi_hit_comb c, unsigned long long int *excluded_samples, unsigned long long int *tumor_matrix, int num_group_samples_tumor ) 
{
   int   s, num_excluded;
   unsigned long long int num_excluded_bits;
   num_excluded = 0;
   for ( s = 0; s < num_group_samples_tumor; s++ )
   {
      num_excluded_bits   = ( ~excluded_samples[s] ) 
	 & tumor_matrix[c.gene1 * num_group_samples_tumor + s]
	 & tumor_matrix[c.gene2 * num_group_samples_tumor + s] 
	 & tumor_matrix[c.gene3 * num_group_samples_tumor + s];
      num_excluded        += count_high_bits_64( num_excluded_bits );
      excluded_samples[s] |= tumor_matrix[c.gene1 * num_group_samples_tumor + s] 
	 & tumor_matrix[c.gene2 * num_group_samples_tumor + s] 
	 & tumor_matrix[c.gene3 * num_group_samples_tumor + s]; 
   }
   return( num_excluded );
}

/***********************************************************************
 * exclude samples contatining gene combinations found
 *
 * Calling Parameters: gene index
 *                     excluded samples list (for update)
 *                     tumor gene-sample matrix
 *                     number of samples
 * 
 ************************************************************************/
int excludeSamples2( int num_genes, multi_hit_comb c) 
{  
   //clock_t begin, end;
   //begin = clock();
   int   s, num_excluded;
   unsigned long long int* num_excluded_bits = new unsigned long long int[current_num_group_samples];
   num_excluded = 0;
   
   for ( s = 0; s < current_num_group_samples; s++ )
   {
      num_excluded_bits[s]   =  
	 tumor_matrix[c.gene1 * current_num_group_samples + s]
	 & tumor_matrix[c.gene2 * current_num_group_samples + s] 
	 & tumor_matrix[c.gene3 * current_num_group_samples + s]
         & tumor_matrix[c.gene4 * current_num_group_samples + s];
      num_excluded        += count_high_bits_64( num_excluded_bits[s] );
   }
  
   //cout << "excludedSamples2: number of excluded samples: " << num_excluded << endl; 
   int num_samples = current_num_samples - num_excluded;
   int num_group_samples = num_samples / NUM_BITS + 1;

   unsigned short* regular_matrix = new unsigned short[num_samples * num_genes];
   for(int i = 0; i < num_samples * num_genes; i++)regular_matrix[i]=0; 
     
   // Remove the columns with covered samples
   int fcount = -1;
   unsigned long long int unit = 1;
   for(int s = 0; s < current_num_group_samples; s++) {
       unsigned long long int mask = 1;
       mask = mask << 63;
          
       unsigned long long int excluded_bits = num_excluded_bits[s];
       
       for(int b = 0; b < NUM_BITS; b++) {
            
           if(!(mask & excluded_bits)) {       
               fcount++;
               if(fcount < num_samples)
               for(int g = 0; g < num_genes; g++) {
                   //if(num_samples * g + s * NUM_BITS + fcount >= num_samples * num_genes)break;
                   regular_matrix[num_samples * g + fcount] = 
                    (tumor_matrix[current_num_group_samples * g + s] & (unit << (NUM_BITS - b - 1))) > 0 ? 1 : 0;   
                }
            }  
            mask = mask >> 1;

       }
   }
   current_num_samples = num_samples;
   current_num_group_samples = num_group_samples;
   tumor_matrix = new unsigned long long int[num_genes * current_num_group_samples];
   //cout << "num samples, group samples: " << num_samples << " and " << num_group_samples << endl;
   for(int i = 0; i < num_genes * current_num_group_samples; i++)tumor_matrix[i] = 0;
   for(int i = 0; i < num_genes * current_num_samples; i++) {
      int row = i / current_num_samples;
      int col = i % current_num_samples;
      int comp_col = col / NUM_BITS;
      int bit_pos = col % NUM_BITS;
      unsigned long long int value = (unsigned long long int)regular_matrix[i];
      tumor_matrix[row * current_num_group_samples + comp_col] |= (value << (NUM_BITS - 1 - bit_pos));
   }
  
   delete[] num_excluded_bits;
   delete[] regular_matrix;
   //end = clock();
   //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
   //cout << "Elapsed time ins exclude samples (s): " << elapsed_secs << endl;
   return( num_excluded );
}
/*************************************************************************
 * parallel reduction to find the maximum f_max in O(logn) time
 * 
 * Calling Parameters: two hit combination structure input
 *                     two hit combination structure output
 *
 ************************************************************************/
__global__ void parallelReduceMax(multi_hit_comb *comb_in, multi_hit_comb *comb_out, int num_genes) {
   extern __shared__ multi_hit_comb shared_data[];
   unsigned int tid = threadIdx.x;
   unsigned int i   = blockIdx.x * blockDim.x + threadIdx.x;
   multi_hit_comb temp;
   temp.f_max = 0.0;
   temp.gene1 = 0;
   temp.gene2 = 0;
   temp.gene3 = 0;
   temp.gene4 = 0;
   if ( i < num_genes )
      temp = comb_in[i];

   shared_data[tid] = temp;
   __syncthreads();
   for ( unsigned int s = blockDim.x/2; s > 0; s >>= 1) {
      if (tid < s) {
	 if( shared_data[tid + s].f_max > shared_data[tid].f_max) {
	    shared_data[tid] = shared_data[tid + s];
	 }
      }
      __syncthreads();
   }

   if (tid == 0 ){
      comb_out[blockIdx.x] = shared_data[0];
   }
}

void print_combinations(int i, multi_hit_comb max_comb, char* gene_id, int num_excluded, int tot_excluded) {
   char *gene1_name, *gene2_name, *gene3_name, *gene4_name;
   gene1_name           = (char *)malloc( NAME_LEN          * sizeof( char ) );
   gene2_name           = (char *)malloc( NAME_LEN          * sizeof( char ) );
   gene3_name           = (char *)malloc( NAME_LEN          * sizeof( char ) );
   gene4_name           = (char *)malloc(NAME_LEN * sizeof(char));
   strcpy( gene1_name, gene_id+((max_comb.gene1)*NAME_LEN) );
   strcpy( gene2_name, gene_id+((max_comb.gene2)*NAME_LEN) );
   strcpy( gene3_name, gene_id+((max_comb.gene3)*NAME_LEN) );
   strcpy( gene4_name, gene_id+((max_comb.gene4)*NAME_LEN) );
   printf( "%d- %s %s %s %s %d %d %d %d F-max = %9.6f , num excluded %d, tot excluded %d \n", 
	 i, gene1_name, gene2_name, gene3_name, gene4_name, max_comb.gene1, max_comb.gene2, max_comb.gene3, max_comb.gene4,
	 max_comb.f_max, num_excluded, tot_excluded);
   free(gene1_name);
   free(gene2_name);
   free(gene3_name);
   free(gene4_name);
}


/***********************************************************************
 * list all combinations found in gene-sample matrix
 *
 * Calling Parameters: tumor gene-sample matrix 
 *                     number of genes
 *                     number of samples (tumor)
 *                     normal gene-sample matrix
 *                     number of samples (normal)
 *                     gene_id (name list)
 *                     tumor samples per gene 
 * 
 ************************************************************************/
int listCombs(int num_genes, int num_samples_tumor, 
      int num_group_samples_tumor, unsigned long long int *normal_matrix, int num_samples_normal, 
      int num_group_samples_normal, char *gene_id, int *tumor_samples_per_gene, 
      int *normal_samples_per_gene, float beta, int ngpus, int world_size, int world_rank, int **schedule)
{
   int   it, num_found, num_excluded, tot_excluded;
   int two_hit_genes =  ceil(num_genes *( num_genes + 1)/2);
   int threadsPerBlock_maxF;
   unsigned long long int  *excluded_samples, **excluded_samples_d;
   int threadsPerBlock, sharedMemSize;
   multi_hit_comb **comb, **comb_d, **comb_d_result, *ecomb;
   unsigned long long int **tumor_matrix_dev;
   unsigned long long int **normal_matrix_dev;
   int          **tumor_samples_per_gene_dev;
   int          **normal_samples_per_gene_dev;
   int numDev;
   cudaGetDeviceCount(&numDev);
   //cout << "Number of seen devices: " << numDev << endl;
   // Forcing to use fixed number of GPUs
   int visible_gpus[] = {0, 1, 2, 3, 4, 5};
  
   
   int total_gpus = ngpus * world_size;    
   numDev = ngpus;
   //int num_threads_per_gpu = (int)ceil(two_hit_genes / total_gpus); 
 
   threadsPerBlock_maxF = 512; 
   threadsPerBlock      = 128; // For parallel redux

   // These parameters need to be set per GPU
   int* num_threads_per_gpu = new int[ngpus];
   int* blocksPerGrid = new int[ngpus];
   int* blocksPerGrid_maxF = new int[ngpus];
   int* two_hit_genes_post = new int[ngpus];
   comb = new multi_hit_comb*[ngpus];
   for(int i = 0; i < ngpus; i++) {
       num_threads_per_gpu[i] = schedule[i][1] - schedule[i][0] + 1; // Inclusive indices
       blocksPerGrid[i]        = (num_threads_per_gpu[i] + threadsPerBlock - 1) / threadsPerBlock;
       blocksPerGrid_maxF[i]   = (num_threads_per_gpu[i] + threadsPerBlock_maxF - 1 ) / threadsPerBlock_maxF;
       two_hit_genes_post[i]   = blocksPerGrid[i];
       //cout << "THGP: " << two_hit_genes_post[i] << endl;
       comb[i] = (multi_hit_comb *)malloc( sizeof(multi_hit_comb) * two_hit_genes_post[i]); 
   }
   
   // These variables don't rely on workload of the GPU
   sharedMemSize        = threadsPerBlock * sizeof(multi_hit_comb);
   ecomb                 = (multi_hit_comb *)malloc( sizeof(multi_hit_comb) * two_hit_genes);
   excluded_samples     = (unsigned long long int  *)malloc( num_group_samples_tumor * sizeof( unsigned long long int ) );
   int dev_i;
   comb_d = new multi_hit_comb*[ngpus];
   comb_d_result = new multi_hit_comb*[ngpus];
   tumor_matrix_dev = new unsigned long long int*[ngpus];
   normal_matrix_dev = new unsigned long long int*[ngpus];
   tumor_samples_per_gene_dev = new int*[ngpus];
   normal_samples_per_gene_dev = new int*[ngpus];
   excluded_samples_d = new unsigned long long int*[ngpus];
  
   if(world_rank == 0) {
       cout << "Inspecting GPU block and thread size parameters:" << endl;
       cout << "threadsPerBlock_maxF: " << threadsPerBlock_maxF << endl;
       cout << "twoHitGenes: " << two_hit_genes << endl;
       cout << "blocksPerGrid: " << blocksPerGrid << endl;
       cout << "two_hit_genes_post:" << two_hit_genes_post << endl;
  
       cout << "NGST: " << num_group_samples_tumor << endl;
   }
   
   double begin, end;
   //begin = clock();
   for( int i = 0; i < ngpus; i++) {
      dev_i = visible_gpus[i];
      cudaSetDevice(dev_i);
      //cout << "Checking the error block with world_rank " << world_rank << " and " << dev_i << "and thgp: " << two_hit_genes_post[dev_i] << endl;
      gpuErrorCheck( cudaMalloc( &comb_d[dev_i] ,       sizeof( multi_hit_comb ) * two_hit_genes) ); 
      gpuErrorCheck( cudaMalloc( &comb_d_result[dev_i], sizeof( multi_hit_comb ) * two_hit_genes_post[dev_i] ) );
      gpuErrorCheck( cudaMalloc( &tumor_matrix_dev[dev_i],  num_genes * num_group_samples_tumor * sizeof( unsigned long long int ) ));
      gpuErrorCheck( cudaMalloc( &normal_matrix_dev[dev_i], num_genes * num_group_samples_normal* sizeof( unsigned long long int ) ));
      gpuErrorCheck( cudaMalloc( &tumor_samples_per_gene_dev[dev_i],  num_genes * sizeof( int ) ) );
      gpuErrorCheck( cudaMalloc( &normal_samples_per_gene_dev[dev_i], num_genes * sizeof( int ) ) );
      gpuErrorCheck( cudaMalloc( &excluded_samples_d[dev_i], num_group_samples_tumor * sizeof( unsigned long long int ) ) );
      gpuErrorCheck( cudaMemcpyAsync( tumor_matrix_dev[dev_i],  tumor_matrix,  num_genes * num_group_samples_tumor  * sizeof( unsigned long long int ), cudaMemcpyHostToDevice ) );
      gpuErrorCheck( cudaMemcpyAsync( normal_matrix_dev[dev_i], normal_matrix, num_genes * num_group_samples_normal * sizeof( unsigned long long int ), cudaMemcpyHostToDevice ) );
      gpuErrorCheck( cudaMemcpyAsync( tumor_samples_per_gene_dev[dev_i],  tumor_samples_per_gene,  num_genes * sizeof( int ), cudaMemcpyHostToDevice ) );
      gpuErrorCheck( cudaMemcpyAsync( normal_samples_per_gene_dev[dev_i], normal_samples_per_gene, num_genes * sizeof( int ), cudaMemcpyHostToDevice ) );
   }

   // Wait till all data transfer is complete
   for(dev_i = 0; dev_i < numDev; dev_i++) {
       cudaSetDevice(dev_i);
       cudaDeviceSynchronize();
   }
   //end = clock();
   //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
   //std::cout << "Time spent in copying data: " << elapsed_secs << std::endl;
   
   if ( excluded_samples == NULL )
   {
      printf( "ERROR: failed to allocate memory for excluded_samples \n" );
      exit( 1 );
   }
   for ( int i = 0; i < num_group_samples_tumor; i++ )
   {
      excluded_samples[i] = 0;
   }
   num_found    = 0;
   tot_excluded = 0;
   it = 1;
   //cout << "Inside main computation loop" << endl;
   //cudaStream_t *streams = new cudaStream_t[numDev];
   while ( tot_excluded < num_samples_tumor )
   {
      //cout << world_rank << ": Executing iteration # " << it << endl;
      
      for ( int i = 0; i < numDev; i++) {
         //begin = clock();
	 dev_i = visible_gpus[i];
         cudaSetDevice(dev_i);
         gpuErrorCheck( cudaMemcpyAsync( tumor_matrix_dev[dev_i],  tumor_matrix,  num_genes * current_num_group_samples  * sizeof( unsigned long long int ), cudaMemcpyHostToDevice ) );
         maxF<<< blocksPerGrid_maxF[i], threadsPerBlock_maxF >>>( tumor_matrix_dev[dev_i],
	       num_genes, two_hit_genes, current_num_samples,
	       current_num_group_samples, normal_matrix_dev[dev_i],
	       num_samples_normal, num_group_samples_normal, 
	       tumor_samples_per_gene_dev[dev_i], 
	       normal_samples_per_gene_dev[dev_i], 
	       excluded_samples_d[dev_i], comb_d[dev_i], beta, i, num_threads_per_gpu[i], schedule[i][0], schedule[i][1] );

	 gpuErrorCheck( cudaGetLastError() );
      }
      // Wait till all F-computation is complete
       for(int i = 0; i < numDev; i++) {
           dev_i = visible_gpus[i];
           cudaSetDevice(dev_i);
	   gpuErrorCheck( cudaDeviceSynchronize() );
       }
       // Perform parallel reduction
       for(int i = 0; i < ngpus; i++) {
           dev_i = visible_gpus[i]; 
           cudaSetDevice(dev_i);    
           //cout << "Performing parallel reduction" << endl;
     
	   /* Parallel Reduction called multiple time */
	   parallelReduceMax <<< blocksPerGrid[dev_i], threadsPerBlock, sharedMemSize>>>(comb_d[dev_i], comb_d_result[dev_i], two_hit_genes);
	   gpuErrorCheck(cudaGetLastError());
	   gpuErrorCheck(cudaDeviceSynchronize());
	   int newSize   = two_hit_genes;
	   while(newSize > 1) {
	       parallelReduceMax<<< blocksPerGrid[dev_i], threadsPerBlock, sharedMemSize>>>( comb_d_result[dev_i], comb_d_result[dev_i], two_hit_genes_post[dev_i] );
	       gpuErrorCheck(cudaGetLastError() );
	       gpuErrorCheck(cudaDeviceSynchronize());
	       newSize = (newSize + threadsPerBlock - 1)/ threadsPerBlock;
               //cout << "NewSize = " << newSize << endl;
	   }

	   parallelReduceMax<<< blocksPerGrid[dev_i], threadsPerBlock, sharedMemSize>>>( comb_d_result[dev_i], comb_d_result[dev_i], two_hit_genes_post[dev_i] );
	   gpuErrorCheck(cudaGetLastError() );
	   gpuErrorCheck(cudaDeviceSynchronize());
	   parallelReduceMax<<< blocksPerGrid[dev_i], threadsPerBlock, sharedMemSize>>>( comb_d_result[dev_i], comb_d_result[dev_i], two_hit_genes_post[dev_i] );
	   gpuErrorCheck(cudaGetLastError() );
       }

      
       // Wait till reduction is complete
       for(int i = 0; i < ngpus; i++) {
           //cudaSetDevice(visible_gpus[i]);
	   gpuErrorCheck(cudaDeviceSynchronize());
       }


      

       //cout << "Reading back from devices" << endl;
       multi_hit_comb max_comb, max_comb_rank;
       max_comb.f_max = -1;
       for(int i = 0; i < ngpus; i++) {
           dev_i = visible_gpus[i]; 
	   cudaSetDevice(dev_i);
           gpuErrorCheck( cudaMemcpy(comb[dev_i], comb_d_result[dev_i], sizeof(multi_hit_comb) * two_hit_genes_post[dev_i], cudaMemcpyDeviceToHost)); //
	   if(comb[dev_i][0].f_max >= max_comb.f_max) 
	   {
	       max_comb = comb[dev_i][0];
	   }
      }

      // Create multi_hit_comb datatype for MPI communication
      int blocklengths[NUM_HITS] = {1,1,1,1,1};
      MPI_Datatype mpi_types[NUM_HITS] = {MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
      MPI_Datatype mpi_multi_hit;
      MPI_Aint     disp[5];
      disp[0] = offsetof(multi_hit_comb, f_max);
      disp[1] = offsetof(multi_hit_comb, gene1);
      disp[2] = offsetof(multi_hit_comb, gene2);
      disp[3] = offsetof(multi_hit_comb, gene3);
      disp[4] = offsetof(multi_hit_comb, gene4);
      //cout << "MPI activities" << endl;

      MPI_Type_create_struct(NUM_HITS, blocklengths, disp, mpi_types, &mpi_multi_hit);
      MPI_Type_commit(&mpi_multi_hit); 
       // Default, main thread is 0. Receives the max from all the sources and finally finds the max
      if(world_rank == 0) {
         //cout << it << "Master worker gathering result" << endl;
         multi_hit_comb source_max;
         float true_fmax = -1;

         int src_i;
         for( src_i = 1; src_i < world_size; src_i++) {
            MPI_Status status;
            MPI_Recv(&source_max, 1, mpi_multi_hit, src_i, 15, MPI_COMM_WORLD, &status);
            if(source_max.f_max > max_comb.f_max) {
               max_comb = source_max;
            }
            //cout << "Received from rank " << src_i << endl;
         }
         int rcv_i;
         for(rcv_i = 1; rcv_i < world_size; rcv_i++) {
             MPI_Status status;
             MPI_Send(&max_comb, 1, mpi_multi_hit, rcv_i, 15, MPI_COMM_WORLD);
         }
      } 
      if(world_rank > 0) {
         multi_hit_comb global_max;
         MPI_Send(&max_comb, 1, mpi_multi_hit, 0, 15, MPI_COMM_WORLD);
         MPI_Status status;
         MPI_Recv(&max_comb, 1, mpi_multi_hit, 0, 15, MPI_COMM_WORLD, &status);
      }
       


      num_excluded = excludeSamples2( num_genes, max_comb); //, excluded_samples, tumor_matrix, num_group_samples_tumor );
      tot_excluded += num_excluded;
      num_found++;
      if(world_rank == 0) {

          print_combinations(it, max_comb, gene_id, num_excluded, tot_excluded);
      }
      
      //if (it == 1 )
      //{
      //    break;
      //}
      
      it++;
         //print_combinations(i, max_comb, gene_id);
 
      //cout << it << ": " << num_excluded << ", " << "tot_excluded" << endl;
      if(num_excluded == 0) {break;}
   }
   //free(comb);
   //free(excluded_samples);
   
   for(int i = 0; i < ngpus; i++) {
      dev_i = visible_gpus[i];
      cudaSetDevice(dev_i);
      cudaFree(comb_d);
      cudaFree(comb_d_result);
      cudaFree(tumor_matrix_dev);
      cudaFree(normal_matrix_dev);
      cudaFree(excluded_samples_d);
      cudaFree(tumor_samples_per_gene_dev);
      cudaFree(normal_samples_per_gene_dev);
      delete comb[dev_i];

   }
   delete num_threads_per_gpu;
   delete blocksPerGrid;
   delete blocksPerGrid_maxF;
   delete two_hit_genes_post;
   delete comb;
   return( num_found );
}

float computeFvalue(char* gene_id, unsigned long long int * tumor_matrix, unsigned long long int * normal_matrix, int g1, int g2, int g3, int g4, int num_genes, int num_tumor_samples, int num_normal_samples, float beta) {
    int covered_tumor_samples = 0;
    int covered_normal_samples = 0;

    for(int i = 0; i < num_tumor_samples; i++) {
        if(tumor_matrix[g1 * num_tumor_samples + i] == 1 && tumor_matrix[g2 * num_tumor_samples + i] == 1 && 
        tumor_matrix[g3 * num_tumor_samples + i] == 1 && tumor_matrix[g4 * num_tumor_samples + i] == 1) {
            covered_tumor_samples++;
        }
    }
    for(int i = 0; i < num_normal_samples; i++) {
        if(normal_matrix[g1 * num_normal_samples + i] == 1 && normal_matrix[g2 * num_normal_samples + i] == 1 && 
        normal_matrix[g3 * num_normal_samples + i] == 1 && normal_matrix[g4 * num_normal_samples + i] == 1) {
            covered_normal_samples++;
        }
    }

    cout << "Covered tumor samples: " << covered_tumor_samples << endl;
    cout << "Covered normal samples: " << covered_normal_samples << endl;
    int TP = covered_tumor_samples;
    int FP = covered_normal_samples;
    int FN = num_tumor_samples - TP;
    int TN = num_normal_samples - FP;

    float f = (float)(beta * (float) TP + (float) TN) / (float)(num_tumor_samples + num_normal_samples);
    cout << "Fvalue for: " << g1 << ", " << g2 << ", " << g3 << ", " << g4 << " : " << f << endl;
    char *gene1_name, *gene2_name, *gene3_name, *gene4_name;
    gene1_name           = (char *)malloc( NAME_LEN          * sizeof( char ) );
    gene2_name           = (char *)malloc( NAME_LEN          * sizeof( char ) );
    gene3_name           = (char *)malloc( NAME_LEN          * sizeof( char ) );
    gene4_name           = (char *)malloc(NAME_LEN * sizeof(char));
    strcpy( gene1_name, gene_id+(g1*NAME_LEN) );
    strcpy( gene2_name, gene_id+(g2*NAME_LEN) );
    strcpy( gene3_name, gene_id+(g3*NAME_LEN) );
    strcpy( gene4_name, gene_id+(g4*NAME_LEN) );
    
    cout << gene1_name << " " << gene2_name << " " << gene3_name << " " << gene4_name << endl;
 
}


void read_schedule(char* fileName, int** schedule) {
    ifstream inFile;
    inFile.open(fileName);
    

    int start, end;
    int i = 0;
    while(inFile >> start >> end) {
        schedule[i][0] = start;
        schedule[i][1] = end;
        i++;
    }

    inFile.close(); 
}



/*****************************************************************************************
 * main()                                                                                *
 *                                                                                       *
 * calculate ratio of occurunce os multi-hit gene combinations                           *
 *                                                                                       *
 * Calling Parameters: list of tumor and normal gene-sample counts,                      *
 *                     list of freq mutated genes                                        *
 *                                                                                       *
 * output: ratio of occurance of multi-hit gene combinations in normal and tumor samples *
 *****************************************************************************************/
int main(int argc, char ** argv)
{
   clock_t begin, end;
   //cout << "Inside main" << endl;
  
   int            		num_genes, num_samples, num_samples_normal;   /* number of genes and samples in tcga maf data */
   int            		num_group_samples, num_group_samples_normal;
   int            		num_comb;                 /* number of 3-hit combinations found in tcga samples */
   //unsigned long long int       *tumor_matrix;            /* matrix of genesxsamples */
   unsigned long long int       *normal_matrix;           /* matrix of genesxsamples */
   unsigned long long int           *tumor_bit_matrix;
   unsigned long long int           *normal_bit_matrix;
   int            		*tumor_samples_per_gene;  /* total True Pos for calulating bound */
   int            		*normal_samples_per_gene;  /* total False Pos for calulating bound */
   FILE           		*fp_tumor_matrix;
   FILE           		*fp_normal_matrix;
   char           		*gene_id;                 /* list of gene ids */
   float          		beta;
   int ngpus;
 
   //cout << "Attempting  MPI_Init" << endl;  
   MPI_Init(NULL, NULL);
   //cout << "Performed MPI_Init" << endl;
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

 
   if (argc != 7)
   {
      printf("ERROR: fComb requires 6 parameters\n");
      printf("       - Tumor gene-sample count matrix file (training)\n");
      printf("       - Normal gene-sample list file (training)\n");
      printf("       - Beta value (Fbeta-score)\n");
      printf("       - GPUs per rank\n");
      printf("       - Schedule\n");
      printf("       - Number of nodes\n");
      exit(1);
   }

   if ( ( fp_tumor_matrix = fopen( argv[1], "r" ) ) == NULL )
   {
      printf( "ERROR: unable to open tumor gene-sample count matrix file %s, \n", argv[2] );
      exit( 1 );
   }

   if ( ( fp_normal_matrix = fopen( argv[2], "r" ) ) == NULL )
   {
      printf( "ERROR: unable to open normal gene-sample count matrix file %s, \n", argv[3] );
      exit( 1 );
   }

   beta = atof( argv[3] );
   ngpus = atoi( argv[4] );
   //cout << "read argv[4]" << endl;

  
   char* schedule_file = argv[5];
   //cout << schedule_file << endl;
   int num_nodes = atoi(argv[6]);
   //cout << num_nodes << endl;

   //cout << num_nodes << " " << schedule_file << endl;
   int ** schedule = new int*[num_nodes * ngpus];
   for(int i = 0; i < num_nodes * ngpus; i++) {
        schedule[i] = new int[2];
   }
    
   read_schedule(schedule_file, schedule); 

   //cout << "read schedule" << endl;
   if(world_rank == 0) begin = clock();
   /* load tumor gene-sample matrix*/
   getNumGenesSamplesTumor( fp_tumor_matrix, &num_genes, &num_samples );
   //printf( "Num Tumor genes = %d tumor samples = %d \n", num_genes, num_samples );
   num_group_samples = ( num_samples / NUM_BITS ) + 1;
   
   current_num_samples = num_samples;
   current_num_group_samples = num_group_samples;
   
   tumor_matrix = (unsigned long long int  *)malloc( num_genes * num_group_samples * sizeof( unsigned long long int ) );
   if ( tumor_matrix == NULL )
   {
      printf( "ERROR: failed to allocate memory for tumor gene_sample_matrix \n" );
      exit( 1 );
   }
   tumor_bit_matrix = (unsigned long long int *)malloc( num_genes * num_samples * sizeof( unsigned long long int ) );
   if ( tumor_bit_matrix  == NULL )
   {
      printf( "ERROR: failed to allocate memory for tumor bit matrix \n" );
      exit( 1 );
   }

   gene_id = (char  *)malloc( num_genes * NAME_LEN * sizeof( char ));
   if ( gene_id == NULL )
   {
      printf( "ERROR: failed to allocate memory for gene_ids \n" );
      exit( 1 );
   }
   tumor_samples_per_gene = (int  *)malloc( num_genes * sizeof( int ) );
   if ( tumor_samples_per_gene == NULL )
   {
      printf( "ERROR: failed to allocate memory for tumor samples per gene \n" );
      exit( 1 );
   }
   loadGeneSampleMatrixTumor( fp_tumor_matrix, num_genes, num_samples, 
	 num_group_samples, tumor_matrix, tumor_bit_matrix, gene_id, tumor_samples_per_gene );

   fclose( fp_tumor_matrix );

   /* load normal gene-sample matrix */
   num_samples_normal = getNumSamplesNormal( fp_normal_matrix );
   //printf( "Num normal samples = %d \n", num_samples_normal );
   num_group_samples_normal = ( num_samples_normal / NUM_BITS ) + 1;

   if(world_rank == 0) {
       cout << "Num genes: " << num_genes << " and num of tumor and normal samples: " << num_samples << ", " << num_samples_normal << endl;
   }
   normal_matrix = ( unsigned long long int *)malloc( num_genes * num_group_samples_normal * sizeof( unsigned long long int ) );
   if ( normal_matrix == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal gene_sample_matrix \n" );
      exit( 1 );
   }
   normal_bit_matrix = (unsigned long long int *)malloc( num_genes * num_samples_normal * sizeof( unsigned long long int ) );
   if ( normal_bit_matrix == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal bit matrix \n");
      exit( 1 );
   }

   normal_samples_per_gene = (int  *)malloc( num_genes * sizeof( int ) );
   if ( normal_samples_per_gene == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal samples per gene \n" );
      exit( 1 );
   }
   loadGeneSampleMatrixNormal( fp_normal_matrix, num_genes, num_samples_normal,
	 num_group_samples_normal, normal_matrix, normal_bit_matrix, gene_id, normal_samples_per_gene );

   fclose( fp_normal_matrix );
   double elapsed_secs;
   if(world_rank == 0) {
       end = clock();
       elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
       cout << "Time spent in reading files (s): " << elapsed_secs << endl;
       begin = clock();
   } 
   //cout << "My rank: " << world_rank << "/" << world_size << endl; 

   /* Check combinations of genes for coverage of samples */
   int** node_schedule = new int*[ngpus];
   for(int i = 0; i < ngpus; i++) {
       node_schedule[i] = new int[2];
       node_schedule[i][0] = schedule[world_rank * ngpus + i][0];
       node_schedule[i][1] = schedule[world_rank * ngpus + i][1];
   }
   num_comb = listCombs( num_genes, num_samples, num_group_samples, 
    	 normal_matrix, num_samples_normal, num_group_samples_normal, gene_id, 
   	 tumor_samples_per_gene, normal_samples_per_gene, beta, ngpus, world_size, world_rank, node_schedule);
   for(int i = 0; i < ngpus; i++) {
       delete node_schedule[i];
   }
   delete node_schedule;

   if(world_rank == 0) {
        printf( "Num 2-hit combinations = %d  (beta = %f )\n", num_comb, beta );

        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "Elapsed time in computation(s): " << elapsed_secs << endl;
    
   }
   // If we are the main process OR PROCESS 0
   // Then we receive everything from all the other processors and find the maximum value between the nodes
   MPI_Finalize(); 
   free( tumor_matrix );
   free( normal_matrix );
   free( gene_id );
   free( tumor_samples_per_gene );
   free( normal_samples_per_gene );
   for(int i = 0; i < num_nodes; i++) {
      delete schedule[i];
   }
   delete schedule;
    
   return( 0 );
}

