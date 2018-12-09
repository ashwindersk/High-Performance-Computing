#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"



// Define output file name
#define OUTPUT_FILE "stencil.pgm"



void stencil(const int nx, const int ny, float * image, float * tmp_image, int rank);
void init_image(const int nx, const int ny, float * image, float * tmp_image);
void output_image(const char * file_name, const int nx, const int ny, float *image);
double wtime(void);
int main(int argc, char *argv[]) {

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments

  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  // Allocate the images
  float *image =_mm_malloc(sizeof(float)*ny*nx,64);

  float *tmp_image = _mm_malloc(sizeof(float)*ny*nx,64);

  // Set the input image
  init_image(nx, ny, image, tmp_image);



  //Figuring out which processors are involved in the computation 
  MPI_Init(&argc, &argv);

  int size;
  int rank;

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int start = ((rank/16)*ny * ny);
  int end = ny* (((rank+1)/16 )*ny -1) + nx-1;

  int numElements = (end - start + 1);
  int numBytes = sizeof(float) * numElements;

  float *slice = malloc(numBytes);
  float *sliceTemp = malloc(numBytes);

  // memcpy(slice, image + start, numBytes);
  // memcpy(sliceTemp, tmp_image + start, numBytes);



  // Call the stencil kernel
  double tic = wtime();
  for (int t = 0; t < niters; ++t) {
    stencil(nx, ny, image, tmp_image, rank);
    stencil(nx, ny, tmp_image, image, rank);
  }
  double toc = wtime();


  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc-tic);
  printf("------------------------------------\n");

  output_image(OUTPUT_FILE, nx, ny, image);
  free(image);
}

void stencil(const int nx, const int ny,  float *restrict image, float *restrict tmp_image, int rank) {

  int size = ny*nx/16;
  float* buffer = malloc(size * sizeof(float));

  if(rank == 0 ){
    MPI_Scatter(image, size, MPI_FLOAT, buffer, size,MPI_FLOAT, 0,MPI_COMM_WORLD);    
  }
  // else if( rank >0 && rank < 15){
  
  
  if( rank < 16){

    printf("Current Rank is %d", rank);
    for( int i =0 ; i< size ; i++){
      //printf("value is %f", buffer[i]);
    }

  }
    // int start = 0 ;
    // int end = nx-1;
    // int numElements = end-start + 1;
    // int numBytes = sizeof(float) * numElements;

    // float *firstRow = malloc(numBytes);
    // memcpy(firstRow, image + start, numBytes);


    // start = (ny-1) * ny ;
    // end = (ny-1) * ny + nx-1;
    // numElements = end-start + 1;
    // numBytes = sizeof(float) * numElements;

    // float *lastRow = malloc(numBytes);
    // memcpy(lastRow, image + start, numBytes);

    // MPI_Send(firstRow, len(firstRow), MPI_FLOAT, rank-1, MPI_COMM_WORLD );

    // MPI_Send(lastRow, len(lastRow), MPI_FLOAT, rank+1, MPI_COMM_WORLD );

    // MPI_Recv()

  





  //manually amending the values of the corners
 tmp_image[0]                   = 0.6f * image[0]                  + 0.1f*image[1 + ny*0]                  + 0.1f*image[0 + ny*1];
 tmp_image[nx-1 + ny*0]         = 0.6f * image[nx-1 + ny*0]        + 0.1f*image[nx-2 + ny*0]               + 0.1f*image[nx-1 + ny*1];
 tmp_image[0 + ny*(nx-1)]       = 0.6f * image[0 + ny*(ny-1)]      + 0.1f*image[0 +ny*(ny-2)]              + 0.1f*image[(1 + ny*(ny-1))];
 tmp_image[nx-1 + (ny)*(ny-1)]  = 0.6f * image[nx-1 + (ny)*(ny-1)] + 0.1f*image[nx-1 + (ny)*(nx-2)]        + 0.1f*image[nx-2 +(nx-1)*(ny)];


  //top row
  for(int j = 1; j<nx-1; ++j){
    tmp_image[j+ny*0] = 0.1f*image[j-1 + ny*0] + 0.6f*image[j+ny*0]  + 0.1f*image[j+1 + ny*0] + 0.1f*image[j+ny*1];
  }

  //first column
  for(int i = 1; i< ny-1 ; ++i){
   tmp_image[0+ny*i] = 0.6f*image[0+ny*i] + 0.1f*image[1+ ny*i] + 0.1f*image[0+ny*(i-1)] + 0.1f*image[0 + ny*(i+1)];
  }

  //editing the values of the (ny-1)*(nx-1) pisxels
  for(int i = 1 ; i<ny-1; ++i){
   for(int j = 1 ; j<nx-1; ++j){
     int base = j+ny*i;
     __assume_aligned(image,64);
     __assume_aligned(tmp_image,64);
     #pragma omp simd
     
     tmp_image[base] = image[base-1]*0.1f   + image[base]*0.6f + image[base+1]*0.1f + image[base -ny]*0.1f + image[base +ny]*0.1f;
   }
  }
  //last column
  for(int i = 1; i< ny-1 ; ++i){
    int base  = nx-1 + ny*i;
    tmp_image[base] = 0.6f*image[base] + 0.1f*image[base-1] + 0.1f*image[base -ny] + 0.1f*image[base + ny];
  }


  //last row
  for(int j = 1; j<nx-1; ++j){
   tmp_image[j + ny*(nx-1)] = 0.6f*image[j+ ny*(nx-1)] + 0.1f*image[(j-1)+ ny*(nx-1)] + 0.1f*image[(j+1)+ ny*(nx-1)] + 0.1f*image[j+ ny*(nx-2)];
  }

 }

// Create the input image
void init_image(const int nx, const int ny, float * image, float * tmp_image) {
  // Zero everything
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
     
     image[j+ny*i] = 0.0;
     tmp_image[j+ny*i] = 0.0;
    }
  }

  // Checkerboard
   
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int ii = i*ny/8; ii < (i+1)*ny/8; ++ii) {
        for (int jj = j*nx/8; jj < (j+1)*nx/8; ++jj) {
          if ((i+j)%2)
          image[jj+ii*ny] = 100.0;

        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, const int nx, const int ny, float *image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  double maximum = 0.0;
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
     if (image[j+i*ny] > maximum)
       maximum = image[j+i*ny];

    }
  }

  // Output image, converting to numbers 0-255
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      //fputc((char)(255.0*image[j+ny*i]/maximum), fp);
     	fputc((char)(255.0*image[j+ny*i]/maximum),fp);
    }
  }

  // Close the file
  fclose(fp);

}

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
