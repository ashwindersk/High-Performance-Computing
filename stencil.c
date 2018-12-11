#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"



// Define output file name
#define OUTPUT_FILE "stencil.pgm"



void stencil(const int nx, const int ny, float * image, float * tmp_image,int rank);
void init_image(const int nx, const int ny, float * image, float * tmp_image);
void output_image(const char * file_name, const int nx, const int ny, float *image);
double wtime(void);
float *extractElements(float *subArray, float *array, int start, int end);
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
  float *image;

  float *tmp_image;

  // Set the input image
  

  
  //Figuring out which processors are involved in the computation 
  MPI_Init(&argc, &argv);

  int size;
  int rank;

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

  int sectionSize = ny*nx/16;
  float * bufferImg = (float *)malloc((ny*nx/16) * sizeof(float));
  float * bufferTempImg = (float *)malloc((ny*nx/16) * sizeof(float));


 
  if(rank ==0){
    image = malloc(sizeof(float)*ny*nx);

    tmp_image = malloc(sizeof(float)*ny*nx);

    init_image(nx,ny,image,tmp_image);

  }
  
  MPI_Scatter(image, sectionSize, MPI_FLOAT, bufferImg, sectionSize, MPI_FLOAT,0,MPI_COMM_WORLD);


  


  

  ny = ny/16;

  // Call the stencil kernel
  double tic = wtime();
  for (int t = 0; t < niters; ++t) {
    stencil(nx, ny,bufferImg, bufferTempImg, rank);
    stencil(nx, ny, bufferTempImg, bufferImg, rank);
  }
  double toc = wtime();

  MPI_Finalize();
  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc-tic);
  printf("------------------------------------\n");

  output_image(OUTPUT_FILE, nx, ny, image);
  free(image);

}


float *extractElements(float *subArray, float *array, int start, int end)
{

  for (int i = start; i <= end; i++)
  {
    subArray[i - start] = array[i];
  }

  return subArray;
}

void stencil(const int nx, const int ny,  float *restrict image, float *restrict tmp_image, int rank) {

  if(rank==0){
    int start = (ny-1) * nx;
    int end   = (ny-1) * nx + (nx-1);

    float * lastRowSend = (float * ) malloc(nx*sizeof(float));
    float * lastRowRecv = (float * ) malloc(nx*sizeof(float));
    lastRowSend = extractElements(lastRowSend, image, start, end);
    

    MPI_Status *status;
    //MPI_Sendrecv(lastRowSend, nx, MPI_FLOAT, rank +1, 0, lastRowRecv, nx, MPI_FLOAT, rank+1,0, MPI_COMM_WORLD, status);
    MPI_Send(lastRowSend, nx ,MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD );
    MPI_Recv(lastRowRecv, nx ,MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD,status );





  }
  if(rank==1){
    int start = (ny-1) * nx;
    int end   = (ny-1) * nx + (nx-1);

    float * lastRowSend = (float * ) malloc(nx*sizeof(float));
    float * lastRowRecv = (float * ) malloc(nx*sizeof(float));
    lastRowSend = extractElements(lastRowSend, image, start, end);
    

    MPI_Status *status;
    //MPI_Sendrecv(lastRowSend, nx, MPI_FLOAT, rank -1, 0, lastRowRecv, nx, MPI_FLOAT, rank-1,0, MPI_COMM_WORLD, status);
    MPI_Send(lastRowSend, nx ,MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD );
    MPI_Recv(lastRowRecv, nx ,MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD,status );


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
