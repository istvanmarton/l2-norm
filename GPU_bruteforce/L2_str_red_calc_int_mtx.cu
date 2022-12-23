#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define length 4096
#define NUM_OF_THREADS 16384
#define NUM_OF_BITS 8 * sizeof(unsigned long long int)

__global__ void func(int* d_mtx_to_vec, unsigned long long int steps, unsigned long long int steps_remainder, int *d_L2_vector, int *d_L2_strategy, int iRows, int iCols){
	int i, l;
	int temp_negative[length], temp_positive[length], vect[NUM_OF_BITS], product, L2;
	unsigned long long int number, index, iMax, iMin, iNumofZeros, iNum_temp;

	int logical;
	index = blockIdx.x * blockDim.x + threadIdx.x;

	iMax = (index + 1) *(steps) - 1;
	iMin = index * (steps);
	if(index < (steps_remainder) ) iMax += index + 1;
	else iMax += (steps_remainder);
	if(index <= (steps_remainder)) iMin += index;
	else iMin += (steps_remainder);
//	if(index == 0) printf("index: %llu, iMin: %llu, iMax: %llu\n", index, iMin, iMax);
		 number = iMin;
		 for(l=0; l < iCols; l++) {temp_negative[l] = 0; temp_positive[l] = d_mtx_to_vec[(iRows-1) * iCols + l];}
		 product = 0;
			for(i = 0 ; (iRows - 1) > i; i++){
				iNumofZeros=(unsigned long long int) 1 << i;
				iNum_temp = (unsigned long long int) iNumofZeros << 1;
				logical = ((number+ iNumofZeros)/iNum_temp) % 2;
				vect[i] = (int) 2 * logical - 1;
//if(index == 2) printf("%d, ",vect[i]);
					if(vect[i] > 0){for(l=0; l < iCols; l++){temp_positive[l] += d_mtx_to_vec[i * iCols + l]; }}
					else {for(l=0; l < iCols; l++){temp_negative[l] += d_mtx_to_vec[i * iCols + l]; }}				
			}
			for(l= 0; l < iCols; l++) {product += abs(temp_negative[l]) + abs(temp_positive[l]);}
			L2 = product;
			for(l=0; l<(iRows - 1); l++){d_L2_strategy[index * (iRows - 1) + l] = vect[l];}

     for(number=iMin + 1; number <= iMax; number++){
//if(index == 2) printf("\n");
		 product = 0;
			for(i = 0 ; (iRows - 1) > i; i++){
				iNumofZeros=(unsigned long long int) 1 << i;
				iNum_temp = (unsigned long long int) iNumofZeros << 1;
				if( ((number+ iNumofZeros) % iNum_temp) ==0 ) {vect[i]=-vect[i] ;					
					if(vect[i] > 0){for(l=0; l < iCols; l++){temp_positive[l] += d_mtx_to_vec[i * iCols + l]; temp_negative[l] -= d_mtx_to_vec[i * iCols + l]; }}
					else {for(l=0; l < iCols; l++){temp_positive[l] -= d_mtx_to_vec[i * iCols + l]; temp_negative[l] += d_mtx_to_vec[i * iCols + l]; }}
				break;
				}
//if(index == 2) printf("%d, ",vect[i]);
            		}
	     for(l = 0; l < (iCols ); l++) {product += abs(temp_negative[l]) + abs(temp_positive[l]);}
	     if(product > L2) {L2 = product;
		for(l=0; l<(iRows - 1); l++){d_L2_strategy[index * (iRows - 1) + l] = vect[l];}
		}
     }
d_L2_vector[index] = L2;
}

int** mtx_read(int *iRows, int *iCols, char* fileName){
	int i = 0,j = 0, k = 0;
	int *row, **mtx, value;
	
	mtx = NULL;
	row = NULL;
	
	char g, cNum[256];
	
	FILE *fp;
	fp = fopen(fileName,"r");
	
	do{
		g = fgetc(fp);	
		if((((g - '0') < 10) && ((g - '0') >= 0)) || (g == 'e') || ( g == 'E') || (g == '.') || (g == '+') || (g == '-')) {cNum[i] = g; i++;}
		else {
			cNum[i] = '\0'; 
			if(cNum[0] != '\0') {sscanf(cNum, "%d", &value); j++; i = 0;  row = (int*) realloc(row, j * sizeof(int)); row[j-1] = value;}
			if( ((g == '\n') || (g == EOF)) && (j > 0)){*iCols = j; j = 0; k++; mtx = (int**) realloc(mtx, k * sizeof(int*)); mtx[k-1] = row; row = NULL;}
		}
		
	}while(!feof(fp));
	*iRows = k;
printf("rows: %d, cols: %d\n",*iRows, *iCols); 
	fclose(fp);
return mtx;
}

void fileN(char *fileName, char** argv, int *argc){
	int r;
	FILE *fp;
	if((*argc) < 2) {
		do{
			printf("Please give me a filename: "); 
			r = scanf("%s",fileName);
		}while(r != 1);
	}
	else sprintf(fileName,"%s", argv[1]);

	fp = fopen(fileName, "r");
	if(fp == NULL) {
		do{
			printf("Please give me a filename that exist within this directory: ");
			r = scanf("%s",fileName);
			fp = fopen(fileName, "r");
		}while(fp == NULL);
	}
	fclose(fp);
}

int main(int argc, char *argv[]){
cudaDeviceProp devProp;
cudaGetDeviceProperties(&devProp, 0);
     char fileName[1024];
     fileN(fileName, argv, &argc);     
     int i, j, iMax, iRows, iCols, **mtx, *mtx_to_vec, *d_mtx_to_vec;

     mtx = mtx_read(&iRows, &iCols, fileName);
     mtx_to_vec = (int*)calloc(iRows * iCols, sizeof(int));

		for(i = 0; i < iRows; i++){
			for(j = 0; j < iCols; j++){
				mtx_to_vec[i * iCols + j] = mtx[i][j];
			}
		}

	if(iRows > (NUM_OF_BITS)) {printf("Matrix is too big. The number of rows or columns can not be more than %lu.\n", NUM_OF_BITS); return 0;}
	if(iCols > length) {printf("Matrix is too big. The length variable %d should be bigger or equal than %d.\n", length, iCols); return 0;}
	cudaMalloc((void**)&d_mtx_to_vec, iRows * iCols * sizeof(int));
     	cudaMalloc((void**)&d_mtx_to_vec, iRows * iCols * sizeof(int));
     	unsigned long long int steps, steps_remainder, Inner_num = (unsigned long long int) 1 << (iRows - 1), copyNum;
	copyNum = NUM_OF_THREADS > Inner_num ? Inner_num : NUM_OF_THREADS;
	int *L2_vector, *d_L2_vector, L2_max = 0, *L2_strategy, *d_L2_strategy, num_ofBlock = (int) ceil((float)copyNum/(float)devProp.warpSize), num_ofThread = copyNum < devProp.warpSize ? copyNum : devProp.warpSize;

	steps=Inner_num/copyNum; steps_remainder = Inner_num % copyNum;

	L2_vector = (int*) malloc(copyNum * sizeof(int));
	L2_strategy = (int*) malloc(copyNum * (iRows - 1) * sizeof(int));

	cudaMalloc((void**)&d_L2_vector, copyNum * sizeof(int));
	cudaMalloc((void**)&d_L2_strategy, copyNum * (iRows - 1) * sizeof(int));
	
	cudaMemcpy(d_mtx_to_vec, mtx_to_vec, iRows * iCols * sizeof(int), cudaMemcpyHostToDevice);
printf("num_ofBlock: %d, num_ofThread: %d\n",num_ofBlock,num_ofThread);
	 func<<<num_ofBlock,num_ofThread>>>(d_mtx_to_vec, steps, steps_remainder, d_L2_vector, d_L2_strategy, iRows, iCols);
	cudaMemcpy(L2_vector, d_L2_vector, copyNum * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(L2_strategy, d_L2_strategy, copyNum * (iRows - 1) * sizeof(int), cudaMemcpyDeviceToHost);
	 
	for(i = 0; i < copyNum; i++){
		if(L2_max < L2_vector[i]) {L2_max = L2_vector[i]; iMax = i;}
	}
FILE *fp;
fp = fopen("strategy_L2.txt", "w");
	for(i=0; i<(iRows - 1); i++) {fprintf(fp, "%d\n", L2_strategy[iMax * (iRows - 1) + i]);}
	fprintf(fp,"1\n");
fclose(fp);

	printf("L2 is: %d\n",L2_max);

	free(L2_vector);
	free(L2_strategy);
	free(mtx_to_vec);

	cudaFree(d_L2_vector);
	cudaFree(d_L2_strategy);
	cudaFree(d_mtx_to_vec);

     return 0;     
}
