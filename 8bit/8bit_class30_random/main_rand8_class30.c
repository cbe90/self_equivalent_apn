/* Compile with: g++ -O2 main_rand8_class30.c -o rand8_class30
Authors: C. Beierle, G. Leander  -- Sep 2020

This program randomly searches for 8-bit APN permutations F s.t. AF = FB with A, B corresponding to class 30, i.e.,
A = B =
[1 0|0 0|0 0|0 0]
[0 1|0 0|0 0|0 0]
[---+---+---+---]
[0 0|0 1|0 0|0 0]
[0 0|1 0|0 0|0 0]
[---+---+---+---]
[0 0|0 0|0 1|0 0]
[0 0|0 0|1 0|0 0]
[---+---+---+---]
[0 0|0 0|0 0|0 1]
[0 0|0 0|0 0|1 0]
up to affine equivalence. The program is optimized for that particular class.
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <random>

// dimension of the APN permutation to find
#define  N 8
#define RE_SHUFFLE 60	// defines the seconds after a fresh random state is generated

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15, 32, 33, 34, 35, 40, 41, 42, 43, 36, 37, 38, 39, 44, 45, 46, 47, 16, 17, 18, 19, 24, 25, 26, 27, 20, 21, 22, 23, 28, 29, 30, 31, 48, 49, 50, 51, 56, 57, 58, 59, 52, 53, 54, 55, 60, 61, 62, 63, 128, 129, 130, 131, 136, 137, 138, 139, 132, 133, 134, 135, 140, 141, 142, 143, 160, 161, 162, 163, 168, 169, 170, 171, 164, 165, 166, 167, 172, 173, 174, 175, 144, 145, 146, 147, 152, 153, 154, 155, 148, 149, 150, 151, 156, 157, 158, 159, 176, 177, 178, 179, 184, 185, 186, 187, 180, 181, 182, 183, 188, 189, 190, 191, 64, 65, 66, 67, 72, 73, 74, 75, 68, 69, 70, 71, 76, 77, 78, 79, 96, 97, 98, 99, 104, 105, 106, 107, 100, 101, 102, 103, 108, 109, 110, 111, 80, 81, 82, 83, 88, 89, 90, 91, 84, 85, 86, 87, 92, 93, 94, 95, 112, 113, 114, 115, 120, 121, 122, 123, 116, 117, 118, 119, 124, 125, 126, 127, 192, 193, 194, 195, 200, 201, 202, 203, 196, 197, 198, 199, 204, 205, 206, 207, 224, 225, 226, 227, 232, 233, 234, 235, 228, 229, 230, 231, 236, 237, 238, 239, 208, 209, 210, 211, 216, 217, 218, 219, 212, 213, 214, 215, 220, 221, 222, 223, 240, 241, 242, 243, 248, 249, 250, 251, 244, 245, 246, 247, 252, 253, 254, 255};

// Inverse of A
int A_inv[(1<<N)] = {0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7, 12, 13, 14, 15, 32, 33, 34, 35, 40, 41, 42, 43, 36, 37, 38, 39, 44, 45, 46, 47, 16, 17, 18, 19, 24, 25, 26, 27, 20, 21, 22, 23, 28, 29, 30, 31, 48, 49, 50, 51, 56, 57, 58, 59, 52, 53, 54, 55, 60, 61, 62, 63, 128, 129, 130, 131, 136, 137, 138, 139, 132, 133, 134, 135, 140, 141, 142, 143, 160, 161, 162, 163, 168, 169, 170, 171, 164, 165, 166, 167, 172, 173, 174, 175, 144, 145, 146, 147, 152, 153, 154, 155, 148, 149, 150, 151, 156, 157, 158, 159, 176, 177, 178, 179, 184, 185, 186, 187, 180, 181, 182, 183, 188, 189, 190, 191, 64, 65, 66, 67, 72, 73, 74, 75, 68, 69, 70, 71, 76, 77, 78, 79, 96, 97, 98, 99, 104, 105, 106, 107, 100, 101, 102, 103, 108, 109, 110, 111, 80, 81, 82, 83, 88, 89, 90, 91, 84, 85, 86, 87, 92, 93, 94, 95, 112, 113, 114, 115, 120, 121, 122, 123, 116, 117, 118, 119, 124, 125, 126, 127, 192, 193, 194, 195, 200, 201, 202, 203, 196, 197, 198, 199, 204, 205, 206, 207, 224, 225, 226, 227, 232, 233, 234, 235, 228, 229, 230, 231, 236, 237, 238, 239, 208, 209, 210, 211, 216, 217, 218, 219, 212, 213, 214, 215, 220, 221, 222, 223, 240, 241, 242, 243, 248, 249, 250, 251, 244, 245, 246, 247, 252, 253, 254, 255};

int pos[143] = {4,5,6,7,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,52,53,54,55,64,65,66,67,68,69,70,71,72,73,74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89, 90, 91, 92, 93,  94,  95,  96,  97,  98,  99,  100,  101,  102,  103,  104,  105, 106,  107,  108,  109,  110,  111,  112,  113,  114,  115,  116,  117,  118,  119,  120,  121,  122,  123,  124,  125,  126,  127,  196,  197,  198,  199,  208, 209,  210,  211,  212,  213,  214,  215,  216,  217,  218,  219,  220,  221,  222,  223,  244,  245,  246,  247,1,2,3,12,13,14,15,48,49,50,51,60,61,62,63,192,193,194,195,204,205,206,207,240,241,242,243,252,253,254,255}; // one element of each cycle within the cycle structure of A (without 0 and fixpoints last). The next free position x will be selected from this list of positions

int repr[(1<<N)] = {0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 24, 25, 26, 27, 20, 21, 22, 23, 28, 29, 30, 31, 48, 49, 50, 51, 52, 53, 54, 55, 52, 53, 54, 55, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 64, 65, 66, 67, 72, 73, 74, 75, 68, 69, 70, 71, 76, 77, 78, 79, 96, 97, 98, 99, 104, 105, 106, 107, 100, 101, 102, 103, 108, 109, 110, 111, 80, 81, 82, 83, 88, 89, 90, 91, 84, 85, 86, 87, 92, 93, 94, 95, 112, 113, 114, 115, 120, 121, 122, 123, 116, 117, 118, 119, 124, 125, 126, 127, 192, 193, 194, 195, 196, 197, 198, 199, 196, 197, 198, 199, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 208, 209, 210, 211, 216, 217, 218, 219, 212, 213, 214, 215, 220, 221, 222, 223, 240, 241, 242, 243, 244, 245, 246, 247, 244, 245, 246, 247, 252, 253, 254, 255}; // the cycle representatives of each element from 0-255.

// numbers of even Hamming weight up to 255 (needed for more efficient APN check)
int evens[127] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126, 129, 130, 132, 135, 136, 139, 141, 142, 144, 147, 149, 150, 153, 154, 156, 159, 160, 163, 165, 166, 169, 170, 172, 175, 177, 178, 180, 183, 184, 187, 189, 190, 192, 195, 197, 198, 201, 202, 204, 207, 209, 210, 212, 215, 216, 219, 221, 222, 225, 226, 228, 231, 232, 235, 237, 238, 240, 243, 245, 246, 249, 250, 252, 255};

int P[(1<<N)];	// contains the random initial state

int solutions;
long long iterations;
int proceed;

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

// random number generators to generate integers uniformly in [0,(1<<N))
std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_int_distribution<int> dis(0,(1<<N)-1);

// applys the Fisher-Yates shuffle to the array arr to get a uniform permutation
void shuffle_array(int* arr) {
    uint32_t j = 0;
    uint32_t tmp = 0;
    for (int i=(1<<N)-1; i!=0; i--) {
        decltype(dis.param()) nrange (0, i);
        dis.param(nrange);
        j = dis(gen);
        // swap
        tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }
}

void printArray(int A[1<<N]) {
	for (int i=0;i<(1<<N);i++)
		printf("%2x ",i);
	printf("\n");
	for (int i=0;i<(1<<N);i++)
		if (A[i]!=-1) printf("%2x ",A[i]);
			else printf("-- ");
	printf("\n\n");
}

void printSkippedArray(int A[128]) {
	for (int i=0;i<32;i++)
		printf("%d ",A[i]);
	printf("\n\n");
}

void printArray2(int A[1<<N]) {
	printf("{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		printf(",0x%02x",A[i]);
	printf("};\n");
}

void fprintArray2(FILE *fp, int A[1<<N]) {
	fprintf(fp, "{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		fprintf(fp, ",0x%02x",A[i]);
	fprintf(fp, "};\n");
}

void printDDT(int A[1<<N][1<<N]) {
    for (int i=0; i<(1<<N); i++) {
        for (int j=0; j<(1<<N); j++) {
            printf("%d ",A[i][j]);
        }
        printf("\n");
    }
}

// checks if A does not contain -1
int isComplete(int A[1<<N]) {
    int val = A[(1<<N)-1];
    int i;
    if (val==-1) return 0;
    A[(1<<N)-1] = -1;
    for (i=0; A[i] != -1; i++);
    A[(1<<N)-1] = val;
	return i==((1<<N)-1);

}

// multiplicative order of x for matrix
int order_mat(int A[1<<N], int x) {
	int y=A[x];
	int i=1;
	while(y!=x) {
		y=A[y];
		i++;
	}
	return i;
}

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for y=0, since it is not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<143; i++) {
		if (repr[sbox[pos[i]]] == repr[y])
			return 0;
	}
	return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<127;i++) {
		a = evens[i];	// we only need to check those input differences corresponding to vectors with even Hamming Weight
                // if sbox[c + a] is set, increase DDT entry by 2
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]+=2;
                    // if set above 2, return that it cannot be APN anymore
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]>2) return 0;
                }

	}
	return 1;
}

// reducing DDT entries for new removed point sbox[c]. It reduces all entries that were increased during the APN check of addDDTInformation(c).
void removeDDTInformation(int c) {
	int a;
	for (int i=0;i<127;i++) {
		a = evens[i];
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]-=2;
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]==2) break; // important because addDDTInformation(c) did not increment any values above that point.
                }

	}
}

int maxDepth=0;

// recursive construction of the sbox.
void nextValue(int depth, char* filename) {
	int xS,yS,j,y;
	iterations++;

	if (depth>maxDepth) maxDepth=depth; // to show maxDepth
	if ((time(NULL)-start)>=RE_SHUFFLE) { // display the status every RE_SHUFFLE seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli \n",iterations);
		printf("running time: %li sec\n", time(NULL)-runtime);
		fflush(stdout);
		proceed = 0;
	}

	//complete
	if (isComplete(sbox)) {
			solutions++;
			printf("found a new apn permutation: #%d\n",solutions);
			printArray2(sbox);
			fflush(stdout);
			FILE *fp = fopen(filename, "a");
			if (fp == NULL)
            {
                printf("Error opening file!\n");
		fflush(stdout);
            }
            else
            {
                fprintArray2(fp, sbox);
            }
            fclose(fp);
		return;
	}

	//not complete

	int x=pos[depth];	// get next free position in the look-up table
	for (int z=0;z<(1<<N);z++) {
		y = P[z];
		if (y!=0) { // do not include y=0 since this value is already set in the beginning. The function isNotTaken(y) would not work for 0.
		if (isNotTaken(y) && order_mat(A,x)==order_mat(A,y)) { // if the orders are different, it cannot be of the desired form
			xS = x;
			yS = y;
			for (int i=0; i<order_mat(A,x)-1; i++) {
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO; // undo if it can't be APN anymore
				xS=A[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO; // undo if it can't be APN anymore

                    	if(proceed)	nextValue(depth+1, filename);

			//undo the changes (i.e., delete the set points of sbox from this step and undo the changes in the DDT)
			UNDO:
			// we remove the points in reverse order to undo the DDT changes correctly
			do {
				removeDDTInformation(xS);
				sbox[xS] = -1;
				xS = A_inv[xS];
			} while (xS != A_inv[x]);
		}
		}
	}
}

void test_for_matrix(char* filename, int max_runtime) {
	if (N!=8) {
		printf("invalid N\n");
		exit(-1);
	}
    int xS, yS;
    solutions=0;
    iterations=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);
    printf("Search S-boxes invariant for matrix: \n");
    printArray2(A);
    printf("\n\n\n");
    fflush(stdout);
    for (int x=0;x<(1<<N);x++) {
	P[x] = x;
    }

    while(time(NULL)-runtime < max_runtime) {
	memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
	for (int x=0;x<(1<<N);x++) {
        	sbox[x]=-1;
    	}
	// w.l.o.g., we can set sbox[0]=0
    	sbox[0]=0;
    	if (!addDDTInformation(0)) exit(0);

    	// shuffle the state
    	shuffle_array(P);
    	printf("Re-shuffled intitial state ..\n");
    	fflush(stdout);

    	// start the recursive search
    	proceed = 1;
    	nextValue(0,filename);
    }
}


int main(int argc, char* argv[])
{
	if (argc != 3) {
        	printf("Usage: %s <outfile> <max_runtime> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);

	test_for_matrix(argv[1],atoi(argv[2]));
	printf("\nFinished. Total running time: %li sec\n", time(NULL)-runtime);
	exit(0);
}
