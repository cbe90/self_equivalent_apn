/* Compile with: g++ -O2 main_rand7_class22.c -o rand7_class22
Authors: C. Beierle, G. Leander  -- Feb 2020

This program randomly searches for 7-bit APN permutations F s.t. AF = FB with A, B corresponding to class 22, i.e.,
A = B = I_1 \oplus Comp(X^3+1) \oplus Comp(X^3+1)
up to affine equivalence. It may return multiple representatives in one affine-equivalence class. 
*/


#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <random>

// dimension of the APN permutation to find
#define  N 7
#define RE_SHUFFLE 60	// defines the seconds after a fresh random state is generated

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15, 32, 33, 36, 37, 40, 41, 44, 45, 34, 35, 38, 39, 42, 43, 46, 47, 64, 65, 68, 69, 72, 73, 76, 77, 66, 67, 70, 71, 74, 75, 78, 79, 96, 97, 100, 101, 104, 105, 108, 109, 98, 99, 102, 103, 106, 107, 110, 111, 16, 17, 20, 21, 24, 25, 28, 29, 18, 19, 22, 23, 26, 27, 30, 31, 48, 49, 52, 53, 56, 57, 60, 61, 50, 51, 54, 55, 58, 59, 62, 63, 80, 81, 84, 85, 88, 89, 92, 93, 82, 83, 86, 87, 90, 91, 94, 95, 112, 113, 116, 117, 120, 121, 124, 125, 114, 115, 118, 119, 122, 123, 126, 127};

// Inverse of A
int A_inv[(1<<N)] = {0, 1, 8, 9, 2, 3, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15, 64, 65, 72, 73, 66, 67, 74, 75, 68, 69, 76, 77, 70, 71, 78, 79, 16, 17, 24, 25, 18, 19, 26, 27, 20, 21, 28, 29, 22, 23, 30, 31, 80, 81, 88, 89, 82, 83, 90, 91, 84, 85, 92, 93, 86, 87, 94, 95, 32, 33, 40, 41, 34, 35, 42, 43, 36, 37, 44, 45, 38, 39, 46, 47, 96, 97, 104, 105, 98, 99, 106, 107, 100, 101, 108, 109, 102, 103, 110, 111, 48, 49, 56, 57, 50, 51, 58, 59, 52, 53, 60, 61, 54, 55, 62, 63, 112, 113, 120, 121, 114, 115, 122, 123, 116, 117, 124, 125, 118, 119, 126, 127};

int pos[40] = {2,3,6,7,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,114,115,118,119}; // one element of each cycle within the cycle structure of A (without the fixed points). The next free position x will be selected from this list of positions (long cycles indexed first)

int repr[(1<<N)] = {0, 1, 2, 3, 2, 3, 6, 7, 2, 3, 6, 7, 6, 7, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 24, 25, 18, 19, 26, 27, 20, 21, 28, 29, 22, 23, 30, 31, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 16, 17, 20, 21, 24, 25, 28, 29, 18, 19, 22, 23, 26, 27, 30, 31, 48, 49, 52, 53, 56, 57, 60, 61, 50, 51, 54, 55, 58, 59, 62, 63, 48, 49, 56, 57, 50, 51, 58, 59, 52, 53, 60, 61, 54, 55, 62, 63, 112, 113, 114, 115, 114, 115, 118, 119, 114, 115, 118, 119, 118, 119, 126, 127}; // the cycle representatives of each element from 0-127.

// numbers of even Hamming weight up to 127 (needed for more efficient APN check)
int evens[63] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126};

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

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for the fixed points, since they are not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<40; i++) {
		if (repr[sbox[pos[i]]] == repr[y])
			return 0;
	}
	return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<63;i++) {
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
	for (int i=0;i<63;i++) {
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
		if ((y!=0) && (y!=1) && (y!=14) && (y!=15) && (y!=112) && (y!=113) && (y!=126) && (y!=127)) { // fixed points are already set in the beginning
		if (isNotTaken(y)) { // we know the order of each element is 3
			xS = x;
			yS = y;
			for (int i=0; i<2; i++) {	// we know that the order of x is 3
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
	if (N!=7) {
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
	// w.l.o.g., we can set the sbox on the fixed points 0,1,14,15,112,113,126,127 to the only known APN permutation in 3 bit
    	sbox[0]=0;
    	if (!addDDTInformation(0)) exit(0);
    	sbox[1]=1;
    	if (!addDDTInformation(1)) exit(0);
    	sbox[14]=14;
    	if (!addDDTInformation(14)) exit(0);
    	sbox[15]=112;
    	if (!addDDTInformation(15)) exit(0);
    	sbox[112]=15;
    	if (!addDDTInformation(112)) exit(0);
    	sbox[113]=126;
    	if (!addDDTInformation(113)) exit(0);
    	sbox[126]=127;
    	if (!addDDTInformation(126)) exit(0);
    	sbox[127]=113;
    	if (!addDDTInformation(127)) exit(0);
    
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
