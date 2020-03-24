/* Compile with: g++ -O2 main_det7_class17.c -o det7_class17
Authors: C. Beierle, G. Leander  -- Feb 2020

This program constructs and returns all 7-bit APN permutations F s.t. AF = FB with A, B corresponding to class 17, i.e.,
A = 
0 0 1 0 0 0 0
1 0 1 0 0 0 0
0 1 0 0 0 0 0
0 0 0 0 0 0 1
0 0 0 1 0 0 0
0 0 0 0 1 0 1
0 0 0 0 0 1 1,
B =
0 0 0 0 0 0 1
1 0 0 0 0 0 0
0 1 0 0 0 0 0
0 0 1 0 0 0 0
0 0 0 1 0 0 0
0 0 0 0 1 0 0
0 0 0 0 0 1 0
up to affine equivalence. It may return multiple representatives in one affine-equivalence class. 
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include "commuting_matrices_A.h"
#include "commuting_matrices_B.h"

// dimension of the APN permutation to find
#define  N 7

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 2, 4, 6, 3, 1, 7, 5, 16, 18, 20, 22, 19, 17, 23, 21, 32, 34, 36, 38, 35, 33, 39, 37, 48, 50, 52, 54, 51, 49, 55, 53, 64, 66, 68, 70, 67, 65, 71, 69, 80, 82, 84, 86, 83, 81, 87, 85, 96, 98, 100, 102, 99, 97, 103, 101, 112, 114, 116, 118, 115, 113, 119, 117, 104, 106, 108, 110, 107, 105, 111, 109, 120, 122, 124, 126, 123, 121, 127, 125, 72, 74, 76, 78, 75, 73, 79, 77, 88, 90, 92, 94, 91, 89, 95, 93, 40, 42, 44, 46, 43, 41, 47, 45, 56, 58, 60, 62, 59, 57, 63, 61, 8, 10, 12, 14, 11, 9, 15, 13, 24, 26, 28, 30, 27, 25, 31, 29};

// the matrix B
int B[(1<<N)] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127};

// inverse of B
int B_inv[(1<<N)] = {0, 64, 1, 65, 2, 66, 3, 67, 4, 68, 5, 69, 6, 70, 7, 71, 8, 72, 9, 73, 10, 74, 11, 75, 12, 76, 13, 77, 14, 78, 15, 79, 16, 80, 17, 81, 18, 82, 19, 83, 20, 84, 21, 85, 22, 86, 23, 87, 24, 88, 25, 89, 26, 90, 27, 91, 28, 92, 29, 93, 30, 94, 31, 95, 32, 96, 33, 97, 34, 98, 35, 99, 36, 100, 37, 101, 38, 102, 39, 103, 40, 104, 41, 105, 42, 106, 43, 107, 44, 108, 45, 109, 46, 110, 47, 111, 48, 112, 49, 113, 50, 114, 51, 115, 52, 116, 53, 117, 54, 118, 55, 119, 56, 120, 57, 121, 58, 122, 59, 123, 60, 124, 61, 125, 62, 126, 63, 127};

// an ordering of the indices mod 2^n. It orders cycle-wise with respect to the cycles of B
int lex[(1<<N)] = {0,1, 2, 4, 8, 16, 32, 64,3, 6, 12, 24, 48, 96, 65,5, 10, 20, 40, 80, 33, 66,7, 14, 28, 56, 112, 97, 67,9, 18, 36, 72, 17, 34, 68,11, 22, 44, 88, 49, 98, 69,13, 26, 52, 104, 81, 35, 70,15, 30, 60, 120, 113, 99, 71, 19, 38, 76, 25, 50, 100, 73, 21, 42, 84, 41, 82, 37, 74, 23, 46, 92, 57, 114, 101, 75, 27, 54, 108, 89, 51, 102, 77, 29, 58, 116, 105, 83, 39, 78, 31, 62, 124, 121, 115, 103, 79, 43, 86, 45, 90, 53, 106, 85, 47, 94, 61, 122, 117, 107, 87, 55, 110, 93, 59, 118, 109, 91, 63, 126, 125, 123, 119, 111, 95,127};

int pos[18] = {1,3,5,7,9,11,13,15,19,21,23,27,29,31,43,47,55,63}; // one element of each cycle within the cycle structure of B (without 0 and 127). The next free position x will be selected from this list of positions

int repr[(1<<N)] = {0, 1, 1, 1, 1, 1, 1, 1, 8, 9, 10, 11, 12, 13, 14, 15, 8, 13, 9, 12, 10, 15, 11, 14, 24, 25, 26, 27, 28, 29, 30, 31, 8, 15, 13, 10, 9, 14, 12, 11, 24, 30, 31, 25, 29, 27, 26, 28, 24, 29, 25, 28, 26, 31, 27, 30, 8, 12, 11, 15, 14, 10, 13, 9, 8, 14, 15, 9, 13, 11, 10, 12, 24, 28, 27, 31, 30, 26, 29, 25, 24, 27, 30, 29, 31, 28, 25, 26, 88, 89, 89, 89, 89, 89, 89, 89, 24, 31, 29, 26, 25, 30, 28, 27, 8, 11, 14, 13, 15, 12, 9, 10, 8, 10, 12, 14, 11, 9, 15, 13, 24, 26, 28, 30, 27, 25, 31, 29}; // the cycle representatives of cycles in A of each element from 0-127.

// numbers of even Hamming weight up to 127 (needed for more efficient APN check)
int evens[63] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126};

int solutions;
long long iterations;
int n_skipped;
int skipped[128];

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

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


// returns 1 if A is the smallest permutation up to conjugation with matrices in commuting_matrices. We use the lexicographic order of the lookup tables.
int is_smallest_in_class(int A[1<<N]){
    int E;
    for (int i=0;i<N_COMMUTING_MATRICES_B;i++) {
        for (int j=0; j<N_COMMUTING_MATRICES_A;j++) {
            for (int x=0;x<(1<<N);x++){
                if (A[lex[x]]==-1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                if (A[commuting_matrices_B[i][lex[x]]] == -1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                else {
                    E = commuting_matrices_A[j][A[commuting_matrices_B[i][lex[x]]]];
                    if (E < A[lex[x]]) return 0; // E is smaller, so A cannot be the smallest in class
                    if (E > A[lex[x]]) goto NEXT_MAT; // A is smaller, so consider next matrix
                } // consider next position
            }
            NEXT_MAT:;
        }
    }
    // A is the smallest of all, so return 1
    return 1;
}

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for y=0 (and y=88), since 0 and 127 are not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<18; i++) {
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
void nextValue_with_filter(int depth, char* filename) {
	int xS,yS,j;
	iterations++;

	if (depth>maxDepth) maxDepth=depth; // to show maxDepth
	if ((time(NULL)-start)>=300) { // display the status every 300 seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli  skipped: %d \nDepth skipped: \n",iterations,n_skipped);
		printSkippedArray(skipped);
		printf("running time: %li sec\n", time(NULL)-runtime);
		fflush(stdout);
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
	for (int y=(1<<N)-1;y>=1;y--) {		// do not include y=0 and y=88 since those values are already set in the beginning. The function isNotTaken(y) would not work for 0.
	if (y != 88) {
		if (isNotTaken(y)) { // if the orders are different, it cannot be of the desired form (we know the order of each element is 7)
			xS = x;
			yS = y;
			for (int i=0; i<6; i++) {	// we know the order of each element is 7
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
				xS=B[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
			
                	if (depth <= 4) { // if it is too deep, it is faster to omit the check for smallest representative.
                        	if (is_smallest_in_class(sbox)) nextValue_with_filter(depth+1, filename);
                        	else {n_skipped++; skipped[depth]++;}
                	}
                	else { // omit check for smallest representative
                    		nextValue_with_filter(depth+1, filename);
                	}
            
			//undo the changes (i.e., delete the set points of sbox from this step and undo the changes in the DDT)
			UNDO_WITH_FILTER:
			// we remove the points in reverse order to undo the DDT changes correctly
			do { 
				removeDDTInformation(xS);
				sbox[xS] = -1;
				xS = B_inv[xS];
			} while (xS != B_inv[x]);
		}
	}
	}
}

void test_for_matrix(char* filename, int first) {
	if (N!=7) {
		printf("invalid N\n");
		exit(-1);
	}
    int xS, yS;
    solutions=0;
    iterations=0;
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int i=0; i<128; i++) skipped[i]=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);
    printf("Search S-boxes invariant for matrices: \n");
    printArray2(B);
    printArray2(A);
    printf("\n\n\n");
    fflush(stdout);
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox[0]=0
    sbox[0]=0;
    if (!addDDTInformation(0)) exit(0);
    // w.l.o.g., we can set sbox[127]=88
    sbox[127]=88;
    if (!addDDTInformation(127)) exit(0);

    sbox[1]=first;
    if (!addDDTInformation(1)) exit(0);
    xS=B[1];
    yS=A[first];
    // set the orbits accordingly
    for (int i=0;i<6;i++) {
                    sbox[xS]=yS;
                    if (!addDDTInformation(xS)) exit(0);
                    xS=B[xS];
                    yS=A[yS];
    } 
    if (!is_smallest_in_class(sbox)) {
	exit(0);
    } 

    // start the recursive search
    nextValue_with_filter(1, filename);
}

int main(int argc, char* argv[])
{
	if (argc != 3) {
        	printf("Usage: %s <outfile> <SBox[1]> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);
	
	
	if ((atoi(argv[2]) == 0) or (atoi(argv[2]) == 88)) {
		printf("No permutation ..\n");
		return -1;
	}
	test_for_matrix(argv[1],atoi(argv[2]));
	printf("\nFinished. Total running time: %li sec\n", time(NULL)-runtime);
	exit(0);
}
