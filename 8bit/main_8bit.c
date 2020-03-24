/* Compile with: g++ -O2 main_8bit.c -o 8bit
Authors: C. Beierle, G. Leander  -- Feb 2020

This program tries all possible tuples for (A,B) and searches for APN permutations F with AF = FB.
It outputs the tuples for which there can be solutions by checking each of them for TIMEOUT seconds.
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include "AllTuples.h"
#include "AllTuples_inv.h"

// dimension of the APN permutation to find
#define  N 8

#define TIMEOUT 900

time_t start;

int C[(1<<N)];	// will store matrix A
int D[(1<<N)];	// will store matrix B
int D_inv[(1<<N)];	// will store inverse of matrix B

int arr[(1<<N)]; // stores the order of the initial state which determines which positions in the look-up table are set first 

// numbers of even Hamming weight up to 255 (needed for more efficient APN check)
int evens[127] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126, 129, 130, 132, 135, 136, 139, 141, 142, 144, 147, 149, 150, 153, 154, 156, 159, 160, 163, 165, 166, 169, 170, 172, 175, 177, 178, 180, 183, 184, 187, 189, 190, 192, 195, 197, 198, 201, 202, 204, 207, 209, 210, 212, 215, 216, 219, 221, 222, 225, 226, 228, 231, 232, 235, 237, 238, 240, 243, 245, 246, 249, 250, 252, 255};

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

int is_possible = 0; // stores whether the current tuple is possible (=1) or impossible (=0)

// sets the state to [0,1,...,2^N-1]
void set_standard_order() {
	for (int i=0; i<(1<<N); i++) {
		arr[i] = i;
	}
}

// sets the state to [s,s+1,...,2^N-1,0,1,...,s-1]
void set_offset(int s) {
	for (int i=0; i<(1<<N); i++) {
		arr[i] = (i+s) % (1<<N);
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

void printArray2(int A[1<<N]) {
	printf("{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		printf(",0x%02x",A[i]);
	printf("};\n");
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

// should only by run if sbox is not complete (i.e., contains -1)
int getNextFreePosition(){
	int i;
	for (i=0; sbox[arr[i]]!=-1;i++);
	return arr[i];
}

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not)
int isNotTaken(int y) {
	int val = sbox[(1<<N)-1];
    int i;
    if (val==y) return 0;
    sbox[(1<<N)-1] = y;
    for (i=0; sbox[i] != y; i++);
    sbox[(1<<N)-1] = val;
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

// recursive construction of the sbox.
void nextValue(int depth) {
	int xS,yS;

	if ((time(NULL)-start)>TIMEOUT) { // skip this tuple after the timeout is reached
		printf("TIMEOUT! Go to next tuple\n");
		is_possible = 1;

		return;
	}

	//complete
	if (isComplete(sbox)) {
			printf("found a new apn permutation:\n");
			printArray2(sbox);
			is_possible = 1;

		return;
	}

	//not complete

	int x=getNextFreePosition();	// get next free position in the look-up table
	for (int y=0;y<(1<<N);y++) {		
		if (isNotTaken(y) && order_mat(D,x)==order_mat(C,y)) { 
			xS = x;
			yS = y;
			for (int i=0; i<order_mat(D,x)-1; i++) {	
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO; // undo if it can't be APN anymore
				xS=D[xS];
                    		yS=C[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO; // undo if it can't be APN anymore
			
                    	if (!is_possible) nextValue(depth+1);
            
			//undo the changes (i.e., delete the set points of sbox from this step and undo the changes in the DDT)
			UNDO:
			// we remove the points in reverse order to undo the DDT changes correctly
			do { 
				removeDDTInformation(xS);
				sbox[xS] = -1;
				xS = D_inv[xS];
			} while (xS != D_inv[x]);
		}
	}
}


void all8() {
	if (N!=8) {
		printf("invalid N\n");
		exit(-1);
	}
	for (int t=0;t<N_TUPLES;t++) {
		set_standard_order();
		if (t==30) {
			set_offset(16);
		}
		if (t==26) {
			set_offset(32);
		}
		if (t==23) {
			set_offset(16);
		}
		if (t==22) {
			set_offset(16);
		}
		if (t==12) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
		if (t==17) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
		if (t==18) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
		if (t==19) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
		if (t==20) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
		if (t==24) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
		if (t==25) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
		if (t==27) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
		if (t==31) {
			// tuple can never yield an APN permutation because there is an i s.t. dim Ord(A,i) in {2,4,N-1}
			printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
			goto NEXT_TUPLE;
		}
        	is_possible = 0;
		memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
		for (int x=0;x<(1<<N);x++) {
			sbox[x]=-1;
			C[x]=AllTuples[t][x];
			D[x]=AllTuples[t][x+(1<<N)];
			D_inv[x]=AllTuples_inv[t][x+(1<<N)];
		}
		sbox[0]=0;
		if (!addDDTInformation(0)) goto NEXT_TUPLE;

        	start=time(NULL);
        	nextValue(0);
        	if (is_possible) printf("t = %d might be a possible tuple and needs further consideration\n", t+1); 
		else printf("t = %d is impossible\n", t+1); // in the paper, class index starts from 1 (therefore t+1)
        	NEXT_TUPLE:;
	}
}


int main(int argc, char* argv[])
{
	srand (time(NULL));
	start= time(NULL);
	all8();
	exit(0);
}
