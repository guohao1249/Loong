#include <iostream>
#include <stdio.h>
#include "stdlib.h"
#include "time.h"
#include <math.h>
#include <random>

using namespace std;

// S-box
const int S[16] =
{0Xc, 0xa, 0xd, 0x3, 0xe, 0xb, 0xf, 0x7, 0x9, 0x8, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6};

//MC, MR
int M1_M[16] = {1,4,9,13, 4,1,13,9, 9,13,1,4, 13,9,4,1};
int M2_M[16] = {13,9,4,1, 9,13,1,4, 4,1,13,9, 1,4,9,13};

const int RC9 = 1;
const int RC10 = 2;
const int RC11 = 4;

const int RC12[33] =
{0x0, 0x0, 0x0, 0x1, 0x3, 0x7, 0x7, 0x7, 0x6, 0x5, 0x3,
 0x7, 0x7, 0x6, 0x4, 0x1, 0x3, 0x7, 0x6, 0x5, 0x2, 0x5,
 0x3, 0x6, 0x4, 0x0, 0x0, 0x1, 0x2, 0x5, 0x3, 0x7, 0x6 };

const int RC13[33] =
{0x1, 0x3, 0x7, 0x7, 0x7, 0x6, 0x5, 0x3, 0x7, 0x7, 0x6,
 0x4, 0x1, 0x3, 0x7, 0x6, 0x5, 0x2, 0x5, 0x3, 0x6, 0x4,
 0x0, 0x0, 0x1, 0x2, 0x5, 0x3, 0x7, 0x6, 0x4, 0x0, 0x1 };

const int RC14[33] =
{0x0, 0x0, 0x0, 0x1, 0x3, 0x7, 0x7, 0x7, 0x6, 0x5, 0x3,
 0x7, 0x7, 0x6, 0x4, 0x1, 0x3, 0x7, 0x6, 0x5, 0x2, 0x5,
 0x3, 0x6, 0x4, 0x0, 0x0, 0x1, 0x2, 0x5, 0x3, 0x7, 0x6 };

const int RC15[33] =
{0x1, 0x3, 0x7, 0x7, 0x7, 0x6, 0x5, 0x3, 0x7, 0x7, 0x6,
 0x4, 0x1, 0x3, 0x7, 0x6, 0x5, 0x2, 0x5, 0x3, 0x6, 0x4,
 0x0, 0x0, 0x1, 0x2, 0x5, 0x3, 0x7, 0x6, 0x4, 0x0, 0x1 };

//Multiplication over finite fields with irreducible polynomial x^4 + x + 1
int Gmul_4(int a, int b)
{
    int px = 0x3;
    int p = 0;
    int high;

    for (int i =0; i < 4; i++)
    {
        if ( (b & 1) == 1 )
            p ^= a;

        high = a & 0x8;
        a <<= 1;
        a &= 0xf;
        if (high == 0x8)
            a ^= px;
        b >>= 1;
    }
    return p;
}

//Matrix multiplication over finite fields
int M_mul(int A[16], int B[16], int res[16] )
{
    for (int i = 0; i < 16; i++)
        res[i] = 0;
    for (int i = 0; i < 4; i ++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
                res[4*i + j] = res[4*i + j]  ^ Gmul_4(A[4*i + k],B[4*k + j]);
        }
    }
}

int Encrypt_128(int Stt[16], int K0[16],int K1[16], int r)
{
    int res[16] ;

    //XOR white key and RC
    for (int i = 0; i < 16; i++)
    {
        Stt[i] = Stt[i] ^ K0[i];
    }

    Stt[9] ^= RC9;
    Stt[10] ^= RC10;
    Stt[11] ^= RC11;

    Stt[12] ^= RC12[0];
    Stt[13] ^= RC13[0];
    Stt[14] ^= RC14[0];
    Stt[15] ^= RC15[0];

     //Round function
    for (int round = 0; round < r; round++)
    {
        //S-layer
        for (int i = 0; i < 16; i++)
            Stt[i] = S[Stt[i]];
        //L-layer
        M_mul(Stt, M1_M, res);
        for (int i = 0; i < 16; i++)
            Stt[i] = res[i];
        M_mul(M2_M, Stt, res);
        for (int i = 0; i < 16; i++)
            Stt[i] = res[i];
        //S-layer and XOR key and RC
        for (int i = 0; i < 16; i++)
        {
            Stt[i] = S[Stt[i]];
            if (round%2 == 0)
                Stt[i] ^= K1[i];
            if (round%2 == 1)
                Stt[i] ^= K0[i];
        }

        Stt[9] ^= RC9;
        Stt[10] ^= RC10;
        Stt[11] ^= RC11;

        Stt[12] ^= RC12[round + 1];
        Stt[13] ^= RC13[round + 1];
        Stt[14] ^= RC14[round + 1];
        Stt[15] ^= RC15[round + 1];
    }
}




int Test_Condition(int C[16])
{
     // Loong-80
     if (C[ 0] == 0x2 and C[ 1] == 0x2 and C[ 2] == 0xa and C[ 3] == 0x3 and
         C[ 4] == 0x0 and C[ 5] == 0x0 and C[ 6] == 0x0 and C[ 7] == 0x0 and
         C[ 8] == 0x0 and C[ 9] == 0x0 and C[10] == 0x0 and C[11] == 0x0 and
         C[12] == 0x0 and C[13] == 0x0 and C[14] == 0x0 and C[15] == 0x0)
     {
          return 1;
     }
     else return 0;
}

int main()
{
    int P0[16], P1[16], K0[16], K1[16], Diff_Out[16], T[16];
    int correct_pair;
    double Pr = 0;
    double Pr_sum = 0;

    //Differential for Loong-80
     int Diff[16] = {0x2, 0x2, 0xa, 0x3,
                    0x0, 0x0, 0x0, 0x0,
                    0x0, 0x0, 0x0, 0x0,
                    0x0, 0x0, 0x0, 0x0 };

    int R = 2;  // Number of round
    double d = pow(2, 20);      // Number of pairs
    int m_test = 6; //Number of experiments

    mt19937_64 mt_rand(time(0));

    for (int m; m < m_test; m++)
    {
        //random choose key
        for (int i = 0; i < 16; i++)
            K0[i] = mt_rand()%16;
        for (int i = 0; i < 16; i++)
            K1[i] = mt_rand()%16;

        //weak key for Loong-80
        K0[0] = 0;
        K0[1] = 0;
        K0[2] = 0;
        K0[3] = 0;
        K0[5] = 0;
        K1[1] = 0;
        K0[9] = 6;

        //key schedule of Loong-80
        for (int i = 4; i < 16; i++)
        {
            K1[i] = K0[i - 4];
        }

        cout << "Key_80: " << endl;
        for (int i = 0; i < 16; i++)
            cout << hex << K0[i] ;
        cout << endl;
        for (int i = 0; i < 16; i++)
            cout << hex << K1[i] ;
        cout << endl;

        cout << "P0: " << endl;

        //random choose data
        for (double data = 0; data < d; data++)
        {
            for (int i = 0; i < 16; i++)
                P0[i] = mt_rand()%16;

            for (int i = 0; i < 16; i++)
                P1[i] = P0[i] ^ Diff[i];

            for (int i = 0; i < 16; i++)
                T[i] = P0[i];

            Encrypt_128( P0,  K0, K1, R);
            Encrypt_128( P1,  K0, K1, R);

            for (int i = 0; i < 16; i++)
                Diff_Out[i] = P0[i] ^ P1[i];

            if (Test_Condition(Diff_Out))
            {
                correct_pair = correct_pair + 1;

                for(int i = 0; i < 16; i++)
                    cout << hex << T[i];
                cout << endl;
            }
        }

        Pr = float(correct_pair)/float(d);
        Pr_sum = Pr_sum + Pr;

        cout << "right pairs : " ;
        cout << correct_pair << endl;

        cout << "pro : " ;
        cout << log2 (Pr) << endl;

    }
    cout << "average pro: " << log2(Pr_sum/m_test);

    return 0;
}
