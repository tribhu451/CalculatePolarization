#include<iostream>
#include<cmath>

using namespace std;

int Sgn(double n){ return (n > 0.) - (n < 0.); }


int main(int argc, char *argv[]){

    int a[4];
    
    for (int ia1 = 0; ia1 < 4; ia1++) 
        for (int ia2 = 0; ia2 < 4; ia2++) 
            for (int ia3 = 0; ia3 < 4; ia3++) 
                for (int ia4 = 0; ia4 < 4; ia4++)
                {
                   a[0] = ia1;
                   a[1] = ia2;
                   a[2] = ia3;
                   a[3] = ia4;
                   int val = 1;
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 3; j > i; j--)
                            {
                                val *= Sgn(a[j] - a[i]);
                            }
        
                        }
                if(val != 0)std :: cout  <<"eps_{"<< ia1 << ia2 << ia3 << ia4 << "} = " << val << std ::endl;
                }
    return 0;

}