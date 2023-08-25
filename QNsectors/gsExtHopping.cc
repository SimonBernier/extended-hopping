//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int main(int argc, char* argv[]){

    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_exthubbard",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");

    auto nsweeps = input.getInt("nsweeps");
    auto g = input.getReal("g",0.);
    auto V = input.getReal("V",1.);
    auto Jz = input.getReal("Jz",1.);
    auto h = input.getReal("h",1.);
    auto silent = input.getYesNo("silent",false);

    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
    println(sweeps);

    // write results to file
    char schar[128];
    int n1 = std::sprintf(schar,"N_%d_g_%0.1f_V_%0.1f_Jz_%0.1f_h_%0.1f.dat",N,g,V,Jz,h);
    std::string s1(schar);
    std::ofstream datafile;
    datafile.open(s1); // opens the file
    if( !datafile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    datafile << "Npart" << " " << "Sz" << " " << "energy" << " " << "variation" << std::endl; 

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = tJ(N);

    //
    // Create the Hamiltonian using AutoMPO
    //
    auto ampo = AutoMPO(sites);
    for(int b = 1; b <= N; b++){
            if(b%2==0){
                ampo += -h,"Sz",b;
            }
            else{
                ampo += +h,"Sz",b;
            }
        }
    for(int b = 1; b < N; b++){
        ampo += -1.,"Cdagup",b,"Cup",b+1;
        ampo += -1.,"Cdagup",b+1,"Cup",b;
        ampo += -1.,"Cdagdn",b,"Cdn",b+1;
        ampo += -1.,"Cdagdn",b+1,"Cdn",b;
        ampo += V, "Ntot",b,"Ntot",b+1;
        ampo += Jz, "Sz",b,"Sz",b+1;
    }
    for(int b = 1; b < N-1; b++){
        ampo += -1.,"Cdagup",b,"Cup",b+2;
        ampo += -1.,"Cdagup",b+2,"Cup",b;
        ampo += -1.,"Cdagdn",b,"Cdn",b+2;
        ampo += -1.,"Cdagdn",b+2,"Cdn",b;
        ampo += V, "Ntot",b,"Ntot",b+2;
        ampo += V, "Jz",b,"Jz",b+2;
    }
    for(int b = 1; b < N-1; b++){
        ampo += g,"Cdagup",b,"Ntot",b+1,"Cup",b+2;
        ampo += g,"Cdagup",b+2,"Ntot",b+1,"Cup",b;
        ampo += g,"Cdagdn",b,"Ntot",b+1,"Cdn",b+2;
        ampo += g,"Cdagdn",b+2,"Ntot",b+1,"Cdn",b;
    }
    auto H = toMPO(ampo);
    
    // loop through all states
    for( int j = 0; j<=N; j++){
        printfln("N = %d particles\n",j);

        for( int k = 0; k <= j; k++){

            //
            // Set the initial wavefunction matrix product state
            //
            auto state = InitState(sites);
            int ndn = j-k;
            int p = j;
            for(int i = N; i >= 1; --i) 
                {
                if(p > i)
                    {
                    state.set(i,"UpDn");
                    p -= 2;
                    ndn -= 1;
                    }
                else
                if(p > 0)
                    {
                    state.set(i,(ndn > 0 ? "Dn" : "Up"));
                    p -= 1;
                    ndn -= 1;
                    }
                else
                    {
                    state.set(i,"Emp");
                    }
                }

            auto psi0 = MPS(state);

            Print(totalQN(psi0));

            //
            // Begin the DMRG calculation
            //
            auto [energy,psi] = dmrg(H,psi0,sweeps,{"Silent",silent});
            auto var = sqrt( abs( inner(H,psi,H,psi) - energy*energy ) );

            /*
            //
            // Measure spin densities
            //
            Vector upd(N),dnd(N);
            for(int j = 1; j <= N; ++j)
                {
                psi.position(j);
                upd(j-1) = elt(dag(prime(psi(j),"Site"))*op(sites,"Nup",j)*psi(j));
                dnd(j-1) = elt(dag(prime(psi(j),"Site"))*op(sites,"Ndn",j)*psi(j));
                }

            println("Up Density:");
            for(int j = 0; j < N; ++j)
                printfln("%d %.10f",1+j,upd(j));
            println();

            println("Dn Density:");
            for(int j = 0; j < N; ++j)
                printfln("%d %.10f",1+j,dnd(j));
            println();

            println("Total Density:");
            for(int j = 0; j < N; ++j)
                printfln("%d %.10f",1+j,(upd(j)+dnd(j)));
            println();

            */
            //
            // Print the final energy reported by DMRG
            //
            printfln("Ground State Energy = %.10f\n",energy);

            datafile << j << " " << 2*k-j << " " << energy << " " << var << std::endl;
            
        }
    }

    datafile.close();

    return 0;
    }
