//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

int main(int argc, char* argv[]){
    std::clock_t tStart = std::clock();

    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_extGap",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");
    auto V = input.getReal("V",0.1);
    auto Jz = input.getReal("Jz",0.);
    auto h = input.getReal("h",1.);

    auto nsweeps = input.getInt("nsweeps");
    auto silent = input.getYesNo("silent",false);

    auto table1 = InputGroup(input,"sweeps1");
    auto sweeps1 = Sweeps(nsweeps,table1);
    auto table2 = InputGroup(input,"sweeps2");
    auto sweeps2 = Sweeps(nsweeps,table2);
    //println(sweeps);

    // write results to file
    char schar[128];
    int n1 = std::sprintf(schar,"N_%d_V_%0.1f_Jz_%0.1f_h_%0.1f_gap.dat",N,V,Jz,h);
    std::string s1(schar);
    std::ofstream datafile;
    datafile.open(s1); // opens the file
    if( !datafile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    datafile << "g" << " " << "en0" << " " << "var0" << " " << "maxDim0" << " " 
                << "gap01" << " " << "var1" << " " << "maxDim1" << " " << "overlap01" << " "
                << "gap12" << " " << "var2" << " " << "maxDim2" << " " << "overlap12" << " " << std::endl; 

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = tJ(N);
    
    //
    // Set the initial wavefunction matrix product state
    //
    auto state = InitState(sites);
    int p = N/2;
    for(int i = N; i >= 1; --i){
        if(p > 0){
            state.set(i,(i%2 == 0 ? "Dn" : "Up"));
            p -= 1;
        }
        else{
            state.set(i,"Emp");
        }
    }

    std::vector<double> g = {0., 0.5, 0.75, 0.9, 1., 1.1, 1.25, 1.5, 2.0};

    printfln("\ng = 0");
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
        ampo += Jz, "Sz",b,"Sz",b+2;
    }
    auto H = toMPO(ampo);

    //
    // calculate ground state
    //
    auto [en0,psi0] = dmrg(H,MPS(state),sweeps1,{"Silent",silent});
    auto var = sqrt( abs( inner(H,psi0,H,psi0) - en0*en0 ) );
    auto maxDim = maxLinkDim(psi0);
    datafile << g[0] << " " << en0 << " " << var << " " << maxDim << " " ;

    //
    // calculate first excited state
    //
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;
    auto [en1,psi1] = dmrg(H,wfs,MPS(state),sweeps2,{"Silent=",silent,"Weight=",20.0});
    var = sqrt( abs( inner(H,psi1,H,psi1) - en1*en1 ) );
    maxDim = maxLinkDim(psi1);
    auto overlap = abs(inner(psi1,psi0));
    datafile << en1-en0 << " " << var << " " << maxDim << " " << overlap << " ";

    //
    // calculate second excited state
    //
    wfs.push_back(psi1);
    auto [en2,psi2] = dmrg(H,wfs,MPS(state),sweeps2,{"Silent=",silent,"Weight=",20.0});
    var = sqrt( abs( inner(H,psi2,H,psi2) - en2*en2 ) );
    maxDim = maxLinkDim(psi2);
    overlap = abs(inner(psi2,psi1));
    datafile << en2-en1 << " " << var << " " << maxDim << " " << overlap << " " << std::endl;

    printfln("Energy gaps = %.4f, %0.4f\n",en1-en0, en2-en1);

    for(int i=1; i<int(g.size()); i++)
        {
        printfln("g = %0.1f", g[i]);

        //
        // Create the Hamiltonian using AutoMPO
        //
        ampo = AutoMPO(sites);
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
            ampo += Jz, "Sz",b,"Sz",b+2;
        }
        for(int b = 1; b < N-1; b++){
            ampo += g[i],"Cdagup",b,"Ntot",b+1,"Cup",b+2;
            ampo += g[i],"Cdagup",b+2,"Ntot",b+1,"Cup",b;
            ampo += g[i],"Cdagdn",b,"Ntot",b+1,"Cdn",b+2;
            ampo += g[i],"Cdagdn",b+2,"Ntot",b+1,"Cdn",b;
        }
        H = toMPO(ampo);

        //
        // calculate ground state
        //
        en0 = dmrg(psi0, H,sweeps1,{"Silent",silent});
        var = sqrt( abs( inner(H,psi0,H,psi0) - en0*en0 ) );
        maxDim = maxLinkDim(psi0);
        datafile << g[i] << " " << en0 << " " << var << " " << maxDim << " " ;

        //
        // calculate first excited state
        //
        wfs = std::vector<MPS>(1);
        wfs.at(0) = psi0;
        en1 = dmrg(psi1,H,wfs,sweeps2,{"Silent=",silent,"Weight=",20.0});
        var = sqrt( abs( inner(H,psi1,H,psi1) - en1*en1 ) );
        maxDim = maxLinkDim(psi1);
        overlap = abs(inner(psi1,psi0));
        datafile << en1-en0 << " " << var << " " << maxDim << " " << overlap << " ";

        wfs.push_back(psi1);
        en2 = dmrg(psi2,H,wfs,sweeps2,{"Silent=",silent,"Weight=",20.0});
        var = sqrt( abs( inner(H,psi2,H,psi2) - en2*en2 ) );
        maxDim = maxLinkDim(psi2);
        overlap = abs(inner(psi2,psi1));
        datafile << en2-en1 << " " << var << " " << maxDim << " " << overlap << " " << std::endl;
            
        printfln("Energy gaps = %.4f, %0.4f\n",en1-en0, en2-en1);

    }

    datafile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}
