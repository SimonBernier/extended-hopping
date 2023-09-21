//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

//calculate Sz Sz, SxSx+SySy and density-density correlators
std::tuple<double, double, double> correlators(int,int,MPS,SiteSet);

int main(int argc, char* argv[]){
    std::clock_t tStart = std::clock();

    //Parse the input file
    if(argc != 2) { printfln("Usage: %s inputfile_ext",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N",8);
    auto V = input.getReal("V",0.1);
    auto Jz = input.getReal("Jz",0.1);

    auto nsweeps = input.getInt("nsweeps");
    auto silent = input.getYesNo("silent",false);

    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
    //println(sweeps);

    auto table2 = InputGroup(input,"sweepsExcited");
    auto sweeps2 = Sweeps(nsweeps,table);

    // write results to file
    char schar[128];
    int n1 = std::sprintf(schar,"N_%d_V_%0.1f_Jz_%0.1f.dat",N,V,Jz);
    std::string s1(schar);
    std::ofstream datafile;
    datafile.open(s1); // opens the file
    if( !datafile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    datafile << "g" << " " << "en0" << " " << "var0" << " " << "maxDim0" << " " 
                << "gap01" << " " << "var1" << " " << "maxDim1" << " " << "overlap01" << " "
                << "gap12" << " " << "var2" << " " << "maxDim2" << " " << "overlap12" << " "
                << "gap23" << " " << "var3" << " " << "maxDim3" << " " << "overlap23" << " " 
                << "upd" << " " << "dnd" << " " << "F(r)" << " " << "G(r)" << " " << "Nd(r)" << " " << "Stot2" << std::endl; 

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = tJ(N);
    
    //
    // Set the initial wavefunction matrix product state
    //
    auto state = InitState(sites);
    int p = N/2;
    for(int i = N; i >= 1; --i) 
        {
        if(p > 0)
            {
            state.set(i,(i%2 == 0 ? "Dn" : "Up"));
            p -= 1;
            }
        else
            {
            state.set(i,"Emp");
            }
        }

    std::vector<double> g = {0., 0.25, 0.5, 0.75, 0.9, 1., 1.1, 1.25, 1.5, 1.75, 2.0};

    for(int i=0; i<int(g.size()); i++){
        printfln("\ng = %0.2f", g[i]);

        //
        // Create the Hamiltonian using AutoMPO
        //
        auto ampo = AutoMPO(sites);
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
        auto H = toMPO(ampo);

        //
        // calculate ground state
        //
        auto [en0,psi0] = dmrg(H,MPS(state),sweeps,{"Silent",silent});
        auto var = sqrt( abs( inner(H,psi0,H,psi0) - en0*en0 ) );
        auto maxDim = maxLinkDim(psi0);
        datafile << g[i] << " " << en0 << " " << var << " " << maxDim << " " ;

        //
        // calculate 3 first excited states
        //
        auto wfs = std::vector<MPS>(1);
        wfs.at(0) = psi0;
        auto [en1,psi1] = dmrg(H,wfs,psi0,sweeps2,{"Silent=",silent,"Weight=",20.0});
        var = inner(H,psi1,H,psi1) - en1*en1;
        maxDim = maxLinkDim(psi1);
        auto overlap = abs(inner(psi1,psi0));
        datafile << en1-en0 << " " << var << " " << maxDim << " " << overlap << " ";

        wfs.push_back(psi1);
        auto [en2,psi2] = dmrg(H,wfs,psi1,sweeps2,{"Silent=",silent,"Weight=",20.0});
        var = inner(H,psi2,H,psi2) - en2*en2;
        maxDim = maxLinkDim(psi2);
        overlap = abs(inner(psi2,psi1));
        datafile << en2-en1 << " " << var << " " << maxDim << " " << overlap << " ";

        wfs.push_back(psi2);
        auto [en3,psi3] = dmrg(H,wfs,psi2,sweeps2,{"Silent=",silent,"Weight=",20.0});
        var = inner(H,psi3,H,psi3) - en3*en3;
        maxDim = maxLinkDim(psi3);
        overlap = abs(inner(psi3,psi2));
        datafile << en3-en2 << " " << var << " " << maxDim << " " << overlap << " ";

        printfln("Energy gaps = %.4f, %0.4f, %0.4f\n",en1-en0, en2-en1, en3-en2);

        // single site correlator
        std::vector<double> upd(N,0.),dnd(N,0.);
        for(int j = 1; j <= N; j++){
            psi0.position(j);
            upd[j-1] = elt(dag(prime(psi0(j),"Site"))*op(sites,"Nup",j)*psi0(j));
            dnd[j-1] = elt(dag(prime(psi0(j),"Site"))*op(sites,"Ndn",j)*psi0(j));
        }

        println("Up Density:");
        for(int j = 0; j < N; j++){
            printf(" %.4f ",upd[j]);
            datafile << upd[j] << " ";
        }
        println();

        println("Dn Density:");
        for(int j = 0; j < N; j++){
            printf(" %.4f ",dnd[j]);
            datafile << dnd[j] << " ";
        }
        println();

        // two-site correlators
        std::vector<double> szsz(N,0.),sxysxy(N,0.),ndnd(N,0.);
        double Stot2 = 0.25*double(N);
        for(int j = 1; j <= N; j++){
            for(int k = 1; k <= N; k++){
                int r = abs(j-k);
                auto [zz,tt,dd] = correlators(j,k,psi0,sites);
                szsz[r] += zz;
                sxysxy[r] += tt;
                ndnd[r] += dd;
                if(r>0){
                    Stot2 += tt;
                }
            }
        }

        println("Sz-Sz Correlator:");
        for(int j = 0; j < N; j++){
            printf(" %.4f ",szsz[j]);
            datafile << szsz[j] << " ";
        }
        println();

        println("Transverse Spin Correlator:");
        for(int j = 0; j < N; j++){
            printf(" %.4f ",sxysxy[j]);
            datafile << sxysxy[j] << " ";
        }
        println();

        println("Density-Density Correlator:");
        for(int j = 0; j < N; j++){
            printf(" %.4f ",ndnd[j]);
            datafile << ndnd[j] << " ";
        }
        println();

        datafile << Stot2 << " " << std::endl;
            
    }// for g

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    datafile.close();

    return 0;
}// main

//calculate spin-spin correlator
std::tuple<double, double, double> correlators(int center, int b, MPS psi, SiteSet sites){
    
    double corrZ, corrT, corrD;

    psi.position(b);
    if(b>center){ //bring site b next to the center from right
        for(int n=b-1; n>center; n--){
            auto g = BondGate(sites,n,n+1);
            auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi.position(g.i1()); //move orthogonality center to the left 
        }
        auto ket = psi(center)*psi(center+1);
        auto SzSz = sites.op("Sz",center)*sites.op("Sz",center+1);
        corrZ = elt( dag(prime(ket,"Site")) * SzSz * ket);
        auto SpSm = 0.5*sites.op("S+",center)*sites.op("S-",center+1);
        auto SmSp = 0.5*sites.op("S-",center)*sites.op("S+",center+1);
        corrT = elt( dag(prime(ket,"Site")) * SpSm * ket);
        corrT += elt( dag(prime(ket,"Site")) * SmSp * ket);
        auto nn = sites.op("Ntot",center)*sites.op("Ntot",center+1);
        corrD = elt( dag(prime(ket,"Site")) * nn * ket);
    }
    else if(b<center){ //bring site b next to the center from left
        for(int n=b; n<center-1; n++){
          auto g = BondGate(sites,n,n+1);
          auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
          AA.replaceTags("Site,1","Site,0");
          psi.svdBond(g.i1(), AA, Fromleft); //svd from the right
          psi.position(g.i2()); //move orthogonality center to the right 
        }
        auto ket = psi(center-1)*psi(center);
        auto SzSz = sites.op("Sz",center-1)*sites.op("Sz",center);
        corrZ = elt( dag(prime(ket,"Site")) * SzSz * ket);
        auto SpSm = 0.5*sites.op("S+",center-1)*sites.op("S-",center);
        auto SmSp = 0.5*sites.op("S-",center-1)*sites.op("S+",center);
        corrT = elt( dag(prime(ket,"Site")) * SpSm * ket);
        corrT += elt( dag(prime(ket,"Site")) * SmSp * ket);
        auto nn = sites.op("Ntot",center-1)*sites.op("Ntot",center);
        corrD = elt( dag(prime(ket,"Site")) * nn * ket);
    }
    else{
        corrZ = 0.25;
        corrT = 0.5;
        corrD = 1.;
    }

    return {corrZ, corrT, corrD};

}//correlations
