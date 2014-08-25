import std.stdio;
import std.getopt;
import std.string;
import std.conv;
import std.algorithm;
import std.array;
import mainProb;
import mainLogl;
import mainMaxl;
import mainMCMC;
import model;
import params;

Params_t p;

version(unittest) {
    void main() {
        writeln("all tests passed");
    }
}
else {
    void main(string[] argv) {
        
        try {
            readParams(argv);
        }
        catch(Exception e) {
            writeln(e.msg);
            printHelp();
            return;
        }
        
        if(argv.length < 2) {
            stderr.writeln("need to specify subprogram");
            printHelp();
            return;
        }
        if(argv[1] == "prob") {
            mainProb.mainProb(argv[1 .. $], p);
        }
        else if(argv[1] == "logl") {
            mainLogl.mainLogl(argv[1 .. $], p);
        }
        else if(argv[1] == "maxl") {
            mainMaxl.mainMaxl(argv[1 .. $], p);
        }
        else if(argv[1] == "mcmc") {
            mainMCMC.mainMCMC(argv[1 .. $], p);
        }
        else {
            stderr.writeln("unknown subprogram: ", argv[1]);
            printHelp();
        }
    }
    
    void readParams(ref string[] argv) {
        void handleJoins(string option, string str) {
            auto fields = str.split(",");
            auto t = fields[0].to!double();
            auto k = fields[1].to!int();
            auto l = fields[2].to!int();
            auto popsize = 1.0;
            if(fields.length == 4)
                popsize = fields[3].to!double();
            p.joins ~= Join_t(t, k, l, popsize);
        }
    
        void handlePopsize(string option, string str) {
            p.popsizeVec = str.split(",").map!"a.to!double()"().array();
        }
        
        p = Params_t();
        p.mu = 1.25e-8;
        p.n0 = 20000;
        p.lingen = 400;
        p.tMax = 20.0;
            
        getopt(argv, std.getopt.config.passThrough,
            "mu"        , &p.mu,
            "n0"        , &p.n0,
            "lingen"    , &p.lingen,
            "tMax"      , &p.tMax,
            "join|j"    , &handleJoins,
            "popsize|p" , &handlePopsize);
    }
    
    void printHelp() {
        writeln("rarecoal [Options] <subprogram>
Options:
        --mu                          mutation rate per generation [1.25e-8]
        --n0                          effective population size, used for scaling time [20000]
        --lingen                      nr. of generations to propagate linearly, before crossing over
                                      to exponential intervals [400]
        --tMax                        maximum time to which to run, in 2n0 generations [20.0]
        --join, -j <t,k,l,popsize=1>  add a join at time t from population l to k, optionally setting the new popsize 
        --popsize, -p <p1,p2,...>     initial population sizes

Subprograms:        
        prob         compute the probability of a single configuration of derived alleles
        logl         compute the likelihood of a whole data set given parameters
        maxl         find maximum likelihood parameters for a given data set
        mcmc         Monte Carlo Markov Chain on parameters
");
    }
}