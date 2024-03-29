import std.stdio;
import std.getopt;
import std.string;
import std.conv;
import std.algorithm;
import std.parallelism;
import std.array;
import mainProb;
import mainLogl;
import mainMaxl;
import mainMCMC;
import mainFindsplit;
import model;
import params;

Params_t p;
uint nrThreads;

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
        reportParams();
        
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
        else if(argv[1] == "findsplit") {
            mainFindsplit.mainFindsplit(argv[1 .. $], p);
        }
        else {
            stderr.writeln("unknown subprogram: ", argv[1]);
            printHelp();
        }
    }
}
    
void readParams(ref string[] argv) {
    void handleMigrations(string option, string str) {
        auto fields = str.split(",");
        auto k = fields[0].to!int();
        auto l = fields[1].to!int();
        auto r = fields[2].to!double();
        p.migrations ~= Migration_t(k, l, r);
    }

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

    void handleLeafTimes(string option, string str) {
        p.leaf_times = str.split(",").map!"a.to!double()"().array();
    }

    void handlePopsize(string option, string str) {
        p.popsizeVec = str.split(",").map!"a.to!double()"().array();
    }

    void handleIndices(string option, string str) {
        p.indices = str.split(",").map!"a.to!int()"().array();
    }
    
    void handleExcludes(string option, string str) {
        p.exclude_list = str.split(";").map!(w => w.split(",").map!"a.to!int()"().array()).array();
    }
    
    p = Params_t();
    p.mu = 1.25e-8;
    p.n0 = 20000;
    p.lingen = 400;
    p.tMax = 20.0;
    p.minFreq = 1;
    p.max_af = 10;
    p.nrCalledBases = 2064554803;
        
    getopt(argv, std.getopt.config.passThrough,
        "mu"           , &p.mu,
        "n0"           , &p.n0,
        "lingen"       , &p.lingen,
        "tMax"         , &p.tMax,
        "join|j"       , &handleJoins,
        "migration|g"  , &handleMigrations,
        "leaf_times"   , &handleLeafTimes,
        "popsize|p"    , &handlePopsize,
        "minFreq"      , &p.minFreq,
        "exclude"      , &handleExcludes,
        "nrThreads"    , &nrThreads,
        "indices"      , &handleIndices,
        "max_af|m"     , &p.max_af,
        "nrCalledBases", &p.nrCalledBases);

    if(nrThreads)
      std.parallelism.defaultPoolThreads(nrThreads);

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
    --migration, -g <k, l, rate>  add a bi-directional migration rate between branch k and l (up to the first coalescence)
    --leaf_times <t1,t2,t3...>    for ancient samples choose t > 0.                   
    --popsize, -p <p1,p2,...>     initial population sizes
    --minFreq                     minimum frequency to evaluate
    --exclude                     exclude these patterns (semicolon-separated patterns of comma-separated patters)
    --nrThreads                   nr of Threads to use, default: nr of CPUs
    --indices i1,i2,...           If given, use only those populations in the input file specified by the indices.
    --max_af, -m                  Maximum allele frequency to read in.
    --nrCalledBases               nr of called bases in the genome [default: 2064554803]

Subprograms:        
    prob         compute the probability of a single configuration of derived alleles
    logl         compute the likelihood of a whole data set given parameters
    maxl         find maximum likelihood parameters for a given data set
    mcmc         Monte Carlo Markov Chain on parameters
    findsplit    Systematically evaluate likelihood of adding a new branch anywhere on an existing tree
");
}

void reportParams() {
    stderr.writeln("mu:              ", p.mu);
    stderr.writeln("n0:              ", p.n0);
    stderr.writeln("lingen:          ", p.lingen);
    stderr.writeln("tMax:            ", p.tMax);
    stderr.writeln("minFreq:         ", p.minFreq);
    stderr.writeln("exclude:         ", p.exclude_list);
    stderr.writeln("joins:           ", p.joins);
    stderr.writeln("migrations:      ", p.migrations);
    stderr.writeln("population sizes:", p.popsizeVec);
    stderr.writeln("leaf times:      ", p.leaf_times);
    stderr.writeln("indices:         ", p.indices);
    stderr.writeln("max_af:          ", p.max_af);
    stderr.writeln("nrCalledBases:   ", p.nrCalledBases);
}