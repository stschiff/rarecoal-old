import stepper;
import coal_state;
import model;
import std.getopt;
import std.string;
import std.algorithm;
import std.array;
import std.conv;
import std.stdio;
import std.range;
import std.math;
import std.exception;

int[] nVec, mVec;
double[] popsizeVec;
Join_t[] joins;
int nrSteps = 10000;
double alpha=0.001, tMax=20.0;

void mainProb(string[] argv) {
    
    try {
        readParams(argv);
    }
    catch(Exception e) {
        printHelp(e);
        return;
    }
    auto model = new Model(nVec, popsizeVec, joins);
    auto state = new CoalState(model, mVec);
    auto m = mVec.reduce!"a+b"();
    auto stepper = new Stepper(nrSteps, tMax, alpha);
    auto factor = zip(nVec, mVec).map!(x => binom(x[0], x[1])).reduce!"a*b"();
    stepper.run(state);
    auto result = m ? model.mu * factor * state.d : exp(-model.mu * state.e);
    writeln(result);

}

void readParams(string[] argv) {
    void handleJoins(string option, string str) {
        auto fields = str.split(",");
        auto t = fields[0].to!double();
        auto k = fields[1].to!int();
        auto l = fields[2].to!int();
        auto popsize = fields[3].to!double();
        joins ~= Join_t(t, k, l, popsize);
    }
    
    void handlePopsize(string option, string str) {
        popsizeVec = str.split(",").map!"a.to!double()"().array();
    }
    
    getopt(argv, std.getopt.config.caseSensitive,
           "nrSteps|N", &nrSteps,
           "alpha|a"  , &alpha,
           "Tmax|T"   , &tMax,
           "join|j"   , &handleJoins,
           "popsize|p" , &handlePopsize);

    enforce(argv.length == 3, "need more arguments");
    nVec = argv[1].split(",").map!"a.to!int()"().array();
    mVec = argv[2].split(",").map!"a.to!int()"().array();
    if(popsizeVec.length == 0) {
        popsizeVec = new double[nVec.length];
        popsizeVec[] = 1.0;
    }
    
    enforce(nVec.length == mVec.length && nVec.length == popsizeVec.length);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal prob [OPTIONS] <n1>[,<n2>[,...]] <m1>[,<m2>[,...]]
Options:
    --join, -j <t,k,l,p>        add a join at time t from population l to k, setting the new population size to p
    --nrSteps, -N <NR>          nr of steps in the stepper [10000]
    --alpha, -a <A>             time scale of transitioning from linear to log scale time intervals [0.001]
    --Tmax, -T <T>              maximum time interval boundary
    --popsize, -p <p1,p2,...>   initial population sizes");
}

double binom(int m, int k) {
  auto prod = 1.0;
  foreach(i; 1 .. k + 1) {
    prod *= to!double(m - (k - i)) / i;
  }
  return prod;
}

