import std.getopt;
import std.string;
import std.algorithm;
import std.array;
import std.conv;
import std.stdio;
import std.range;
import std.math;
import std.exception;
import stepper;
import coal_state;
import model;
import data;
import logl;
import powell;
import minfunc;
import mcmc;

Data input_data;
Join_t[] joins;
double[] popsizeVec;
bool fixedPopSize = false;
int max_af = 10;
double theta;
int burnin_cycles = 100;
int main_cycles = 1000;
int burnin_steps, main_steps;
auto tracefile = "/dev/null";

void mainMCMC(string[] argv, double mu, int n0, int lingen, double tMax) {
    try {
        readParams(argv);
    }
    catch(Exception e) {
        printHelp(e);
        return;
    }
    theta = 2.0 * mu * n0;
    auto init_model = new Model(input_data.nVec, popsizeVec, joins);
    auto stepper = Stepper.make_stepper(n0, lingen, tMax);
    auto minFunc = new MinFunc(init_model, input_data, stepper, fixedPopSize, theta);
    auto init_params = minFunc.model_to_params(init_model);
    burnin_steps = burnin_cycles * minFunc.nrParams;
    main_steps = main_cycles * minFunc.nrParams;
    auto mcmc = new MCMC!MinFunc(minFunc, init_params);
    mcmc.run(burnin_steps + main_steps);
    reportMain(mcmc);
    reportTraces(mcmc, tracefile);
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
           "join|j"         , &handleJoins,
           "popsize|p"      , &handlePopsize,
           "fixedPopSize|f" , &fixedPopSize,
           "max_af|m"       , &max_af,
           "burnin_cycles|b", &burnin_cycles,
           "main_cycles|s"  , &main_cycles,
           "tracefile|t"    , &tracefile);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1], max_af);
    if(popsizeVec.length == 0) {
        popsizeVec = new double[input_data.nVec.length];
        popsizeVec[] = 1.0;
    }
    enforce(popsizeVec.length == input_data.nVec.length);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal mcmc [OPTIONS] <input_file>
Options:
    --join, -j <t,k,l,p>        add a join at time t from population l to k, setting new popsize to p
    --popsize, -p <p1,p2,...>   initial population sizes
    --fixedPopSize, -f          keep population sizes fixed during maximization
    --max_af, -m                maximum allele frequency to use [10]
    --burnin_steps, -b          MCMC burnin steps [100]
    --main_steps, -s            MCMC main steps [1000]
    --tracefile, -t             File to write the MCMC trace to [/dev/null]");
}

void reportMain(MCMC!MinFunc mcmc) {
    writeln("Param\tMaxL\tMedian\tCI_lower\tCI_higher");
    auto names = mcmc.minfunc.paramNames();
    auto minIndex = main_steps - mcmc.fvalues[burnin_steps .. $].minPos().length.to!int();
    // stderr.writeln(burnin_steps, " ", main_steps, " ", minIndex);
    assert(minIndex >= 0 && minIndex < main_steps);
    foreach(i; 0 .. mcmc.N) {
        write(names[i], "\t");
        auto singleTrace = mcmc.trace[burnin_steps .. $].map!(t => t[i]).array();
        // stderr.writeln(singleTrace.length);
        write(singleTrace[minIndex], "\t");
        auto median = singleTrace.sort().drop(main_steps / 2).front();
        auto ci_lower = singleTrace.sort().drop(to!int(main_steps * 0.025)).front();
        auto ci_higher = singleTrace.sort().dropBack(to!int(main_steps * 0.025)).back();
        writefln("%s\t%s\t%s", median, ci_lower, ci_higher);
    }
}

void reportTraces(MCMC!MinFunc mcmc, string tracefile) {
    auto f = File(tracefile, "w");
    f.write(mcmc.minfunc.paramNames().join("\t"));
    f.writeln("\tscore");
    foreach(i, point; mcmc.trace) {
        f.write(point.map!"text(a)"().join("\t"));
        f.writefln("\t%20.2f", mcmc.fvalues[i]);
    }
    f.close();
}