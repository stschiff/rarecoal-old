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
import params;

Data input_data;
bool fixedPopSize = false;
int max_af = 10;
double theta;
int burnin_cycles = 100;
int main_cycles = 1000;
int burnin_steps, main_steps;
auto tracefile = "/dev/null";
Params_t p;

void mainMCMC(string[] argv, Params_t params_) {
    try {
        readParams(argv);
    }
    catch(Exception e) {
        printHelp(e);
        return;
    }
    p = params_;
    theta = 2.0 * p.mu * p.n0;
    if(p.popsizeVec.length == 0) {
        p.popsizeVec = new double[input_data.nVec.length];
        p.popsizeVec[] = 1.0;
    }
    enforce(p.popsizeVec.length == input_data.nVec.length);

    auto init_model = new Model(input_data.nVec, p.popsizeVec, p.joins, p.migrations);
    auto stepper = Stepper.make_stepper(p.n0, p.lingen, p.tMax);
    auto minFunc = new MinFunc(init_model, input_data, stepper, fixedPopSize, theta, -1, p.minFreq, p.exclude_list);
    auto init_params = minFunc.model_to_params(init_model);
    burnin_steps = burnin_cycles * minFunc.nrParams;
    main_steps = main_cycles * minFunc.nrParams;
    auto mcmc = new MCMC!MinFunc(minFunc, init_params, tracefile);
    mcmc.run(burnin_steps + main_steps);
    reportMain(mcmc);
}

void readParams(string[] argv) {
    getopt(argv, std.getopt.config.caseSensitive,
           "fixedPopSize|f" , &fixedPopSize,
           "max_af|m"       , &max_af,
           "burnin_cycles|b", &burnin_cycles,
           "main_cycles|s"  , &main_cycles,
           "tracefile|t"    , &tracefile);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1], max_af);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal mcmc [OPTIONS] <input_file>
Options:
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
