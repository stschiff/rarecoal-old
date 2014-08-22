import std.getopt;
import std.string;
import std.algorithm;
import std.array;
import std.conv;
import std.stdio;
import std.range;
import std.math;
import std.exception;
import std.typecons;
import stepper;
import coal_state;
import model;
import data;
import logl;
import powell;
import minfunc;

Data input_data;
Join_t[] joins;
double[] popsizeVec;
bool fixedPopSize = false;
int max_af = 10;
double theta;

void mainMaxl(string[] argv, double mu, int n0, int lingen, double tMax) {
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
    auto max_res = maximize(init_model, stepper, input_data, fixedPopSize);
    auto max_model = max_res[0];
    auto logl = max_res[1];
    stderr.writeln(max_model.joins);
    stderr.writeln(max_model.popsizeVec);
    report(max_model, logl);
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
           "join|j"        , &handleJoins,
           "popsize|p"     , &handlePopsize,
           "fixedPopSize|f", &fixedPopSize,
           "max_af|m"      , &max_af);
    
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
    writeln("./rarecoal maxl [OPTIONS] <input_file>
Options:
    --join, -j <t,k,l,p>   add a join at time t from population l to k, setting new popsize to p
    --popsize, -p <p1,p2,...>   initial population sizes
    --fixedPopSize, -f  keep population sizes fixed during maximization
    --max_af, -m        maximum allele frequency to use [10]");

}

Tuple!(Model, double) maximize(Model init_model, Stepper stepper, Data input_data, bool fixedPopSize) {
    
    auto minFunc = new MinFunc(init_model, input_data, stepper, fixedPopSize, theta);
    auto powell = new Powell!MinFunc(minFunc);

    auto init_params = minFunc.model_to_params(init_model);
    auto min_params = powell.minimize(init_params);
    double logl = minFunc(min_params);

    auto min_model = minFunc.params_to_model(min_params);
    return tuple(min_model, logl);
}

void report(Model model, double logl) {
    writefln("Log Likelihood: %.2f", logl);
    writefln("Population sizes\t%s", model.popsizeVec.map!"text(a)"().join(","));
    foreach(j; model.joins)
        writefln("Join\t%s,%s,%s,%s", j.t, j.k, j.l, j.popsize);
    char[] cmdopt;
    cmdopt ~= format("-p %s", model.popsizeVec.map!"text(a)"().join(","));
    foreach(j; model.joins)
        cmdopt ~= format(" -j %s,%s,%s,%s", j.t, j.k, j.l, j.popsize);
    writefln("command line options\t%s", cmdopt);
}


