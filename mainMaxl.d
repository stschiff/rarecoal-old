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
import params;

Data input_data;
bool fixedPopSize = false;
double theta;
int singleJoin = -1;
double maxTime = 1.0;
Params_t p;

void mainMaxl(string[] argv, Params_t params_) {
    p = params_;
    try {
        readParams(argv);
    }
    catch(Exception e) {
        printHelp(e);
        return;
    }
    theta = 2.0 * p.mu * p.n0;
    if(p.popsizeVec.length == 0) {
        p.popsizeVec = new double[input_data.nVec.length];
        p.popsizeVec[] = 1.0;
    }
    if(p.leaf_times.length == 0) {
        p.leaf_times = new double[input_data.nVec.length];
        p.leaf_times[] = 0.0;
    }
    enforce(p.popsizeVec.length == input_data.nVec.length);
    auto init_model = new Model(input_data.nVec, p.popsizeVec, p.joins, p.migrations, p.leaf_times);
    auto stepper = Stepper.make_stepper(p.n0, p.lingen, p.tMax);
    auto max_res = maximize(init_model, stepper, input_data, fixedPopSize, singleJoin);
    auto max_model = max_res[0];
    auto logl = max_res[1];
    stderr.writeln(max_model.joins);
    stderr.writeln(max_model.migrations);
    stderr.writeln(max_model.popSizeVec);
    report(max_model, logl);
}

void readParams(string[] argv) {
    getopt(argv, std.getopt.config.caseSensitive,
           "fixedPopSize|f", &fixedPopSize,
           "singleJoin"    , &singleJoin,
           "max_time"      , &maxTime);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1], p.max_af, p.indices, p.nrCalledBases);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal maxl [OPTIONS] <input_file>
Options:
    --singleJoin        only maximize single join
    --fixedPopSize, -f  keep population sizes fixed during maximization");
}

Tuple!(Model, double) maximize(Model init_model, Stepper stepper, Data input_data, bool fixedPopSize, int singleJoin) {
    
    auto minFunc = new MinFunc(init_model, input_data, stepper, fixedPopSize, theta, singleJoin,
                               p.minFreq, p.exclude_list, maxTime);
    auto powell = new Powell!MinFunc(minFunc);

    auto init_params = minFunc.model_to_params(init_model);
    auto min_params = powell.minimize(init_params);
    double logl = minFunc(min_params);

    auto min_model = minFunc.params_to_model(min_params);
    return tuple(min_model, logl);
}

void report(Model model, double logl) {
    writefln("Log Likelihood: %.2f", logl);
    writefln("Population sizes\t%s", model.popSizeVec.map!"text(a)"().join(","));
    foreach(j; model.joins)
        writefln("Join\t%s,%s,%s,%s", j.t, j.k, j.l, j.popsize);
    foreach(m; model.migrations)
        writefln("Migration\t%s,%s,%s", m.k, m.l, m.r);
    char[] cmdopt;
    cmdopt ~= format("-p %s", model.popSizeVec.map!"text(a)"().join(","));
    foreach(j; model.joins)
        cmdopt ~= format(" -j %s,%s,%s,%s", j.t, j.k, j.l, j.popsize);
    foreach(m; model.migrations)
        cmdopt ~= format(" -g %s,%s,%s", m.k, m.l, m.r);
    writefln("command line options\t%s", cmdopt);
}


