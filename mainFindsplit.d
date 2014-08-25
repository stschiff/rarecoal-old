import std.stdio;
import std.getopt;
import std.conv;
import std.exception;
import params;
import data;
import model;
import stepper;
import logl;


int newBranch;
Params_t p;
int max_af = 10;
Data input_data;
double eval_dt = 0.0005;
double max_eval_time = 0.025;

void mainFindsplit(string[] argv, Params_t params_) {
    p = params_;
    try {
        readParams(argv);
    }
    catch (Exception e) {
        printHelp(e);
        return;
    }
    auto model = new Model(input_data.nVec, p.popsizeVec, p.joins);
    auto stepper = Stepper.make_stepper(p.n0, p.lingen, p.tMax);
    auto logl = totalLikelihood(model, input_data, stepper, p.mu * 2.0 * p.n0);
    auto eval_nr = to!int(max_eval_time / eval_dt);
    auto P = p.popsizeVec.length;
    foreach(k; 0 .. P) {
        for(auto t = eval_dt; t < max_eval_time; t += eval_dt) {
            auto new_join = Join_t(t, to!int(k), newBranch, 1.0);
            auto new_model = new Model(input_data.nVec, p.popsizeVec, p.joins ~ new_join);
            auto logl_ = 0.0;
            try {
                logl_ = totalLikelihood(new_model, input_data, stepper, p.mu * 2.0 * p.n0);
            }
            catch (IllegalModelException e) {
                break;
            }
            writefln("%s\t%s\t%s", k, t, logl_);
        }
    }
}

void readParams(string[] argv) {
    getopt(argv, std.getopt.config.caseSensitive,
           "max_af|m"       , &max_af,
           "newBranch|n"    , &newBranch,
           "eval_dt|d"      , &eval_dt,
           "max_eval_time|t", &max_eval_time);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1], max_af);
    enforce(p.popsizeVec.length == input_data.nVec.length);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal findsplit [OPTIONS] <input_file>
Options:
    --max_af, -m                maximum allele frequency to use [10]
    --newBranch, -n             the new branch to test
    --eval_dt, -d               time interval for evaluation
    --max_eval_time, -t,        maximum time for evaluation
");
}
