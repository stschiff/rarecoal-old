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
auto evalFileName = "/dev/null";

void mainFindsplit(string[] argv, Params_t params_) {
    p = params_;
    try {
        readParams(argv);
    }
    catch (Exception e) {
        printHelp(e);
        return;
    }
    auto model = new Model(input_data.nVec, p.popsizeVec, p.joins, p.migrations);
    auto stepper = Stepper.make_stepper(p.n0, p.lingen, p.tMax);
    auto eval_nr = to!int(max_eval_time / eval_dt);
    auto P = p.popsizeVec.length;
    auto minLogL = -double.infinity;
    auto minK = 0;
    auto minT = 0.0;
    auto evalFile = File(evalFileName, "w");
    foreach(k; 0 .. P) {
        if(k == newBranch)
            continue;
        for(auto t = eval_dt; t < max_eval_time; t += eval_dt) {
            // stderr.writefln("trying k=%s, t=%s", k, t);
            auto new_join = Join_t(t, to!int(k), newBranch, 1.0);
            auto new_model = new Model(input_data.nVec, p.popsizeVec, p.joins ~ new_join);
            try {
                auto logl_ = totalLikelihood(new_model, input_data, stepper, p.mu * 2.0 * p.n0);
                if(logl_ > minLogL) {
                    minK = to!int(k);
                    minT = t;
                    minLogL = logl_;
                }
                evalFile.writefln("%s\t%s\t%.2f", k, t, logl_);
                stderr.writefln("%s\t%s\t%.2f", k, t, logl_);
            }
            catch (IllegalModelException e) {
                stderr.writeln(e.msg);
                break;
            }
        }
    }
    writefln("branch:\t%s", minK);
    writefln("time:\t%s", minT);
    writefln("logL:\t%.2f", minLogL);
}

void readParams(string[] argv) {
    getopt(argv, std.getopt.config.caseSensitive,
           "max_af|m"       , &max_af,
           "newBranch|n"    , &newBranch,
           "eval_dt|d"      , &eval_dt,
           "max_eval_time|t", &max_eval_time,
           "evalFile|f"     , &evalFileName);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1], max_af);
    if(p.popsizeVec.length == 0) {
        p.popsizeVec = new double[input_data.nVec.length];
        p.popsizeVec[] = 1.0;
    }
    enforce(p.popsizeVec.length == input_data.nVec.length);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal findsplit [OPTIONS] <input_file>
Options:
    --max_af, -m                maximum allele frequency to use [10]
    --newBranch, -n             the new branch to test
    --eval_dt, -d               time interval for evaluation [0.0005]
    --max_eval_time, -t,        maximum time for evaluation [0.025]
    --evalFile, -f              file to write full evaluation to [/dev/null]
");
}
