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

Data input_data;
Join_t[] joins;
double[] popsizeVec;
int max_af = 10;
string spectrumfile = "/dev/null";

void mainLogl(string[] argv, double mu, int n0, int lingen, double tMax) {
    
    try {
        readParams(argv);
    }
    catch(Exception e) {
        printHelp(e);
        return;
    }
    auto model = new Model(input_data.nVec, popsizeVec, joins);
    auto stepper = Stepper.make_stepper(n0, lingen, tMax);
    auto logl = totalLikelihood(model, input_data, stepper, mu * 2.0 * n0);
    writeSpectrum(model, input_data, stepper, mu * 2.0 * n0, spectrumfile);
    writefln("%.2f", logl);

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
           "join|j"    , &handleJoins,
           "popsize|p" , &handlePopsize,
           "max_af|m"  , &max_af,
           "spectrumfile|s", &spectrumfile);
    
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
    writeln("./rarecoal prob [OPTIONS] <input_file>
Options:
    --join, -j <t,k,l,popsize>  add a join at time t from population l to k, setting the new popsize
    --popsize, -p <p1,p2,...>   initial population sizes
    --max_af, -m                maximum allele frequency to use [10]");
}

void writeSpectrum(Model model, in Data input_dat, in Stepper stepper, double theta, string spectrumfile)
in {
    assert(model.nVec == input_dat.nVec);
}
body {
    auto f = File(spectrumfile, "w");
    auto L = input_dat.counts.values.reduce!"a+b"() + input_dat.higher;
    foreach(key, count; input_dat.counts) {
        if(key.reduce!"a+b"() == 0)
            continue;
        model.reset();
        auto state = new CoalState(model, key);
        auto factor = zip(input_dat.nVec, key).map!(x => binom(x[0], x[1])).reduce!"a*b"();
        stepper.run(state);
        auto pred = theta * state.d * factor * L;
        f.writefln("%s %s %s", key.map!"text(a)"().join(","), count, pred);
    }
}
