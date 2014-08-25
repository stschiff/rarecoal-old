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
import params;

Data input_data;
int max_af = 10;
string spectrumfile = "/dev/null";
Params_t p;

void mainLogl(string[] argv, Params_t params_) {
    
    try {
        readParams(argv);
    }
    catch(Exception e) {
        printHelp(e);
        return;
    }
    p = params_;
    auto model = new Model(input_data.nVec, p.popsizeVec, p.joins);
    auto stepper = Stepper.make_stepper(p.n0, p.lingen, p.tMax);
    auto logl = totalLikelihood(model, input_data, stepper, p.mu * 2.0 * p.n0);
    writeSpectrum(model, input_data, stepper, p.mu * 2.0 * p.n0, spectrumfile);
    writefln("%.2f", logl);

}

void readParams(string[] argv) {
    getopt(argv, std.getopt.config.caseSensitive,
           "max_af|m"  , &max_af,
           "spectrumfile|s", &spectrumfile);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1], max_af);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal prob [OPTIONS] <input_file>
Options:
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
