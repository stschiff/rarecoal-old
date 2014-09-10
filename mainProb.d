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
import params;

int[] nVec, mVec;
Params_t p;

void mainProb(string[] argv, Params_t params_) {
    
    try {
        readParams(argv);
    }
    catch(Exception e) {
        printHelp(e);
        return;
    }
    p = params_;
    auto model = new Model(nVec, p.popsizeVec, p.joins, p.migrations);
    auto state = new CoalState(model, mVec);
    auto m = mVec.reduce!"a+b"();
    auto stepper = Stepper.make_stepper(p.n0, p.lingen, p.tMax);
    auto factor = zip(nVec, mVec).map!(x => binom(x[0], x[1])).reduce!"a*b"();
    stepper.run(state);
    auto theta = 2.0 * p.n0 * p.mu;
    auto result = theta * factor * state.d;
    writeln(result);

}

void readParams(string[] argv) {
    enforce(argv.length == 3, "need more arguments");
    nVec = argv[1].split(",").map!"a.to!int()"().array();
    mVec = argv[2].split(",").map!"a.to!int()"().array();
    if(p.popsizeVec.length == 0) {
        p.popsizeVec = new double[nVec.length];
        p.popsizeVec[] = 1.0;
    }
    enforce(nVec.length == mVec.length && nVec.length == p.popsizeVec.length);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal prob <n1>[,<n2>[,...]] <m1>[,<m2>[,...]]
");
}

double binom(int m, int k) {
  auto prod = 1.0;
  foreach(i; 1 .. k + 1) {
    prod *= to!double(m - (k - i)) / i;
  }
  return prod;
}

