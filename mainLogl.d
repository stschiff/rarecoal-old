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
    writeln(logl);

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
           "join|j"   , &handleJoins,
           "popsize|p", &handlePopsize);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1]);
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
    --join, -j <t,k,l,popsize>   add a join at time t from population l to k, setting the new popsize
    --popsize, -p <p1,p2,...>   initial population sizes");
}

