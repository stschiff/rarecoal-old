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
    auto max_model = maximize(init_model, input_data, fixedPopSize);
    stderr.writeln(max_model.joins);
    stderr.writeln(max_model.popsizeVec);
    report(max_model);
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

Model maximize(Model init_model, Stepper stepper, Data input_data, bool fixedPopSize) {
    
    auto minFunc = new MinFunc(init_model, input_data, stepper, fixedPopSize);
    auto powell = new Powell!MinFunc(minFunc);

    auto init_params = minFunc.model_to_params(init_model);
    auto min_params = powell.minimize(init_params);

    auto min_model = minFunc.params_to_model(min_params);
    return min_model;
}

void report(Model model) {
    writefln("Population sizes\t%s", model.popsizeVec.map!"text(a)"().join(","));
    foreach(j; model.joins)
        writefln("Join\t%s,%s,%s,%s", j.t, j.k, j.l, j.popsize);
    char[] cmdopt;
    cmdopt ~= format("-p %s", model.popsizeVec.map!"text(a)"().join(","));
    foreach(j; model.joins)
        cmdopt ~= format(" -j %s,%s,%s,%s", j.t, j.k, j.l, j.popsize);
    writefln("command line options\t%s", cmdopt);
}

class MinFunc {
    const Model init_model;
    const Data input_data;
    double penalty = 1.0e20;
    const Stepper stepper_;
    bool fixedPopSize;

    this(in Model init_model, in Data input_data, in Stepper stepper_, bool fixedPopSize) {
        this.init_model = init_model;
        this.input_data = input_data;
        this.stepper_ = stepper_;
        this.fixedPopSize = fixedPopSize;
        totalLikelihood(new Model(init_model.nVec, init_model.popsizeVec, init_model.joins), input_data, stepper_); // this just serves to check for any exceptions with the initial model
    }
    
    double opCall(double[] params)
    {
        if(invalid(params))
            return penalty;
        double l;
        try {
            auto new_model = params_to_model(params);
            l = totalLikelihood(new_model, input_data, stepper_, theta);
        }
        catch(IllegalModelException e) {
            return penalty;
        }
        return -l;
    }
    
    bool invalid(in double[] params) {
        return params.any!"a<0.0"();
    }

    Model params_to_model(in double[] params)
    in {
        if(fixedPopSize)
            assert(params.length == init_model.joins.length);
        else
            assert(params.length == 2 * init_model.joins.length + init_model.P);
    }
    body {
        auto joins = init_model.joins.dup;
        auto nj = joins.length;
        auto popsizeVec = init_model.popsizeVec.dup;
        auto K = init_model.P;
        foreach(i, p; params[0 .. nj])
            joins[i].t = p;
        
        if(!fixedPopSize) {
            foreach(i, p; params[nj .. nj + K])
                popsizeVec[i] = p;
            foreach(i, p; params[nj + K .. $])
                joins[i].popsize = p;
        }
        return new Model(init_model.nVec, popsizeVec, joins);
    }

    double[] model_to_params(in Model model) {
        auto ret = model.joins.map!(j => j.t.to!double())().array();
        if(!fixedPopSize) {
            ret ~= model.popsizeVec.dup;
            ret ~= model.joins.map!(j => j.popsize.to!double())().array();
        }
        return ret;
    }
}


