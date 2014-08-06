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
size_t nrSteps = 10000;
double alpha=0.001, tMax=20.0;
Stepper stepper_;

void mainMaxl(string[] argv) {
    try {
        readParams(argv);
    }
    catch(Exception e) {
        printHelp(e);
        return;
    }
    auto init_model = new Model(input_data.nVec, joins);
    auto max_model = maximize(init_model, input_data);
    writeln(max_model.joins);
}

void readParams(string[] argv) {
    void handleJoins(string option, string str) {
        auto fields = str.split(",");
        auto t = fields[0].to!double();
        auto k = fields[1].to!size_t();
        auto l = fields[2].to!size_t();
        joins ~= Join_t(t, k, l);
    }
    
    getopt(argv, std.getopt.config.caseSensitive,
           "nrSteps|N", &nrSteps,
           "alpha|a"  , &alpha,
           "Tmax|T"   , &tMax,
           "join|j"   , &handleJoins);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1]);
    stepper_ = new Stepper(nrSteps, tMax, alpha);
}

void printHelp(Exception e) {
    writeln(e.msg);
    writeln("./rarecoal maxl [OPTIONS] <input_file>
Options:
    --join, -j <t,k,l>   add a join at time t from population l to k
    --nrSteps, -N <NR>  nr of steps in the stepper [10000]
    --alpha, -a <A>     time scale of transitioning from linear to log scale time intervals [0.001]
    --Tmax, -T <T>      maximum time interval boundary");
}

Model maximize(Model init_model, Data input_data) {
    
    auto minFunc = new MinFunc(init_model, input_data, stepper_);
    auto powell = new Powell!MinFunc(minFunc);

    auto init_params = model_to_params(init_model);
    auto min_params = powell.minimize(init_params);

    auto min_model = params_to_model(min_params, init_model);
    return min_model;
}

class MinFunc {
    const Model model;
    const Data input_data;
    double penalty = 1.0e20;
    const Stepper stepper_;

    this(in Model model, in Data input_data, in Stepper stepper_) {
        this.model = model;
        this.input_data = input_data;
        this.stepper_ = stepper_;
    }
    
    double opCall(double[] params)
    in {
        assert(params.length == model.joins.length);
    }
    body {
        if(invalid(params))
            return penalty;
        auto new_model = params_to_model(params, model);
        auto l = totalLikelihood(new_model, input_data, stepper_);
        // writeln(l);
        return -l;
    }
    
    bool invalid(in double[] params) {
        return params.any!"a<0.0"();
    }
}

double[] model_to_params(in Model model) {
    return model.joins.map!(j => j.t)().array().dup;
}

Model params_to_model(in double[] params, in Model init_model) {
    auto joins = init_model.joins.dup;
    foreach(i, p; params)
        joins[i].t = p;
    return new Model(init_model.nVec, joins);
}