import std.conv;
import std.algorithm;
import std.array;
import std.string;
import stepper;
import data;
import model;
import logl;

class MinFunc {
    const Model init_model;
    const Data input_data;
    double penalty = 1.0e20;
    const Stepper stepper_;
    bool fixedPopSize;
    double theta;

    this(in Model init_model, in Data input_data, in Stepper stepper_, bool fixedPopSize, double theta) {
        this.init_model = init_model;
        this.input_data = input_data;
        this.stepper_ = stepper_;
        this.fixedPopSize = fixedPopSize;
        this.theta = theta;
        totalLikelihood(new Model(init_model.nVec, init_model.popSizeVec, init_model.joins), input_data, stepper_, theta); // this just serves to check for any exceptions with the initial model
    }
    
    double opCall(in double[] params)
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
    
    @property int nrParams() {
        return fixedPopSize ? to!int(init_model.joins.length) : 2 * to!int(init_model.joins.length) + init_model.P;
    }
    
    string[] paramNames() {
        string[] ret;
        foreach(j; init_model.joins)
            ret ~= format("t_join_%s_%s", j.k, j.l);
        if(!fixedPopSize) {
            foreach(i; 0 .. init_model.P)
                ret ~= format("N_%s", i);
            foreach(j; init_model.joins)
                ret ~= format("N_join_%s_%s", j.k, j.l);
        }
        return ret;
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
        auto popsizeVec = init_model.popSizeVec.dup;
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
            ret ~= model.popSizeVec.dup;
            ret ~= model.joins.map!(j => j.popsize.to!double())().array();
        }
        return ret;
    }
}
