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
    int singleJoin;
    int minFreq;
    double max_time;
    const int[][] exclude_list;

    this(in Model init_model, in Data input_data, in Stepper stepper_, bool fixedPopSize, double theta, int singleJoin,
         int minFreq, in int[][] exclude_list, double max_time) {
        this.init_model = init_model;
        this.input_data = input_data;
        this.stepper_ = stepper_;
        this.fixedPopSize = fixedPopSize;
        this.theta = theta;
        this.singleJoin = singleJoin;
        this.exclude_list = exclude_list;
        this.minFreq = minFreq;
        this.max_time = max_time;
        // totalLikelihood(new Model(init_model.nVec, init_model.popSizeVec, init_model.joins), input_data, stepper_, theta); // this just serves to check for any exceptions with the initial model
        totalLikelihood(init_model, input_data, stepper_, theta, minFreq, exclude_list);
    }
    
    double opCall(in double[] params)
    {
        if(invalid(params))
            return penalty;
        double l;
        try {
            auto new_model = params_to_model(params);
            l = totalLikelihood(new_model, input_data, stepper_, theta, minFreq, exclude_list);
        }
        catch(IllegalModelException e) {
            return penalty;
        }
        return -l;
    }
    
    @property int nrParams() {
        if(singleJoin >= 0) {
            return fixedPopSize ? 1 : 3;
        }
        else {
            auto nrJoins = init_model.joins.length.to!int();
            auto nrMigs = init_model.migrations.length.to!int();
            return fixedPopSize ? nrJoins + nrMigs : 2 * nrJoins + nrMigs + init_model.P;
        }
    }
    
    string[] paramNames() {
        string[] ret;
        if(singleJoin >= 0) {
            auto j = init_model.joins[singleJoin];
            ret ~= format("t_join_%s_%s", j.k, j.l);
            if(!fixedPopSize) {
                ret ~= format("N_%s", j.l);
                ret ~= format("N_join_%s_%s", j.k, j.l);
            }
        }
        else {
            foreach(j; init_model.joins)
                ret ~= format("t_join_%s_%s", j.k, j.l);
            foreach(m; init_model.migrations)
                ret ~= format("mig_%s_%s", m.k, m.l);
            if(!fixedPopSize) {
                foreach(i; 0 .. init_model.P)
                    ret ~= format("N_%s", i);
                foreach(j; init_model.joins)
                    ret ~= format("N_join_%s_%s", j.k, j.l);
            }
        }
        return ret;
    }
    
    bool invalid(in double[] params) {
        auto neg = params.any!"a<0.0"();
        auto max_t = singleJoin >= 0 ? params[0] > max_time : params[0..init_model.joins.length].any!(t => t > max_time)();
        return neg || max_t;
    }

    Model params_to_model(in double[] params)
    in {
        assert(params.length == nrParams);
    }
    body {
        auto joins = init_model.joins.dup;
        auto popsizeVec = init_model.popSizeVec.dup;
        auto migrations = init_model.migrations.dup;
        if(singleJoin >= 0) {
            joins[singleJoin].t = params[0];
            if(!fixedPopSize) {
                joins[singleJoin].popsize = params[1];
                popsizeVec[joins[singleJoin].l] = params[2];
            }
        }
        else {
            auto nj = joins.length;
            auto nm = migrations.length;
            auto K = init_model.P;
            foreach(i, p; params[0 .. nj])
                joins[i].t = p;
            foreach(i, p; params[nj .. nj + nm])
                migrations[i].r = p;
            if(!fixedPopSize) {
                foreach(i, p; params[nj + nm .. nj + nm + K])
                    popsizeVec[i] = p;
                foreach(i, p; params[nj + nm + K .. $])
                    joins[i].popsize = p;
            }
        }
        return new Model(init_model.nVec, popsizeVec, joins, migrations, init_model.leaf_times);
    }

    double[] model_to_params(in Model model) {
        double[] ret;
        if(singleJoin >= 0) {
            ret ~= model.joins[singleJoin].t;
            if(!fixedPopSize) {
                ret ~= model.joins[singleJoin].popsize;
                ret ~= model.popSizeVec[model.joins[singleJoin].l];
            }
        }
        else {
            ret = model.joins.map!(j => j.t.to!double())().array();
            ret ~= model.migrations.map!(m => m.r.to!double())().array();
            if(!fixedPopSize) {
                ret ~= model.popSizeVec.dup;
                ret ~= model.joins.map!(j => j.popsize.to!double())().array();
            }
        }
        return ret;
    }
}
