import std.math;
import coal_state;

class Stepper {
    size_t nr_steps;
    double tMax;
    double alpha;
    double[] time_boundaries;
    
    this(size_t nr_steps=10000, double tMax=20, double alpha=0.001) {
        this.nr_steps = nr_steps;
        this.alpha = alpha;
        this.tMax = tMax;
        time_boundaries = new double[nr_steps];
        foreach(i; 0 .. nr_steps)
            time_boundaries[i] = time_function(i);
    }
    
    double time_function(size_t i) const {
        return alpha * exp(cast(double)(i) / nr_steps * log(1.0 + tMax / alpha)) - alpha;
    }
    
    void run(CoalState state) const {
        foreach(i; 1 .. nr_steps) {
            auto next_t = time_boundaries[i];
            state.step(next_t);
        }
    }
}

unittest {
    import model;
    auto nVec = [500UL, 500, 500];
    auto joins = [Join_t(0.1, 0, 1, 1.0), Join_t(0.2, 0, 2, 1.0)];
    auto m = new Model(nVec, [1.0, 1.0, 1.0], joins);
    
    auto cs = new CoalState(m, [1, 2, 0]);
    
    auto stepper = new Stepper();
    
    stepper.run(cs);
    assert(cs.d > 0.0);
}