import std.conv;
import std.algorithm;
import std.range;
import std.math;
import model;
import data;
import coal_state;
import stepper;

double totalLikelihood(Model model, in Data input_dat, in Stepper stepper, size_t nrSteps=10000, double alpha=0.001, double tMax=20.0)
in {
    assert(model.nVec == input_dat.nVec);
}
body {
    auto total_l = 0.0;
    auto log_l = 0.0;
    auto m = input_dat.max_m;
    foreach(order; input_dat.standardOrder) {
        model.reset();
        auto state = new CoalState(model, order);
        auto f = order.reduce!"a+b"();
        auto factor = zip(input_dat.nVec, order).map!(x => binom(x[0], x[1])).reduce!"a*b"();
        stepper.run(state);
        auto l = f ? model.mu * state.d : exp(-model.mu * state.e);
        // stderr.writeln(order, " ", l);
        total_l += l * factor;
        log_l += log(l) * (order in input_dat.counts ? input_dat.counts[order] : 0.0);
    }
    auto higher_l = 1.0 - total_l;
    // stderr.writeln("HIGHER ", higher_l);
    assert(higher_l > 0.0);
    log_l += log(higher_l) * input_dat.higher;
    return log_l;
}

unittest {
    auto dat = new Data("testDat/testDat.txt");
    auto model = new Model([500UL, 500], [1.0, 1.0], [Join_t(0.4, 0, 1, 1.0)]);
    auto stepper = new Stepper();
    auto l = totalLikelihood(model, dat, stepper);
    assert(l < 0.0 && l > -double.infinity);
}

double binom(size_t m, size_t k) {
  return reduce!"a*b"(1.0, iota(1, k + 1).map!(i => to!double(m - (k - i)) / i));
}

unittest {
    assert(binom(4, 2) == 6);
}