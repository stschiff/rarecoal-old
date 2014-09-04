import std.conv;
import std.algorithm;
import std.range;
import std.math;
import std.parallelism;
import model;
import data;
import coal_state;
import stepper;

double totalLikelihood(in Model model, in Data input_dat, in Stepper stepper, double theta)
in {
    assert(model.nVec == input_dat.nVec);
}
body {
    auto total_l = 0.0;
    auto log_l = 0.0;
    auto m = input_dat.max_m;
    auto hom_type = (new int[model.nVec.length]).idup;
    // foreach(order; input_dat.standardOrder) {
    foreach(order; taskPool.parallel(input_dat.standardOrder)) {
        if(order == hom_type)
            continue;
        auto state = new CoalState(model, order);
        auto f = order.reduce!"a+b"();
        auto factor = zip(input_dat.nVec, order).map!(x => binom(x[0], x[1])).reduce!"a*b"();
        stepper.run(state);
        auto l = theta * state.d;
        assert(!isNaN(l) && l > 0.0, text(l, " ", state.a, " ", state.b, " ", model, " ", order));
        total_l += l * factor;
        log_l += log(l) * (order in input_dat.counts ? input_dat.counts[order] : 0.0);
    }
    auto higher_l = 1.0 - total_l;
    if(higher_l <= 0.0)
        throw new IllegalModelException("likelihood of configurations exceeds one");
    log_l += log(higher_l) * (input_dat.higher + input_dat.counts[hom_type.idup]);
    return log_l;
}

unittest {
    auto dat = new Data("testDat/testDat.txt");
    auto model = new Model([500UL, 500], [1.0, 1.0], [Join_t(0.4, 0, 1, 1.0)]);
    auto stepper = new Stepper();
    auto l = totalLikelihood(model, dat, stepper, 0.0005);
    assert(l < 0.0 && l > -double.infinity);
}

double binom(int m, int k) {
  return reduce!"a*b"(1.0, iota(1, k + 1).map!(i => to!double(m - (k - i)) / i));
}

unittest {
    assert(binom(4, 2) == 6);
}