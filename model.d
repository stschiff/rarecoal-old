import std.typecons;
import std.string;
import std.algorithm;
import std.conv;
import std.array;

alias Tuple!(double, "t", int, "k", int, "l", double, "popsize") Join_t;

class IllegalModelException : Exception {
    this(string msg) {
        super(msg);
    }
}

class Model {
    const int[] nVec;
    const double[] popSizeVec;
    const Join_t[] joins;
    
    this(in int[] nVec, in double[] popsizeVec, in Join_t[] joins=[]) {
        this.nVec = nVec;
        if(popsizeVec.any!(p => p < 0.001)() || joins.any!(j => j.popsize < 0.001)())
            throw new IllegalModelException("population size can't be lower than 0.001");
        if(popsizeVec.any!(p => p > 100.0)() || joins.any!(j => j.popsize > 100.0)())
            throw new IllegalModelException("population size can't be greater than 100");
        this.popSizeVec = popsizeVec;
        assert(popsizeVec.all!"a>0.0"());
        assert(joins.all!"a.t>0.0 && a.popsize>0.0"());
        assert(joins.all!"a.k!=a.l"());
        this.joins = joins;
    }
    
    @property int P() const {
        return nVec.length.to!int();
    }
    
    override string toString() const {
        return format("Model {popsizeVec=%s, joins=%s}", popSizeVec, joins);
    }
    
}
 
class ModelState {
    const Model model;
    const Join_t[] sorted_joins;
    double[] popsizeVec;
    int joins_index;
    
    this(in Model model) {
        this.model = model;
        this.sorted_joins = model.joins.dup.sort!"a.t < b.t"().array();
        this.popsizeVec = model.popSizeVec.dup;
        this.joins_index = 0;
    }
    
    Join_t next_join() {
        if(joins_index >= sorted_joins.length) {
            auto ret = Join_t();
            ret.t = double.infinity;
            return ret;
        }
        return sorted_joins[joins_index];
    }
    
    void join_done() {
        auto j = sorted_joins[joins_index];
        popsizeVec[j.k] = j.popsize;
        joins_index += 1;
    }
    
    double coal_rate(int k) {
        return 1.0 / popsizeVec[k];
    }
}

unittest {
    auto nVec = [500, 500, 500];
    auto joins = [Join_t(0.1, 0, 1, 0.5), Join_t(0.2, 0, 2, 2.0)];
    auto m = new Model(nVec, [1.0, 1.0, 1.0], joins);
    assert(m.P == 3);
    auto ms = new ModelState(m);
    assert(ms.next_join() == Join_t(0.1, 0, 1, 0.5));
    ms.join_done();
    assert(ms.popsizeVec[0] == 0.5);
    assert(ms.next_join() == Join_t(0.2, 0, 2, 2.0));
    ms.join_done();
    assert(ms.popsizeVec[0] == 2.0);
    assert(ms.next_join().t == double.infinity);
}