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
    const double[] initialPopSizeVec;
    double[] popsizeVec;
    const Join_t[] joins;
    int joins_index;
    
    this(in int[] nVec, in double[] popsizeVec, in Join_t[] joins=[]) {
        this.nVec = nVec;
        if(popsizeVec.any!(p => p < 0.001)() || joins.any!(j => j.popsize < 0.001)())
            throw new IllegalModelException("population size can't be lower than 0.001");
        this.initialPopSizeVec = popsizeVec;
        assert(popsizeVec.all!"a>0.0"());
        this.popsizeVec = popsizeVec.dup;
        assert(joins.all!"a.t>0.0 && a.popsize>0.0"());
        this.joins = joins.dup.sort!"a.t < b.t"().array();
        this.joins_index = 0;
    }
    
    @property int P() const {
        return nVec.length.to!int();
    }
    
    override string toString() const {
        return format("Model {popsizeVec=%s, joins=%s}", popsizeVec, joins);
    }
    
    double coal_rate(int k) {
        return 1.0 / popsizeVec[k];
    }
    
    Join_t next_join() {
        if(joins_index >= joins.length) {
            auto ret = Join_t();
            ret.t = double.infinity;
            return ret;
        }
        return joins[joins_index];
    }
    
    void join_done() {
        auto j = joins[joins_index];
        popsizeVec[j.k] = j.popsize;
        joins_index += 1;
    }
    
    void reset() {
        joins_index = 0;
        popsizeVec = initialPopSizeVec.dup;
    }
}

unittest {
    auto nVec = [500, 500, 500];
    auto joins = [Join_t(0.1, 0, 1, 0.5), Join_t(0.2, 0, 2, 2.0)];
    auto m = new Model(nVec, [1.0, 1.0, 1.0], joins);
    assert(m.P == 3);
    assert(m.next_join() == Join_t(0.1, 0, 1, 0.5));
    m.join_done();
    assert(m.popsizeVec[0] == 0.5);
    assert(m.next_join() == Join_t(0.2, 0, 2, 2.0));
    m.join_done();
    assert(m.popsizeVec[0] == 2.0);
    assert(m.next_join().t == double.infinity);
}