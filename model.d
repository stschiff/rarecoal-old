import std.typecons;

alias Tuple!(double, "t", size_t, "k", size_t, "l", double, "popsize") Join_t;

class Model {
    const size_t[] nVec;
    const double[] initialPopSizeVec;
    double[] popsizeVec;
    const Join_t[] joins;
    double mu;
    size_t joins_index;
    
    this(in size_t[] nVec, in double[] popsizeVec, in Join_t[] joins=[], double mu=0.0005) {
        this.nVec = nVec;
        this.initialPopSizeVec = popsizeVec;
        this.popsizeVec = popsizeVec.dup;
        this.joins = joins;
        this.joins_index = 0;
        this.mu = mu;
    }
    
    @property size_t P() {
        return nVec.length;
    }
    
    double coal_rate(size_t k) {
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
    auto nVec = [500UL, 500, 500];
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