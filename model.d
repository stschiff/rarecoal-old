import std.typecons;

alias Tuple!(double, "t", size_t, "k", size_t, "l") Join_t;

class Model {
    const size_t[] nVec;
    const Join_t[] joins;
    double mu;
    size_t joins_index;
    
    this(in size_t[] nVec, in Join_t[] joins=[], double mu=0.0005) {
        this.nVec = nVec;
        this.joins = joins;
        this.joins_index = 0;
        this.mu = mu;
    }
    
    @property size_t P() {
        return nVec.length;
    }
    
    double coal_rate(size_t k, double t) {
        return 1.0;
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
        joins_index += 1;
    }
    
    void reset() {
        joins_index = 0;
    }
}

unittest {
    auto nVec = [500UL, 500, 500];
    auto joins = [Join_t(0.1, 0, 1), Join_t(0.2, 0, 2)];
    auto m = new Model(nVec, joins);
    assert(m.P == 3);
    assert(m.next_join() == Join_t(0.1, 0, 1));
    m.join_done();
    assert(m.next_join() == Join_t(0.2, 0, 2));
    m.join_done();
    assert(m.next_join().t == double.infinity);
}