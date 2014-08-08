import std.range;
import std.array;
import std.algorithm;
import std.conv;
import std.stdio;
import std.string;
import std.math;
import model;

double approx_exp(double arg) {
    if(abs(arg) < 0.05)
        return 1.0 + arg;
    else
        return exp(arg);
}

class CoalState {
    Model model;
    int[] max_m;
    int m;
    double[] a, a_buf;
    double[][] b, b_buf;
    double d, e, t;
    
    this(Model model, in int[] init)
    in {
        assert(model.P == init.length);
    }
    body {
        this.model = model;
        max_m = init.dup;
        m = max_m.reduce!"a+b"();
        a = zip(model.nVec, init).map!"to!double(a[0] - a[1])"().array();
        a_buf = new double[model.P];
        b = new double[][](model.P, m + 1);
        b_buf = new double[][](model.P, m + 1);
        foreach(ref bb; b)
            bb[] = 0.0;
        d = 0.0; // This is the cumulative version of c
        e = 0.0; // This is the cumulative version of a
        t = 0.0;
        foreach(k, s; init)
            b[k][s] = 1.0;
    }
    
    void step(double new_t)
    in {
        assert(new_t > t);
    }
    body {
        auto t_delta = new_t - t;
        update_a(t_delta);
        update_b(t_delta);
        update_e(t_delta);
        update_d(t_delta);
        swap(a, a_buf);
        swap(b, b_buf);
        auto jn = model.next_join();
        if(new_t > jn.t) {
            perform_join(jn.k, jn.l);
            model.join_done();
        }
        t = new_t;
    }
    
    void update_a(double t_delta) {
        foreach(int k, aa; a) {
            auto lambda_ = model.coal_rate(k);
            a_buf[k] = aa * approx_exp(-(aa - 1.0) / 2.0 * lambda_ * t_delta);
        }
    }
    
    void update_b(double t_delta) {
        foreach(k; 0 .. model.P) {
            auto lambda_ = model.coal_rate(k);
            foreach(i; 0 .. max_m[k] + 1) {
                b_buf[k][i] = b[k][i] * approx_exp(-(i * (i - 1.0) / 2.0 + i * a[k]) * lambda_ * t_delta);
                if(i < max_m[k])
                    b_buf[k][i] += b[k][i + 1] * (1.0 - approx_exp(-i * (i + 1.0) / 2.0 * lambda_ * t_delta));
            }
        }
    }

    void update_d(double t_delta) {
        foreach(k; 0 .. model.P)
            d += compute_c(k) * t_delta;
    }
    
    void update_e(double t_delta) {
        e += a.reduce!"a+b"() * t_delta;
    }
    
    void perform_join(int k, int l) {
        if(a[l] + max_m[l] == 0)
            throw new IllegalModelException(format("tried to merge %s into %s at time %s", l, k, t));
        a[k] += a[l];
        a[l] = 0.0;
        auto new_max_m = max_m[k] + max_m[l];
        auto new_b = new double[m + 1];
        new_b[] = 0.0;
        foreach(i; 0 .. new_max_m + 1)
            foreach(j; 0 .. i + 1)
                if(j <= max_m[k] && i - j <= max_m[l])
                    new_b[i] += b[k][j] * b[l][i - j];
        b[k] = new_b;
        b[l][] = 0.0;
        b[l][0] = 1.0;
        max_m[k] = new_max_m;
        max_m[l] = 0;
    }
    
    double compute_c(int k) {
        if(max_m[k] == 0)
            return 0.0;
        auto ret = b[k][1];
        foreach(l; 0 .. model.P)
            if(l != k)
                ret *= b[l][0];
        return ret;
    }
}

unittest {
    auto nVec = [500, 500, 500];
    auto joins = [Join_t(0.1, 0, 1, 1.0), Join_t(0.2, 0, 2, 1.0)];
    auto m = new Model(nVec, [1.0, 1.0, 1.0], joins);
    
    auto cs = new CoalState(m, [1, 2, 0]);
    assert(cs.a == [499, 498, 500]);
    assert(cs.b[0][1] == 1.0);
    cs.perform_join(0, 1);
    assert(cs.b[1][2] == 0.0);
    assert(cs.b[0][3] == 1.0);

    cs = new CoalState(m, [1, 2, 0]);
    assert(cs.a == [499, 498, 500]);
    assert(cs.b[0][1] == 1.0);
    cs.step(0.15);
    assert(cs.b[1][2] == 0.0);
    assert(cs.b[0][3] > 0.0);
}