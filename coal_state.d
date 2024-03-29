import std.range;
import std.array;
import std.algorithm;
import std.conv;
import std.stdio;
import std.string;
import std.math;
import std.exception;
import model;

double approx_exp(double arg) {
    return exp(arg);
    // if(abs(arg) < 0.05)
    //     return 1.0 + arg;
    // else
    //     return exp(arg);
}

class CoalState {
    const Model model;
    ModelState modelState;
    int[] max_m;
    int m;
    double[] a, a_buf;
    double[][] b, b_buf;
    double d, t;
    
    this(in Model model, in int[] init)
    in {
        assert(model.P == init.length);
    }
    body {
        this.model = model;
        this.modelState = new ModelState(model);
        max_m = init.dup;
        m = max_m.reduce!"a+b"();
        a = zip(model.nVec, init).map!"to!double(a[0] - a[1])"().array();
        enforce(a.all!"a>=0.0"(), "nr of derived alleles exceeds haplotypes");
        a_buf = new double[model.P];
        b = new double[][](model.P, m + 1);
        b_buf = new double[][](model.P, m + 1);
        foreach(k; 0 .. model.P) {
            b[k][] = 0.0;
            b_buf[k][] = 0.0;
        }
        d = 0.0; // This is the cumulative version of c
        t = 0.0;
        foreach(k, s; init)
            b[k][s] = 1.0;
    }
    
    void step(double new_t)
    in {
        assert(new_t > t);
    }
    body {
        while(new_t > modelState.next_join().t) {
            auto jn = modelState.next_join();
            auto t_delta = jn.t - t;
            update_all(t_delta);
            perform_join(jn.k, jn.l);
            modelState.join_done();
            t = jn.t;
        }
        auto t_delta = new_t - t;
        update_all(t_delta);
        t = new_t;
    }
    
    void update_all(double t_delta) {
        update_a(t_delta);
        update_b(t_delta);
        foreach(mig; modelState.migrations) {
            perform_migration(mig.k, mig.l, mig.r, t_delta);
            perform_migration(mig.l, mig.k, mig.r, t_delta);
        }
    }
    
    void update_a(double t_delta) {
        foreach(int k, aa; a) {
            if(t + t_delta > model.leaf_times[k]) {
                auto lambda_ = modelState.coal_rate(k);
                a_buf[k] = aa * approx_exp(-(aa - 1.0) / 2.0 * lambda_ * t_delta);
            }
            else
                a_buf[k] = aa;
        }
    }
    
    void update_b(double t_delta) {
        foreach(k; 0 .. model.P) {
            if(t + t_delta > model.leaf_times[k]) {
                auto lambda_ = modelState.coal_rate(k);
                b_buf[k][] = 0.0;
                foreach(i; 0 .. max_m[k] + 1) {
                    b_buf[k][i] = b[k][i] * approx_exp(-(i * (i - 1.0) / 2.0 + i * a[k]) * lambda_ * t_delta);
                    if(i < max_m[k])
                        b_buf[k][i] += b[k][i + 1] * (1.0 - approx_exp(-i * (i + 1.0) / 2.0 * lambda_ * t_delta));
                }
            }
            else
                b_buf[k][] = b[k][];
        }
        foreach(k; 0 .. model.P)
            d += compute_c(k) * t_delta;
        swap(a, a_buf);
        swap(b, b_buf);
    }

    void perform_join(int k, int l) {
        if(empty(l))
            throw new IllegalModelException(format("merge (%s=>%s) at time %s: empty source population", l, k, t));
        if(empty(k))
            throw new IllegalModelException(format("merge (%s=>%s) at time %s: empty target population", l, k, t));
        if(t <= model.leaf_times[k])
            throw new IllegalModelException(format("merge (%s=>%s) at time %s: too young for target population", l, k, t));
        if(t <= model.leaf_times[l])
            throw new IllegalModelException(format("merge (%s=>%s) at time %s: too young for source population", l, k, t));
        a[k] += a[l];
        a[l] = 0.0;
        // stderr.writeln("before join: ", b);
        auto new_max_m = min(m, max_m[k] + max_m[l]);
        auto new_b = new double[m + 1];
        new_b[] = 0.0;
        foreach(i; 0 .. new_max_m + 1)
            foreach(j; 0 .. i + 1)
                if(j <= max_m[k] && i - j <= max_m[l])
                    new_b[i] += b[k][j] * b[l][i - j];
        b[k] = new_b;
        b[l][] = 0.0;
        b[l][0] = 1.0;
        // stderr.writefln("after join %s->%s: %s", l, k, b);
        max_m[k] = new_max_m;
        max_m[l] = 0;
    }
    
    bool empty(int k) {
        return a[k] + reduce!"a+b"(0.0, b[k][1..$]) == 0;
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
    
    void perform_migration(int k, int l, double r, double t_delta) {
        if(empty(k) || empty(l) || r == 0 || t <= model.leaf_times[k] || t <= model.leaf_times[l])
            return;
        //print "before migration, b={}".format(self.b)
        auto new_ak = a[k] + a[l] * (1.0 - approx_exp(-r * t_delta));
        auto new_al = a[l] * approx_exp(-r * t_delta);
        a[k] = new_ak;
        a[l] = new_al;
        double[] new_bkiVec;
        double[] new_bliVec;
        foreach(i; 0 .. m + 1) {
            auto tot_rate = reduce!"a+b"(0.0, iota(1, max_m[l] + 1).map!(j => b[l][j] * j * r * t_delta)());
            auto new_bki = b[k][i] * approx_exp(-tot_rate);
            if(i > 0)
                new_bki += b[k][i - 1] * (1.0 - approx_exp(-tot_rate));
            auto new_bli = b[l][i] * approx_exp(-i * r * t_delta);
            if(i < max_m[l])
                new_bli += b[l][i + 1] * (1.0 - approx_exp(-(i + 1) * r * t_delta));
            new_bkiVec ~= new_bki;
            new_bliVec ~= new_bli;
        }
        
        b[k] = new_bkiVec;
        b[l] = new_bliVec;
        // print "after migration, b={}".format(self.b)
        auto new_max_m = min(m, max_m[k] + max_m[l]);
        max_m[k] = new_max_m;
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
    cs.step(0.11);
    assert(cs.b[1][2] == 0.0, text(cs.b));
    assert(cs.b[0][3] > 0.0);
    
    
}