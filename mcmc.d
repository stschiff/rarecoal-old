import std.math;
import std.random;
import std.range;
import std.mathspecial;
import std.exception;
import std.stdio;

class MCMC(MinfuncType) {
    MinfuncType minfunc;
    int N;
    double[][] trace;
    double[] fvalues;
    bool[] success_trace;
    int nrSuccesses;
    double[] currentPoint;
    double[] successRate;
    double[] stepWidths;
    double f;
    double fnext;
    string tracefile;
    
    this(MinfuncType minfunc, double[] init, string tracefile) {
        this.minfunc = minfunc;
        this.tracefile = tracefile;
        N = minfunc.nrParams;
        enforce(init.length == N);
        successRate = new double[N];
        successRate[] = 0.44;
        currentPoint = init.dup;
        f = minfunc(currentPoint);
        stepWidths = init.map!"a*0.01"().array();
    }
    
    void run(int nrSteps) {
        auto trials = 0;
        while(trace.length < nrSteps) {
            foreach(i; iota(N).randomCover()) {
                auto prop = proposal(i);
                currentPoint[i] += prop;
                // stderr.writeln("trying ", currentPoint);
                auto fnext = minfunc(currentPoint);
                trials += 1;
                if(accept(fnext, f)) {
                    trace ~= currentPoint.dup;
                    fvalues ~= fnext;
                    success_trace ~= true;
                    nrSuccesses += 1;
                    f = fnext;
                    successRate[i] = 0.99 * successRate[i] + 0.01;
                    if(trace.length >= nrSteps)
                        break;
                }
                else {
                    trace ~= currentPoint.dup;
                    fvalues ~= fnext;
                    success_trace ~= false;
                    currentPoint[i] -= prop;
                    successRate[i] *= 0.99;
                }
                if(trials == 10 || trials == 20 || trials == 50 || trials % 100 == 99) {
                    stderr.writefln("Trial %s, Successful steps %s/%s, adapting step widths", trials, nrSuccesses, nrSteps);
                    stderr.writefln("Current parameters: %s", currentPoint);
                    stderr.writefln("Current function value: %s", f);
                    adaptStepWidths();
                    reportTraces();
                }
            }
        }
    }
    
    bool accept(double fnext, double f) {
        if(fnext < f)
            return true;
        auto rnd = uniform(0.0, 1.0);
        if(rnd < exp(f - fnext))
            return true;
        return false; 
    }
    
    double proposal(int i) {
        auto rnd = uniform(0.0, 1.0);
        return normalDistributionInverse(rnd) * stepWidths[i];
    }
    
    void adaptStepWidths() {
        stderr.writeln("Sucess rates:, ", successRate);
        foreach(i; 0 .. N) {
            if(successRate[i] < 0.29)
                stepWidths[i] /= 1.5;
            if(successRate[i] > 0.59)
                stepWidths[i] *= 1.5;
        }
        stderr.writeln("new step widths: ", stepWidths);
    }
    
    void reportTraces() {
        auto f = File(tracefile, "w");
        f.write(minfunc.paramNames().join("\t"));
        f.write("\tsuccess");
        f.writeln("\tscore");
        foreach(i, point; trace) {
            f.write(point.map!"text(a)"().join("\t"));
            f.writef("\t%s", success_trace[i]);
            f.writefln("\t%20.2f", fvalues[i]);
        }
        f.close();
    }
}
