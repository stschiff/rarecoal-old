import std.stdio;
import std.getopt;
import mainProb;
import mainLogl;
import mainMaxl;
import mainMCMC;

double mu = 1.25e-8;
int n0 = 20000;
int lingen = 400;
double tMax=20.0;

version(unittest) {
    void main() {
        writeln("all tests passed");
    }
}
else {
    void main(string[] argv) {
        
        try {
            readParams(argv);
        }
        catch(Exception e) {
            writeln(e.msg);
            printHelp();
            return;
        }
        
        if(argv.length < 2) {
            printHelp();
            return;
        }
        if(argv[1] == "prob") {
            mainProb.mainProb(argv[1 .. $], mu, n0, lingen, tMax);
        }
        else if(argv[1] == "logl") {
            mainLogl.mainLogl(argv[1 .. $], mu, n0, lingen, tMax);
        }
        else if(argv[1] == "maxl") {
            mainMaxl.mainMaxl(argv[1 .. $], mu, n0, lingen, tMax);
        }
        else if(argv[1] == "mcmc") {
            mainMCMC.mainMCMC(argv[1 .. $], mu, n0, lingen, tMax);
        }
        else {
            printHelp();
        }
    }
    
    void readParams(string[] argv) {
        getopt(argv, std.getopt.config.passThrough, "mu", &mu, "n0", &n0, "lingen", &lingen, "tMax", &tMax);
    }
    
    void printHelp() {
        writeln("rarecoal [Options] <subprogram>
Options:
        --mu         mutation rate per generation [1.25e-8]
        --n0         effective population size, used for scaling time [20000]
        --lingen     nr. of generations to propagate linearly, before crossing over to exponential intervals [400]
        --tMax       maximum time to which to run, in 2n0 generations [20.0]
Subprograms:        
        prob         compute the probability of a single configuration of derived alleles
        logl         compute the likelihood of a whole data set given parameters
        maxl         find maximum likelihood parameters for a given data set
");
    }
}