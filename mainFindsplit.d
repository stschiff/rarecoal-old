import std.stdio;
import std.getopt;
import params;
import data;

int newBranch;
Params_t p;
int max_af = 10;
Data input_data;
string spectrumfile = "/dev/null";

void mainFindsplit(string[] argv, Params_t params_) {
    p = params_;
}

void readParams(string[] argv) {
    getopt(argv, std.getopt.config.caseSensitive,
           "max_af|m"  , &max_af,
           "spectrumfile|s", &spectrumfile);
    
    enforce(argv.length == 2, "need more arguments");
    input_data = new Data(argv[1], max_af);
    enforce(p.popsizeVec.length == input_data.nVec.length);
}