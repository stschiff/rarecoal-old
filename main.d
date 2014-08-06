import std.stdio;
import mainProb;
import mainLogl;
import mainMaxl;


version(unittest) {
    void main() {
        writeln("all tests passed");
    }
}
else {
    void main(string[] argv) {
        if(argv.length == 1) {
            printHelp();
            return;
        }
        if(argv[1] == "prob") {
            mainProb.mainProb(argv[1 .. $]);
        }
        else if(argv[1] == "logl") {
            mainLogl.mainLogl(argv[1 .. $]);
        }
        else if(argv[1] == "maxl") {
            mainMaxl.mainMaxl(argv[1 .. $]);
        }
        else {
            printHelp();
        }
    }
    
    void printHelp() {
        writeln("choose program [prob, logl, maxl]");
    }
}