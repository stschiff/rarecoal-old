import std.stdio;
import std.exception;
import std.string;
import std.algorithm;
import std.array;
import std.conv;

class Data {
    int[] nVec;
    long[immutable(int[])] counts;
    long higher;
    int max_m;
    int[][] standardOrder;

    this(string fn, int max_af=10) {
        char[] line;
        auto f = File(fn, "r");
        f.readln(line);
        enforce(line[0..2] == "N=", "Expect 'N=' in file");
        nVec = line[2..$].strip().split(",").map!"a.to!int()"().array();
        auto P = nVec.length.to!int();
        f.readln(line);
        enforce(line[0..6] == "MAX_M=", "Expect 'N=' in file");
        max_m = line[6..$].strip().to!int();
        if(max_m > max_af)
            max_m = max_af;
        foreach(line_; f.byLine) {
            auto fields = line_.strip().split();
            auto count = fields[1].to!long();
            if(fields[0] == "HIGHER")
                higher += count;
            else {
                auto key = fields[0].split(",").map!"a.to!int()"().array();
                if(key.reduce!"a+b"() > max_af)
                    higher += count;
                else
                    counts[key.idup] = count;
            }
        }
        standardOrder = standardConfigOrder(P, max_m);
    }
}

int[][] standardConfigOrder(int P, int m) {
    bool[immutable(int)[]] lookup;
    auto config = new int[P];
    auto order = [config];
    lookup[config.idup] = true;
    foreach(i; 0 .. m) {
        int[][] new_configs;
        foreach(o; order) {
            foreach(k; 0 .. P) {
                auto new_config = o.dup;
                new_config[k] += 1;
                auto key = new_config.idup;
                if(key !in lookup) {
                    new_configs ~= new_config;
                    lookup[key] = true;
                }
            }
        }
        order ~= new_configs;
    }
    return order;
}

unittest {
    auto order = standardConfigOrder(3, 1);
    assert(order == [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    order = standardConfigOrder(3, 2);
    assert(order == [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                     [2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1],
                     [0, 0, 2]], text(order));
    order = standardConfigOrder(10, 6);
    assert(order.length == 8008);
}


unittest {
    auto fn = "testDat/testDat.txt";
    auto data = new Data(fn);
    assert(data.counts.length == 15);
    assert(data.higher == 5275);
    assert(data.nVec == [500, 500]);
    
    data = new Data(fn, 2);
    assert(data.counts.length == 6);
    assert(data.higher == 5922);
}

