import std.stdio;
import std.exception;
import std.string;
import std.algorithm;
import std.array;
import std.conv;

class Data {
    size_t[] nVec;
    size_t[immutable(size_t[])] counts;
    size_t higher;
    size_t max_m;
    size_t[][] standardOrder;

    this(string fn, size_t max_af=10) {
        char[] line;
        auto f = File(fn, "r");
        f.readln(line);
        enforce(line[0..2] == "N=", "Expect 'N=' in file");
        nVec = line[2..$].strip().split(",").map!"a.to!size_t()"().array();
        auto P = nVec.length;
        f.readln(line);
        enforce(line[0..6] == "MAX_M=", "Expect 'N=' in file");
        max_m = line[6..$].strip().to!size_t();
        if(max_m > max_af)
            max_m = max_af;
        foreach(line_; f.byLine) {
            auto fields = line_.strip().split();
            auto count = fields[1].to!size_t();
            if(fields[0] == "HIGHER")
                higher += count;
            else {
                auto key = fields[0].split(",").map!"a.to!size_t()"().array();
                if(key.reduce!"a+b"() > max_af)
                    higher += count;
                else
                    counts[key.idup] = count;
            }
        }
        standardOrder = standardConfigOrder(P, max_m);
    }
}

size_t[][] standardConfigOrder(size_t P, size_t m) {
    bool[immutable(size_t)[]] lookup;
    auto config = new size_t[P];
    auto order = [config];
    lookup[config.idup] = true;
    foreach(i; 0 .. m) {
        size_t[][] new_configs;
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
    assert(data.nVec == [500UL, 500]);
    
    data = new Data(fn, 2);
    assert(data.counts.length == 6);
    assert(data.higher == 5922);
}

