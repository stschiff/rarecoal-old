import model;

struct Params_t {
    double mu = 1.25e-8;
    int n0 = 20000;
    int lingen = 400;
    double tMax = 20.0;
    Join_t[] joins;
    Migration_t[] migrations;
    double[] popsizeVec;
    double[] leaf_times;
    int minFreq;
    int[][] exclude_list;
    int[] indices;
    int max_af;
    long nrCalledBases;
}
