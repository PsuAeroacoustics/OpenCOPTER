module opencopter.airfoilmodels.c81;

import opencopter.airfoilmodels;
import opencopter.config;
import opencopter.math;
import opencopter.memory;

import std.stdio;
import std.array;
import std.math;
import std.algorithm;
import std.range;
import std.conv;

size_t findindx(double[] array, double query) {
    // find index of the greatest element less then query
    size_t indx;
    foreach(i, ref a; array.enumerate) {
        if (a > query) {
            if(i > 0) {
                indx = i - 1;
            } else {
                indx = 0;
            }
            
            break;
        }
    }
    return indx;
}

final class Coefftable {
    double[] aoa;
    double[] mach;
    double[][] vals;
    double alpha_query;
    double mach_query;

    this(double[] aoa, double[] mach, double[][] vals)
    {
        this.aoa = aoa;
        this.mach = mach;
        this.vals = vals;

        checkdim();
    }

    void checkdim() {
        // checks dimentions  of variables
        if (mach.length < 2) {
            writeln("Atleast two values required for Mach Number");
        } if (aoa.length < 2) {
            writeln("Atleast two values required for Angle of Attack");
        }

        auto cols = mach.length;
        auto rows = aoa.length;

        if (cols != vals.length) {
            writeln("Inconsistant no. of mach");
            debug writeln("Number of Mach number:",cols);
            writeln("Number of Coeffiecient values:",vals.length);
        }

        if (rows != vals.length) {
            writeln("Inconsistant no. of aoa");
            writeln("Number of aoa:",rows);
            writeln("Number of Coeffiecient values:",vals[0].length);
        }
    }

    auto interpolation(double alpha_query, double mach_query) {
        double val_req;
        if (canFind(aoa, alpha_query) && canFind(mach, mach_query)) {
            debug writeln("Requaired Angle of attack and Mach number are in the file");
            immutable alpha_indx = countUntil(aoa[], alpha_query);
            immutable mach_indx = countUntil(mach[], mach_query);
            val_req = vals[alpha_indx][mach_indx];
        } else if (canFind(aoa, alpha_query)) {
            debug writeln("Requaired Angle of attack is in the file");
            
            immutable alpha_indx = countUntil(aoa[], alpha_query);
            debug writeln(aoa[alpha_indx]);
            
            immutable mach_indx = findindx(mach[], mach_query);
            debug writeln(mach[mach_indx], '\t', mach[mach_indx + 1]);

            debug writeln(vals[alpha_indx][mach_indx + 1], "\t", vals[alpha_indx][mach_indx]);

            val_req = ((vals[alpha_indx][mach_indx + 1] - vals[alpha_indx][mach_indx]) /
                    (mach[mach_indx + 1] - mach[mach_indx])) * (mach_query - mach[mach_indx]) + vals[alpha_indx][mach_indx];
        } else if (canFind(mach, mach_query)) {
            debug writeln("Requaired Mach number is in the file");
            
            immutable mach_indx = countUntil(mach[], mach_query);
            debug writeln(this.mach[mach_indx]);
            
            immutable alpha_indx = findindx(aoa, alpha_query);
            debug writeln(aoa[alpha_indx], aoa[alpha_indx + 1]);
            
            val_req = ((vals[alpha_indx + 1][mach_indx] - vals[alpha_indx][mach_indx]) /
                    (aoa[alpha_indx + 1] - aoa[alpha_indx])) * (alpha_query - aoa[alpha_indx]) + vals[alpha_indx][mach_indx];
        } else {
            debug writeln("Interpolating through the data");
            
            immutable alpha_indx = findindx(aoa[], alpha_query);
            debug writeln(aoa[alpha_indx], aoa[alpha_indx + 1]);
            
            immutable mach_indx = findindx(mach[], mach_query);
            debug writeln(mach[mach_indx], mach[mach_indx + 1]);
            
            double val_intr1 = ((vals[alpha_indx + 1][mach_indx] - vals[alpha_indx][mach_indx]) /
                    (aoa[alpha_indx + 1] - aoa[alpha_indx])) * (alpha_query - aoa[alpha_indx]) + vals[alpha_indx][mach_indx];
            double val_intr2 = ((vals[alpha_indx + 1][mach_indx + 1] - vals[alpha_indx][mach_indx + 1]) /
                    (aoa[alpha_indx + 1] - aoa[alpha_indx])) * (alpha_query - aoa[alpha_indx]) + vals[alpha_indx][mach_indx + 1];
            
            val_req = ((val_intr2 - val_intr1) / (mach[mach_indx + 1] - mach[mach_indx])) * (mach_query - mach[mach_indx]) + val_intr1;
        }
        return val_req;
    }

    auto interpolation(Chunk alpha_query, Chunk mach_query) {
        Chunk val_req;
        foreach(i; 0..chunk_size) {
            val_req[i] = interpolation(alpha_query[i], mach_query[i]);
        }
        return val_req;
    }
}

class C81: AirfoilModel {

    string airfoilname;

    Coefftable CL;
    Coefftable CD;
    Coefftable CM;

    this(string airfoilname,
        double[] aoa_L, double[] mach_l, double[][] Cl,
        double[] aoa_D, double[] mach_d, double[][] Cd,
        double[] aoa_M, double[] mach_m, double[][] Cm)
    {

        this.airfoilname = airfoilname;
        this.CL = new Coefftable(aoa_L, mach_l, Cl);
        this.CD = new Coefftable(aoa_D, mach_d, Cd);
        this.CM = new Coefftable(aoa_M, mach_m, Cm);
    }

    override double lift_curve_slope() {
        return 2.0*PI;
    }

    override double zero_lift_aoa() {
        return 0;
    }

    override Chunk get_Cl(Chunk alpha_query, Chunk mach_query) {
        return CL.interpolation(alpha_query, mach_query);
    }

    override Chunk get_Cd(Chunk alpha_query, Chunk mach_query) {
        return CD.interpolation(alpha_query, mach_query);
    }

    override double get_Cl(double alpha_query, double mach_query) {
        return CL.interpolation(alpha_query, mach_query);
    }

    override double get_Cd(double alpha_query, double mach_query) {
        return CD.interpolation(alpha_query, mach_query);
    }

    Chunk get_Cm(Chunk alpha_query, Chunk mach_query) {
        return CM.interpolation(alpha_query, mach_query);
    }

    double get_Cm(double alpha_query, double mach_query) {
        return CM.interpolation(alpha_query, mach_query);
    }
}

auto load_c81_file(string filename) {
    /*reads airfoil tabel from C81 formatted File*/
    bool multilinedata;
    auto file = File(filename, "r");
    writeln("file loaded");
    auto blade_param = to!(double[])(split(file.readln));
    double nfoil = blade_param[0]; // number of airfoils
    double ithick = blade_param[1];
    double[] xfoil = to!(double[])(split(file.readln)); // location of airfoil sections
    if (ithick == 0)
    {
        double[] thickness = to!(double[])(split(file.readln)); // max thickness of airfoil at the section
    }
    auto com1 = file.readln;
    auto com2 = file.readln;
    auto header = file.readln;
    //writeln(header[31..33]);
    string airfoilname = to!string(header[0 .. 30]);
    int nmach_l = to!int(header[30 .. 32]);
    int naoa_L = to!int(header[32 .. 34]);
    int nmach_d = to!int(header[34 .. 36]);
    int naoa_D = to!int(header[36 .. 38]);
    int nmach_m = to!int(header[39 .. 40]);
    int naoa_M = to!int(header[40 .. 42]);

    // lift
    if (nmach_l > 9)
    {
        multilinedata = true;
    }
    // read mach number
    double[] mach_l = to!(double[])(split(file.readln));
    if (multilinedata)
    {
        mach_l = (mach_l ~ to!(double[])(split(file.readln)));
    }
    //writeln(mach_l[0]);
    // read alpha and Cl
    double[] aoa_L = new double[naoa_L];
    double[][] CL = new double[][](naoa_L, nmach_l);
    for (int i = 0; i < naoa_L; i++)
    {
        double[] line = split(file.readln).map!(x => to!double(x)).array;
        if (multilinedata)
        {
            line = (line ~ split(file.readln).map!(x => to!double(x)).array);
        }
        aoa_L[i] = line[0];
        CL[i][0 .. $] = line[1 .. $];
    }
    //writeln(naoa_L);
    //writeln(CL[38][0..$]);
    multilinedata = false;
    writeln("Read the Cl data successfully");
    writeln("aoa table: \n ", aoa_L);
    // Drag
    if (nmach_d > 9)
    {
        multilinedata = true;
    }
    // read mach number
    double[] mach_d = to!(double[])(split(file.readln));
    if (multilinedata)
    {
        mach_d = (mach_d ~ to!(double[])(split(file.readln)));
    }
    // read alpha and Cd
    double[] aoa_D = new double[naoa_D];
    double[][] CD = new double[][](naoa_D, nmach_d);
    for (int i = 0; i < naoa_D; i++)
    {
        double[] line = to!(double[])(split(file.readln));
        if (multilinedata)
        {
            line = (line ~ to!(double[])(split(file.readln)));
        }
        aoa_D[i] = line[0];
        CD[i][0 .. $] = line[1 .. $];
    }
    multilinedata = false;
    writeln(naoa_D);
    writeln(CD.length);

    writeln("Read the Cd data successfully");
    // Moment
    if (nmach_m > 9)
    {
        multilinedata = true;
    }
    // read mach number
    double[] mach_m = to!(double[])(split(file.readln));
    if (multilinedata)
    {
        mach_m = (mach_m ~ to!(double[])(split(file.readln)));
    }
    // read alpha and Cm
    double[] aoa_M = new double[naoa_M];
    double[][] CM = new double[][](naoa_M, nmach_m);
    for (int i = 0; i < naoa_M; i++)
    {
        double[] line = to!(double[])(split(file.readln));
        if (multilinedata)
        {
            line = (line ~ to!(double[])(split(file.readln)));
        }
        aoa_M[i] = line[0];
        CM[i][0 .. $] = line[1 .. $];
    }
    writeln("Read the Cm data successfully");
    return new C81(airfoilname,
        aoa_L, mach_l, CL,
        aoa_D, mach_d, CD,
        aoa_M, mach_m, CM
    );
}
