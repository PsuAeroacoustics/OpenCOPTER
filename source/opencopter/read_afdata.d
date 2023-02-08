module opencopter.read_afdata;

import opencopter.math;
import opencopter.liftmodel_interface;
import opencopter.memory;
import std.stdio;
import std.array;
import std.math;
import std.algorithm;
import std.range;
import std.conv;

int findindx(double[] array, double query)
{
    /*find index of the greatest element less then query*/
    int indx;
    for (int i = 0; i < array.length; i++)
    {
        if (array[i] > query)
        {
            indx = i - 1;
            break;
        }
    }
    return indx;
}

auto interpolation(Chunk alphaQuery, Chunk machQuery, double[] aoa, double[] mach, double[][] vals)
    {   
        double[] val_req;
        for(int i=0; i<8; i++){
            if (canFind(aoa, alphaQuery[i]) && canFind(mach, machQuery[i]))
            {
            debug writeln("Requaired Angle of attack and Mach number are in the file");
            auto alpha_indx = countUntil(aoa[], alphaQuery[i]);
            auto mach_indx = countUntil(mach[], machQuery[i]);
            val_req[i] = vals[alpha_indx][mach_indx];
            }
        else if (canFind(aoa, alphaQuery[i]))
            {
            debug writeln("Requaired Angle of attack is in the file");
            auto alpha_indx = countUntil(aoa[], alphaQuery[i]);
            debug writeln(aoa[alpha_indx]);
            int mach_indx = findindx(mach[], machQuery[i]);
            debug writeln(mach[mach_indx], '\t', mach[mach_indx + 1]);
            debug writeln(vals[alpha_indx][mach_indx + 1], "\t", vals[alpha_indx][mach_indx]);
            val_req[i] = ((vals[alpha_indx][mach_indx + 1] - vals[alpha_indx][mach_indx]) /
                    (mach[mach_indx + 1] - mach[mach_indx])) * (machQuery[i] - mach[mach_indx]) + vals[alpha_indx][mach_indx];
            }
        else if (canFind(mach, machQuery[i]))
            {
            debug writeln("Requaired Mach number is in the file");
            auto mach_indx = countUntil(mach[], machQuery[i]);
            debug writeln(this.mach[mach_indx]);
            int alpha_indx = findindx(aoa, alphaQuery[i]);
            debug writeln(aoa[alpha_indx], aoa[alpha_indx + 1]);
            val_req[i] = ((vals[alpha_indx + 1][mach_indx] - vals[alpha_indx][mach_indx]) /
                    (aoa[alpha_indx + 1] - aoa[alpha_indx])) * (alphaQuery[i] - aoa[alpha_indx]) + vals[alpha_indx][mach_indx];
            }
        else
            {
            writeln("Interpolating through the data");
            int alpha_indx = findindx(aoa[], alphaQuery[i]);
            debug writeln(aoa[alpha_indx], aoa[alpha_indx + 1]);
            int mach_indx = findindx(mach[], machQuery[i]);
            debug writeln(mach[mach_indx], mach[mach_indx + 1]);
            double val_intr1 = ((vals[alpha_indx + 1][mach_indx] - vals[alpha_indx][mach_indx]) /
                    (aoa[alpha_indx + 1] - aoa[alpha_indx])) * (alphaQuery[i] - aoa[alpha_indx]) + vals[alpha_indx][mach_indx];
            double val_intr2 = ((vals[alpha_indx + 1][mach_indx + 1] - vals[alpha_indx][mach_indx + 1]) /
                    (aoa[alpha_indx + 1] - aoa[alpha_indx])) * (alphaQuery[i] - aoa[alpha_indx]) + vals[alpha_indx][mach_indx + 1];
            val_req[i] = ((val_intr2 - val_intr1) / (mach[mach_indx + 1] - mach[mach_indx])) * (machQuery[i] - mach[mach_indx]) + val_intr1;
            }
        }
        return val_req;
    }

    void checkdim(double[] aoa,double[] mach,double[][] vals)
    {
        /* checks dimentions  of variables */
        if (mach.length < 2)
        {
            writeln("Atleast two values required for Mach Number");
        }
        if (aoa.length < 2)
        {
            writeln("Atleast two values required for Angle of Attack");
        }

        auto cols = mach.length;
        auto rows = aoa.length;

        if (cols != vals.length)
        {
            writeln("Inconsistant no. of mach");
            debug writeln("Number of Mach number:",cols);
            writeln("Number of Coeffiecient values:",vals.length);
        }

        if (rows != vals.length)
        {
            writeln("Inconsistant no. of aoa");
            writeln("Number of aoa:",rows);
            writeln("Number of Coeffiecient values:",vals[0].length);
        }
    }


class Coefftable
{
    double[] aoa;
    double[] mach;
    double[][] vals;
    double alphaQuery;
    double machQuery;

    this(double[] aoa, double[] mach, double[][] vals)
    {
        this.aoa = aoa;
        this.mach = mach;
        this.vals = vals;

        checkdim(this.aoa,this.mach,this.vals);
    }


    /*
    auto interpolation(double alphaQuery, double machQuery)
    {
        if (canFind(this.aoa[], alphaQuery) && canFind(this.mach[], machQuery))
        {
            debug writeln("Requaired Angle of attack and Mach number are in the file");
            auto alpha_indx = countUntil(this.aoa[], alphaQuery);
            auto mach_indx = countUntil(this.mach[], machQuery);
            double val_req = this.vals[alpha_indx][mach_indx];
            return val_req;
        }
        else if (canFind(this.aoa, alphaQuery))
        {
            debug writeln("Requaired Angle of attack is in the file");
            auto alpha_indx = countUntil(this.aoa[], alphaQuery);
            debug writeln(this.aoa[alpha_indx]);
            int mach_indx = findindx(this.mach, machQuery);
            debug writeln(this.mach[mach_indx], '\t', this.mach[mach_indx + 1]);
            debug writeln(this.vals[alpha_indx][mach_indx + 1], "\t", this.vals[alpha_indx][mach_indx]);
            double val_req = ((this.vals[alpha_indx][mach_indx + 1] - this.vals[alpha_indx][mach_indx]) /
                    (
                        this.mach[mach_indx + 1] - this.mach[mach_indx])) * (
                machQuery - this.mach[mach_indx]) + this.vals[alpha_indx][mach_indx];
            return val_req;
        }
        else if (canFind(this.mach, machQuery))
        {
            debug writeln("Requaired Mach number is in the file");
            auto mach_indx = countUntil(this.mach[], machQuery);
            debug writeln(this.mach[mach_indx]);
            int alpha_indx = findindx(this.aoa, alphaQuery);
            debug writeln(this.aoa[alpha_indx], this.aoa[alpha_indx + 1]);
            double val_req = ((this.vals[alpha_indx + 1][mach_indx] - this.vals[alpha_indx][mach_indx]) /
                    (
                        this.aoa[alpha_indx + 1] - this.aoa[alpha_indx])) * (
                alphaQuery - this.aoa[alpha_indx]) + this.vals[alpha_indx][mach_indx];
            return val_req;
        }
        else
        {
            writeln("Interpolating through the data");
            int alpha_indx = findindx(this.aoa, alphaQuery);
            debug writeln(this.aoa[alpha_indx], this.aoa[alpha_indx + 1]);
            int mach_indx = findindx(this.mach, machQuery);
            debug writeln(this.mach[mach_indx], this.mach[mach_indx + 1]);
            double val_intr1 = ((this.vals[alpha_indx + 1][mach_indx] - this.vals[alpha_indx][mach_indx]) /
                    (
                        this.aoa[alpha_indx + 1] - this.aoa[alpha_indx])) * (
                alphaQuery - this.aoa[alpha_indx]) + this.vals[alpha_indx][mach_indx];
            double val_intr2 = (
                (this.vals[alpha_indx + 1][mach_indx + 1] - this.vals[alpha_indx][mach_indx + 1]) /
                    (this.aoa[alpha_indx + 1] - this.aoa[alpha_indx])) * (
                alphaQuery - this.aoa[alpha_indx]) + this.vals[alpha_indx][mach_indx + 1];
            double val_req = ((val_intr2 - val_intr1) / (
                    this.mach[mach_indx + 1] - this.mach[mach_indx])) * (
                machQuery - this.mach[mach_indx]) + val_intr1;
            return val_req;
        }
    }*/

}

class C81:Coefficent
{

    string airfoilname;
    double[] aoa_L;
    double[] mach_l;
    Coefftable CL;
    double[] aoa_D;
    double[] mach_d;
    Coefftable CD;
    double[] aoa_M;
    double[] mach_m;
    Coefftable CM;
    double alphaQuery;
    double machQuery;

    this(string airfoilname,
        double[] aoa_L, double[] mach_l, double[][] Cl,
        double[] aoa_D, double[] mach_d, double[][] Cd,
        double[] aoa_M, double[] mach_m, double[][] Cm)
    {

        this.airfoilname = airfoilname;
        this.CL = new Coefftable(aoa_L, mach_l, Cl);
        this.CD = new Coefftable(aoa_D, mach_d, Cd);
        this.CM = new Coefftable(aoa_M, mach_m, Cm);

        /*this.CL.checkdim(Cl);
        this.CD.checkdim(Cd);
        this.CM.checkdim(Cm);*/
    }

    bool increascheck(double[] arr)
    {
        auto a = zip(arr.dropBackOne(), arr.dropOne()).map!"a[1]-a[0]";
        auto b = any!"a<0"(a);
        return b;
    }

    void refreshinterpolate()
    {
        /* don't know if it is nessesory or not. But its good if we have strictly inceasing Mach Number and Angle of attack.*/
        if (increascheck(this.CL.aoa))
        {
            writeln("Angle of attack for Cl is not strictly Increasing");
        }
        if (increascheck(this.CL.mach))
        {
            writeln("Mach Number for Cl is not strictly Increasing");
        }
        if (increascheck(this.CD.aoa))
        {
            writeln("Angle of attack for Cd is not strictly Increasing");
        }
        if (increascheck(this.CD.mach))
        {
            writeln("Mach Number for Cd is not strictly Increasing");
        }
        if (increascheck(this.CM.aoa))
        {
            writeln("Angle of attack for Cm is not strictly Increasing");
        }
        if (increascheck(this.CM.mach))
        {
            writeln("Mach Number for Cm is not strictly Increasing");
        }
    }

    double[] getCl(Chunk alphaQuery, Chunk machQuery)
    {
        return interpolation(alphaQuery, machQuery,this.CL.aoa,this.CL.mach,this.CL.vals);
    }

    double[] getCd(Chunk alphaQuery, Chunk machQuery)
    {
        return interpolation(alphaQuery, machQuery,this.CD.aoa,this.CD.mach,this.CD.vals);
    }

    double[] getCm(Chunk alphaQuery, Chunk machQuery)
    {
        return interpolation(alphaQuery, machQuery,this.CM.aoa,this.CM.mach,this.CM.vals);       
    }

}

auto loadfile(string filename)
{
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
        aoa_M, mach_m, CM);
}
