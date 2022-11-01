module opencopter.read_afdata;

import opencopter.math;

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
    }

    void checkdim(double[][] vals)
    {
        /* checks dimentions  of variables */
        if (this.mach.length < 2)
        {
            writeln("Atleast two values required for Mach Number");
        }
        if (this.aoa.length < 2)
        {
            writeln("Atleast two values required for Angle of Attack");
        }

        auto cols = this.mach.length;
        auto rows = this.aoa.length;

        if (cols != this.vals.length)
        {
            writeln("Inconsistant no. of mach");
            writeln(cols);
            writeln(this.vals.length);
        }

        if (rows != this.vals.length)
        {
            writeln("Inconsistant no. of aoa");
            writeln(rows);
            writeln(this.vals[0].length);
        }
    }

    auto interpolation(double alphaQuery, double machQuery)
    {
        if (canFind(this.aoa[], alphaQuery) && canFind(this.mach[], machQuery))
        {
            writeln("Requaired Angle of attack and Mach number are in the file");
            auto alpha_indx = countUntil(this.aoa[], alphaQuery);
            auto mach_indx = countUntil(this.mach[], machQuery);
            double val_req = this.vals[alpha_indx][mach_indx];
            return val_req;
        }
        else if (canFind(this.aoa, alphaQuery))
        {
            writeln("Requaired Angle of attack is in the file");
            auto alpha_indx = countUntil(this.aoa[], alphaQuery);
            writeln(this.aoa[alpha_indx]);
            int mach_indx = findindx(this.mach, machQuery);
            writeln(this.mach[mach_indx], '\t', this.mach[mach_indx + 1]);
            writeln(this.vals[alpha_indx][mach_indx + 1], "\t", this.vals[alpha_indx][mach_indx]);
            double val_req = ((this.vals[alpha_indx][mach_indx + 1] - this.vals[alpha_indx][mach_indx]) /
                    (
                        this.mach[mach_indx + 1] - this.mach[mach_indx])) * (
                machQuery - this.mach[mach_indx]) + this.vals[alpha_indx][mach_indx];
            return val_req;
        }
        else if (canFind(this.mach, machQuery))
        {
            writeln("Requaired Mach number is in the file");
            auto mach_indx = countUntil(this.mach[], machQuery);
            writeln(this.mach[mach_indx]);
            int alpha_indx = findindx(this.aoa, alphaQuery);
            writeln(this.aoa[alpha_indx], this.aoa[alpha_indx + 1]);
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
            writeln(this.aoa[alpha_indx], this.aoa[alpha_indx + 1]);
            int mach_indx = findindx(this.mach, machQuery);
            writeln(this.mach[mach_indx], this.mach[mach_indx + 1]);
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
    }

}

class C81
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

        this.CL.checkdim(Cl);
        this.CD.checkdim(Cd);
        this.CM.checkdim(Cm);
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

    double getCl(double alphaQuery, double machQuery)
    {
        return this.CL.interpolation(alphaQuery, machQuery);
    }

    double getCd(double alphaQuery, double machQuery)
    {
        return this.CD.interpolation(alphaQuery, machQuery);
    }

    double getCm(double alphaQuery, double machQuery)
    {
        return this.CM.interpolation(alphaQuery, machQuery);
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

/* Airfoil Cl and Cd datt based on AERODAS model*/

class Input_param
{
    double[] alpha; // alpha,cl, and cd from 2D airoil data file
    double[] CL;
    double[] CD;
    double A0; //aoa with CL1=0 for all aspect ratios(AR) (CL1: lift coeff. in pre-stall regime)
    double ACL1_2D; //aoa at maximum pre-stall lift (for infinite AR balde)
    double cl1max_2D; //maximum pre-stall lift coefficient at aoa= ACL1_2D (for infinite AR balde)
    double S1_2D; //slope of linear segment of lift curve slope (for infinite AR balde)
    double CD0; //minimum drag coefficient for aoa=A0 for all ARs.
    double ACD1_2D; //aoa at maximum pre-stall drag (for infinite AR balde)
    double cd1max_2D; //maximum pre-stall drag coefficient at aoa= ACD1_2D (for infinite AR balde)
    int A0_indx;
    this(double[] alpha, double[] CL, double[] CD)
    {
        this.alpha = alpha;
        this.CL = CL;
        this.CD = CD;
        this.A0 = A0;
        this.ACL1_2D = ACL1_2D;
        this.cl1max_2D = cl1max_2D;
        this.S1_2D = S1_2D;
        this.ACD1_2D = ACD1_2D;
        this.CD0 = CD0;
        this.cd1max_2D = cd1max_2D;

        if (this.alpha.length != this.CL.length || this.alpha.length != this.CD.length || this.CL.length != this
            .CD.length)
        {
            writeln("Error: Angle of attack, Cl, and Cd must have the same length");
        }

        if (all(CL) > 0 || all(CL) < 0)
        {
            writeln("Error: CL should contain at leaset one positive and one negative value");
        }

        if (canFind(this.CL, 0))
        {
            this.A0 = this.alpha[countUntil(this.CL, 0)];
        }
        else
        {
            A0_indx = findindx(this.CL, 0);
            writeln("aoa_0",this.alpha[A0_indx]);
            A0 = this.alpha[A0_indx] - ((this.alpha[A0_indx+1]-this.alpha[A0_indx])/(this.CL[A0_indx+1]-this.CL[A0_indx]))*this.CL[A0_indx];
        }

        this.cl1max_2D = maxElement(this.CL);
        this.ACL1_2D = this.alpha[countUntil(this.CL[], this.cl1max_2D)];
        this.S1_2D = (this.CL[A0_indx + 3] - this.CL[A0_indx - 3]) / (
            this.alpha[A0_indx + 3] - this.alpha[A0_indx - 3]);
        this.CD0 = minElement(this.CD);
        this.cd1max_2D = maxElement(this.CD);
        this.ACD1_2D = this.alpha[countUntil(this.CD[], this.cd1max_2D)];
    }

}

class ClCdmodel
{
    double[] alpha; // alpha,cl, and cd from 2D airoil data file
    double[] CL;
    double[] CD;
    double AR; //Aspect ratio of wing
    double tbyc; //thickness to chord ratio
    Input_param ip;
    double alphaquery;
    double CL1; //pre-stall lift coefficient
    double CD1; //pre-stall drag coefficient
    double CL2; //post-stall lift coefficient
    double CD2; //post-stall drag coefficient
    this(double[] alpha, double[] CL, double[] CD, double AR, double tbyc)
    {
        this.alpha = alpha;
        this.CL = CL;
        this.CD = CD;
        this.AR = AR;
        this.tbyc = tbyc;

        this.ip = new Input_param(this.alpha, this.CL, this.CD);
        writeln("A0 is",this.ip.A0);
    }
        double getCL(double alphaquery)
        {
            //finite aspect ratio adjustment
            double ACL1_3D = this.ip.ACL1_2D + this.ip.cd1max_2D * 18.2 * pow(this.AR,-0.9);
            writeln("ACL1 = ",ACL1_3D);
            double S1_3D = this.ip.S1_2D / (1 + this.ip.S1_2D * 18.2 * pow(this.AR,-0.9));
            double ACD1_3D = this.ip.ACD1_2D + this.ip.cl1max_2D * 18.2 * pow(this.AR,-0.9);
            double cd1max_3D = this.ip.cd1max_2D + 0.280 * pow(this.ip.cl1max_2D,2.0) * pow(this.AR,-0.9);
            double cl1max_3D = this.ip.cl1max_2D * (0.67 + 0.33 * exp(-16 / pow(this.AR,2.0)));
            // pre-stall Cl calculation
            double RCL1 = S1_3D * (ACL1_3D - this.ip.A0) - cl1max_3D; //reduction from extension of linear segment of lift curve to CL1max
            double N1 = 1 + cl1max_3D / RCL1;
            if (alphaquery >= this.ip.A0)
            {
                CL1 = S1_3D * (alphaquery - this.ip.A0) - RCL1 * pow(
                    (alphaquery - this.ip.A0) / (ACL1_3D - this.ip.A0),N1);
            }
            else
            {
                CL1 = S1_3D * (alphaquery - this.ip.A0) + RCL1 * pow(
                    (this.ip.A0 - alphaquery) / (ACL1_3D - this.ip.A0),N1);
            }
            //writeln("CL1 = ", CL1);
            //post stall regim

            double getCL2(double aoa)
            {
                double ans;
                double F1 = 1.190 * (1 - pow(this.tbyc,2));
                double F2 = 0.65 + 0.35 * exp(-pow((9.0 / this.AR),2.3));
                double cl2max = F1 * F2;
                double RCL2 = 1.632 - cl2max;
                double N2 = 1 + cl2max / RCL2;
                if (0 < aoa && aoa < ACL1_3D)
                {
                    ans = 0;
                }
                else if (ACL1_3D <= aoa && aoa <= 92.0)
                {
                    ans = -0.032 * (aoa - 92.0) - RCL2 * pow(((92.0 - aoa) / 51.0),N2);
                }
                else if (aoa > 92)
                {
                    ans = -0.032 * (aoa - 92.0) + RCL2 * pow(((aoa - 92.0) / 51.0),N2);
                }
                return ans;
            }

            if (alphaquery >= 0)
            {
                CL2 = getCL2(alphaquery);
                //writeln("CL2 = ",CL2);
            }
            else
            {
                double updated_alphaquery = -alphaquery + 2.0 * this.ip.A0;
                //writeln("alpha_neg = ",updated_alphaquery);
                CL2 = -1.0 * getCL2(updated_alphaquery);
                //writeln("CL2 = ",CL2);
            }

            if (alphaquery >= this.ip.A0)
            {
                return max(CL1, CL2);
            }
            else
            {
                return min(CL1, CL2);
            }
        }

        double getCD(double alphaquery)
        {   
            //finite aspect ratio adjustment
            double ACL1_3D = this.ip.ACL1_2D + this.ip.cd1max_2D * 18.2 * pow(this.AR,-0.9);
            double S1_3D = this.ip.S1_2D / (1 + this.ip.S1_2D * 18.2 * pow(this.AR,-0.9));
            double ACD1_3D = this.ip.ACD1_2D + this.ip.cl1max_2D * 18.2 * pow(this.AR,-0.9);
            double cd1max_3D = this.ip.cd1max_2D + 0.280 * pow(this.ip.cl1max_2D,2.0) * pow(this.AR,-0.9);
            writeln("CDmax_3D = ",cd1max_3D);
            double cl1max_3D = this.ip.cl1max_2D * (0.67 + 0.33 * exp(-16 / pow(this.AR,2.0)));
                        
            //pre-stall Cd calculation
            if ((2 * this.ip.A0 - ACD1_3D) <= alphaquery && alphaquery <= ACD1_3D)
            {
                CD1 = this.ip.CD0 + (cd1max_3D - this.ip.CD0) * pow(
                    (alphaquery - this.ip.A0) / (ACD1_3D - this.ip.A0),2);
            }
            else
            {
                CD1 = 0;
            }
            writeln("CD1 = ",CD1);

            double getCD2(double aoa)
            {
                double ans1;
                double G1 = 2.300 * exp(-1 * pow(0.65 * tbyc,0.9));
                double G2 = 0.52 + 0.48 * exp(-1 * pow((0.65 / this.AR),1.1));
                double cd2max = G1 * G2;
                writeln("cd2max = ",cd2max);
                if ((2 * this.ip.A0 - ACL1_3D) < aoa && aoa < ACL1_3D)
                { //(check this condition(might be typo in the paper pg18.) 
                    ans1 = 0;
                }
                else if (aoa > ACD1_3D)
                {
                    ans1 = cd1max_3D + (cd2max - cd1max_3D) * sin((aoa - ACD1_3D) * PI*0.5 / (
                            90.0 - ACD1_3D));
                }
                return ans1;
            }

            if (alphaquery > (2 * this.ip.A0 - ACD1_3D))
            {
                CD2 = getCD2(alphaquery);
            }
            else
            {
                double updated_alphaquery = -alphaquery + 2 * this.ip.A0;
                CD2 = getCD2(updated_alphaquery);
            }
            writeln("CD2 = ",CD2);

            return max(CD1, CD2);
        }
}

auto load2Dfile(string filename, double AR, double tbyc)
{
    /*reads 2D airfoil Data file*/
    auto file = File(filename, "r");
    writeln("file loaded");
    // not so important things
    auto line = file.readln;
    line = file.readln;
    line = file.readln;
    auto afname = file.readln;
    line = file.readln;
    line = file.readln;
    line = file.readln;
    line = file.readln;
    auto Re = file.readln;
    writeln("Re line is",Re);
    line = file.readln;
    //writeln(line);
    line = file.readln;
    //writeln(line);
    line = file.readln;
    //writeln(line);
    // airfoil data begines
    double[] aoa;
    double[] cl;
    double[] cd;

    while (1)
    {
        double[] line1 = split(file.readln).map!(x => to!double(x)).array;
        //writeln(line1);
        if (line1 is null)
        {
            break;
        }
        else
        {
            aoa = (aoa ~ line1[0]);
            //writeln("alpha = ",aoa);
            cl = (cl ~ line1[1]);
            cd = (cd ~ line1[2]);
        }

    }
    writeln("File reading completed");
    auto afmodel = new ClCdmodel(aoa, cl, cd, AR, tbyc);
    return afmodel;
}
