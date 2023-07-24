module opencopter.airfoilmodels.aerodas;

import opencopter.airfoilmodels;
import opencopter.config;
import opencopter.math;
import opencopter.memory;

import std.stdio;
import std.array;
import std.exception;
import std.math;
import std.algorithm;
import std.range;
import std.conv;


/* Airfoil Cl and Cd based on AERODAS model*/

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

class InputParam {
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

    this(double[] _alpha, double[] _CL, double[] _CD) {
        this.alpha = _alpha;
        this.CL = _CL;
        this.CD = _CD;
        //this.A0 = A0;
        //this.ACL1_2D = ACL1_2D;
        //this.cl1max_2D = cl1max_2D;
        //this.S1_2D = S1_2D;
        //this.ACD1_2D = ACD1_2D;
        //this.CD0 = CD0;
        //this.cd1max_2D = cd1max_2D;

        enforce(alpha.length == CL.length && alpha.length == CD.length && CL.length == CD.length,
            "Angle of attack, Cl, and Cd must have the same length"
        );
        //if (alpha.length != CL.length || alpha.length != CD.length || CL.length != CD.length) {
        //    writeln("Error: Angle of attack, Cl, and Cd must have the same length");
        //}

        //if (all(CL) > 0 || all(CL) < 0) {
        //    writeln("Error: CL should contain at leaset one positive and one negative value");
        //}

        if (canFind(this.CL, 0)) {
            A0_indx = findindx(this.CL, 0);
            this.A0 = this.alpha[countUntil(this.CL, 0)];
        } else {
            A0_indx = findindx(this.CL, 0);
            debug writeln("aoa_0",this.alpha[A0_indx]);
            A0 = this.alpha[A0_indx] - ((this.alpha[A0_indx+1]-this.alpha[A0_indx])/(this.CL[A0_indx+1]-this.CL[A0_indx]))*this.CL[A0_indx];
        }

        this.cl1max_2D = maxElement(this.CL);
        this.ACL1_2D = this.alpha[countUntil(this.CL[], this.cl1max_2D)];
        debug writeln("this.CL: ", this.CL);
        debug writeln("A0_indx: ", A0_indx);
        this.S1_2D = (this.CL[A0_indx + 3] - this.CL[A0_indx - 3]) / (
            this.alpha[A0_indx + 3] - this.alpha[A0_indx - 3]);
        this.CD0 = minElement(this.CD);
        this.cd1max_2D = maxElement(this.CD);
        this.ACD1_2D = this.alpha[countUntil(this.CD[], this.cd1max_2D)];
    }
}

// Function for CL calulation in the post stall regime
double getCL2(double aoa, double tbyc, double AR, double ACL1_3D) {
    double ans;
    double F1 = 1.190 * (1 - pow(tbyc,2));
    double F2 = 0.65 + 0.35 * exp(-pow((9.0 / AR),2.3));
    double cl2max = F1 * F2;
    double RCL2 = 1.632 - cl2max;
    double N2 = 1 + cl2max / RCL2;

    if (0 < aoa && aoa < ACL1_3D) {
        ans = 0;
    } else if (ACL1_3D <= aoa && aoa <= 92.0) {
        ans = -0.032 * (aoa - 92.0) - RCL2 * pow(((92.0 - aoa) / 51.0),N2);
    } else if (aoa > 92) {
        ans = -0.032 * (aoa - 92.0) + RCL2 * pow(((aoa - 92.0) / 51.0),N2);
    }
    return ans;
}
// Function for CD calulation in the post stall regime
double getCD2(double aoa, double A0, double cd1max_3D, double ACL1_3D, double ACD1_3D, double tbyc, double AR) {
    double ans1;
    double G1 = 2.300 * exp(-1 * pow(0.65 * tbyc,0.9));
    double G2 = 0.52 + 0.48 * exp(-1 * pow((0.65 / AR),1.1));
    double cd2max = G1 * G2;
    debug writeln("cd2max = ",cd2max);

    //(check this condition(might be typo in the paper pg18.) 
    if ((2 * A0 - ACL1_3D) < aoa && aoa < ACL1_3D) {
        
        ans1 = 0;
    } else if (aoa > ACD1_3D) {
        ans1 = cd1max_3D + (cd2max - cd1max_3D) * sin(((aoa - ACD1_3D) * 90.0 / (
                90.0 - ACD1_3D))*(PI/180.0));
    }

    return ans1;
}

class AeroDAS: AirfoilModel {
    double[] alpha; // alpha,cl, and cd from 2D airoil data file
    double[] CL;
    double[] CD;
    double AR; //Aspect ratio of wing
    double tbyc; //thickness to chord ratio
    InputParam ip;
    double alpha_query;
    double CL1; //pre-stall lift coefficient
    double CD1; //pre-stall drag coefficient
    double CL2; //post-stall lift coefficient
    double CD2; //post-stall drag coefficient

    this(double[] alpha, double[] CL, double[] CD, double tbyc, double AR = double.infinity) {
        this.alpha = alpha;
        this.CL = CL;
        this.CD = CD;
        this.AR = AR;
        this.tbyc = tbyc;

        this.ip = new InputParam(this.alpha, this.CL, this.CD);
        debug writeln("A0 is",this.ip.A0);
    }

    override double get_Cl(double alpha_query, double mach_query) {
        alpha_query *= 180.0/PI;

        //finite aspect ratio adjustment
        double ACL1_3D = this.ip.ACL1_2D + this.ip.cd1max_2D * 18.2 * pow(this.AR,-0.9);
        debug writeln("ACL1_2d = ", this.ip.ACL1_2D);
        debug writeln("ACL1_3d = ", ACL1_3D);
        
        double S1_3D = this.ip.S1_2D / (1 + this.ip.S1_2D * 18.2 * pow(this.AR,-0.9));
        double ACD1_3D = this.ip.ACD1_2D + this.ip.cl1max_2D * 18.2 * pow(this.AR,-0.9);
        double cd1max_3D = this.ip.cd1max_2D + 0.280 * pow(this.ip.cl1max_2D,2.0) * pow(this.AR,-0.9);
        double cl1max_3D = this.ip.cl1max_2D * (0.67 + 0.33 * exp(-16 / pow(this.AR,2.0)));
        // pre-stall Cl calculation
        double RCL1 = S1_3D * (ACL1_3D - this.ip.A0) - cl1max_3D; //reduction from extension of linear segment of lift curve to CL1max
        double N1 = 1 + cl1max_3D / RCL1;
        double Cl;

        if (alpha_query >= this.ip.A0) {
            CL1 = S1_3D * (alpha_query - this.ip.A0) - RCL1 * pow(
                (alpha_query - this.ip.A0) / (ACL1_3D - this.ip.A0),N1);
        } else {
            CL1 = S1_3D * (alpha_query - this.ip.A0) + RCL1 * pow(
                (this.ip.A0 - alpha_query) / (ACL1_3D - this.ip.A0),N1);
        }
        //post stall regim

        if (alpha_query >= 0) {
            CL2 = getCL2(alpha_query,this.tbyc,this.AR,ACL1_3D);
            //writeln("CL2 = ",CL2);
        } else {
            double updated_alpha_query = -alpha_query + 2.0 * this.ip.A0;
            //writeln("alpha_neg = ",updated_alpha_query);
            CL2 = -1.0 * getCL2(updated_alpha_query,this.tbyc,this.AR,ACL1_3D);
            //writeln("CL2 = ",CL2);
        }

        if (alpha_query >= this.ip.A0) {
            Cl = max(CL1, CL2);
        } else {
            Cl = min(CL1, CL2);
        }

        return Cl;
    }

    override double get_Cd(double alpha_query, double mach_query) {
        alpha_query *= 180.0/PI;

        //finite aspect ratio adjustment
        double ACL1_3D = ip.ACL1_2D + ip.cd1max_2D * 18.2 * pow(this.AR,-0.9);
        //double S1_3D = this.ip.S1_2D / (1 + this.ip.S1_2D * 18.2 * pow(this.AR,-0.9));
        double ACD1_3D = ip.ACD1_2D + ip.cl1max_2D * 18.2 * pow(this.AR,-0.9);
        double cd1max_3D = ip.cd1max_2D + 0.280 * pow(ip.cl1max_2D,2.0) * pow(AR,-0.9);
        debug writeln("CDmax_3D = ",cd1max_3D);
        
        //double cl1max_3D = this.ip.cl1max_2D * (0.67 + 0.33 * exp(-16 / pow(this.AR,2.0)));
        double Cd;

        //pre-stall Cd calculation
        if ((2 * ip.A0 - ACD1_3D) <= alpha_query && alpha_query <= ACD1_3D){
            CD1 = ip.CD0 + (cd1max_3D - ip.CD0) * pow(
                (alpha_query - ip.A0) / (ACD1_3D - ip.A0), 2);
        } else {
            CD1 = 0;
        }
        debug writeln("CD1 = ",CD1);

        if (alpha_query > (2 * ip.A0 - ACD1_3D)){
            CD2 = getCD2(alpha_query, ip.A0, cd1max_3D, ACL1_3D, ACD1_3D, tbyc, AR);
        } else {
            double updated_alpha_query = -alpha_query + 2 * this.ip.A0;
            CD2 = getCD2(updated_alpha_query, ip.A0, cd1max_3D, ACL1_3D, ACD1_3D, tbyc, AR);
        }
        debug writeln("CD2 = ",CD2);
        Cd = max(CD1, CD2);

        return Cd;
    }

    override Chunk get_Cl(Chunk alpha_query, Chunk mach_query) {
        Chunk Cl;
        foreach(i; 0..chunk_size) {
            Cl[i] = get_Cl(alpha_query[i], mach_query[i]);
        }
        return Cl;
    }

    override Chunk get_Cd(Chunk alpha_query, Chunk mach_query) {
        Chunk Cd;
        foreach(i; 0..chunk_size) {
            Cd[i] = get_Cd(alpha_query[i], mach_query[i]);
        }
        return Cd;
    }
}

auto create_aerodas_from_xfoil_polar(string filename, double tbyc) {
    /*reads 2D airfoil Data file*/
    auto file = File(filename, "r");
    debug writeln("file loaded");
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
    debug writeln("Re line is",Re);
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

    while (1) {
        double[] line1 = split(file.readln).map!(x => to!double(x)).array;
        if (line1 is null) {
            break;
        } else {
            aoa ~= line1[0];
            cl ~= line1[1];
            cd ~= line1[2];
        }
    }
    
    debug writeln("File reading completed");
    auto afmodel = new AeroDAS(aoa, cl, cd, tbyc, double.infinity);
    return afmodel;
}
