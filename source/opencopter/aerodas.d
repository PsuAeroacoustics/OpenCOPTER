module opencopter.aerodas;

import opencopter.math;
import opencopter.liftmodel_interface;
import opencopter.memory;
import std.stdio;
import std.array;
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
// Function for CL calulation in the post stall regime
double getCL2(double aoa,double tbyc,double AR,double ACL1_3D){
    double ans;
    double F1 = 1.190 * (1 - pow(tbyc,2));
    double F2 = 0.65 + 0.35 * exp(-pow((9.0 / AR),2.3));
    double cl2max = F1 * F2;
    double RCL2 = 1.632 - cl2max;
    double N2 = 1 + cl2max / RCL2;
    if (0 < aoa && aoa < ACL1_3D){
        ans = 0;
    }
    else if (ACL1_3D <= aoa && aoa <= 92.0){
        ans = -0.032 * (aoa - 92.0) - RCL2 * pow(((92.0 - aoa) / 51.0),N2);
    }
    else if (aoa > 92){
        ans = -0.032 * (aoa - 92.0) + RCL2 * pow(((aoa - 92.0) / 51.0),N2);
    }
    return ans;
}
// Function for CD calulation in the post stall regime
double getCD2(double aoa,double A0,double cd1max_3D,double ACL1_3D,double ACD1_3D,double tbyc,double AR){
    double ans1;
    double G1 = 2.300 * exp(-1 * pow(0.65 * tbyc,0.9));
    double G2 = 0.52 + 0.48 * exp(-1 * pow((0.65 / AR),1.1));
    double cd2max = G1 * G2;
    debug writeln("cd2max = ",cd2max);
    if ((2 * A0 - ACL1_3D) < aoa && aoa < ACL1_3D)
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

class AeroDAS:Coefficent
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
        debug writeln("A0 is",this.ip.A0);
    }
        double[] getCl(Chunk alphaQuery,Chunk machQuery)
        {
            //finite aspect ratio adjustment
            double ACL1_3D = this.ip.ACL1_2D + this.ip.cd1max_2D * 18.2 * pow(this.AR,-0.9);
            debug writeln("ACL1 = ",ACL1_3D);
            double S1_3D = this.ip.S1_2D / (1 + this.ip.S1_2D * 18.2 * pow(this.AR,-0.9));
            double ACD1_3D = this.ip.ACD1_2D + this.ip.cl1max_2D * 18.2 * pow(this.AR,-0.9);
            double cd1max_3D = this.ip.cd1max_2D + 0.280 * pow(this.ip.cl1max_2D,2.0) * pow(this.AR,-0.9);
            double cl1max_3D = this.ip.cl1max_2D * (0.67 + 0.33 * exp(-16 / pow(this.AR,2.0)));
            // pre-stall Cl calculation
            double RCL1 = S1_3D * (ACL1_3D - this.ip.A0) - cl1max_3D; //reduction from extension of linear segment of lift curve to CL1max
            double N1 = 1 + cl1max_3D / RCL1;
            double[] Cl = new double[8];

            for(int i=0;i<8;i++){
                if (alphaQuery[i] >= this.ip.A0)
                {
                    CL1 = S1_3D * (alphaQuery[i] - this.ip.A0) - RCL1 * pow(
                        (alphaQuery[i] - this.ip.A0) / (ACL1_3D - this.ip.A0),N1);
                }
                else
                {
                    CL1 = S1_3D * (alphaQuery[i] - this.ip.A0) + RCL1 * pow(
                        (this.ip.A0 - alphaQuery[i]) / (ACL1_3D - this.ip.A0),N1);
                }
                //post stall regim

                if (alphaQuery[i] >= 0)
                {
                    CL2 = getCL2(alphaQuery[i],this.tbyc,this.AR,ACL1_3D);
                    //writeln("CL2 = ",CL2);
                }
                else
                {
                    double updated_alphaquery = -alphaQuery[i] + 2.0 * this.ip.A0;
                    //writeln("alpha_neg = ",updated_alphaquery);
                    CL2 = -1.0 * getCL2(updated_alphaquery,this.tbyc,this.AR,ACL1_3D);
                    //writeln("CL2 = ",CL2);
                }

                if (alphaQuery[i] >= this.ip.A0)
                {
                    Cl[i] = max(CL1, CL2);
                }
                else
                {
                    Cl[i] = min(CL1, CL2);
                }
            }
            return Cl;
        }

        double[] getCd(Chunk alphaQuery, Chunk machQuery)
        {   
            //finite aspect ratio adjustment
            double ACL1_3D = this.ip.ACL1_2D + this.ip.cd1max_2D * 18.2 * pow(this.AR,-0.9);
            double S1_3D = this.ip.S1_2D / (1 + this.ip.S1_2D * 18.2 * pow(this.AR,-0.9));
            double ACD1_3D = this.ip.ACD1_2D + this.ip.cl1max_2D * 18.2 * pow(this.AR,-0.9);
            double cd1max_3D = this.ip.cd1max_2D + 0.280 * pow(this.ip.cl1max_2D,2.0) * pow(this.AR,-0.9);
            writeln("CDmax_3D = ",cd1max_3D);
            double cl1max_3D = this.ip.cl1max_2D * (0.67 + 0.33 * exp(-16 / pow(this.AR,2.0)));
            double[] Cd = new double [8];

            for(int i=0; i<8 ; i++){
            //pre-stall Cd calculation
                if ((2 * this.ip.A0 - ACD1_3D) <= alphaQuery[i] && alphaQuery[i] <= ACD1_3D){
                    CD1 = this.ip.CD0 + (cd1max_3D - this.ip.CD0) * pow(
                        (alphaQuery[i] - this.ip.A0) / (ACD1_3D - this.ip.A0),2);
                }
                else
                {
                    CD1 = 0;
                }
                writeln("CD1 = ",CD1);

                if (alphaQuery[i] > (2 * this.ip.A0 - ACD1_3D)){
                    CD2 = getCD2(alphaQuery[i],this.ip.A0,cd1max_3D,ACL1_3D,ACD1_3D,this.tbyc,this.AR);
                }
                else{
                    double updated_alphaquery = -alphaQuery[i] + 2 * this.ip.A0;
                    CD2 = getCD2(updated_alphaquery,this.ip.A0,cd1max_3D,ACL1_3D,ACD1_3D,this.tbyc,this.AR);
                }
                debug writeln("CD2 = ",CD2);
                Cd[i]=max(CD1, CD2);
            }
            return Cd;
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
    debug writeln("File reading completed");
    auto afmodel = new AeroDAS(aoa, cl, cd, AR, tbyc);
    return afmodel;
}
