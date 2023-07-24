import opencopter.math;
//import opencopter.read_afdata;
import opencopter.config;
import opencopter.airfoilmodels;
import opencopter.memory;

import std.stdio;
import std.array;
import std.math;
import std.algorithm;
import std.conv;
import std.file;
import std.range;

import plt = matplotlibd.pyplot;

void main(){
    bool haveafdata = false;

    writeln("hello");
    class DummyAF: AirfoilModel {
        string name;
        this(string _name) {
            name = _name;
        }

        Chunk get_Cl(Chunk alpha_query, Chunk mach_query) {
            Chunk c;
            return c;
        }

        Chunk get_Cd(Chunk alpha_query, Chunk mach_query) {
            Chunk c;
            return c;
        }

        double get_Cl(double alpha_query, double mach_query) {
            return 0.0;
        }

        double get_Cd(double alpha_query, double mach_query) {
            return 0.0;
        }
    }

    AirfoilModel af1 = new DummyAF("af1");
    AirfoilModel af2 = new DummyAF("af2");
    AirfoilModel af3 = new DummyAF("af3");
    AirfoilModel af4 = new DummyAF("af4");
    AirfoilModel af5 = new DummyAF("af5");
    AirfoilModel af6 = new DummyAF("af6");
    AirfoilModel af7 = new DummyAF("af7");

    auto af_models = [af1, af2, af3, af4, af5, af6, af7];

    size_t[2][] extents = [
        [0, 9],
        [10, 21],
        [22, 25],
        [26,28],
        [29,30],
        [31, 41],
        [42, 47]
    ];

    auto blade_af = new BladeAirfoil(af_models, extents);

    writeln("blade_af.airfoil_models.length: ", blade_af.airfoil_models.length);
    foreach(blade_af_models; blade_af.airfoil_models) {
        writeln("blade_af_models.length: ", blade_af_models.length);
        foreach(af_model; blade_af_models) {
            writeln("af_model.name: ", af_model.to!DummyAF.name);
        }
        writeln;
    }

    if(haveafdata){
        //auto file= File("vsppropaf.C81");
        writeln("loading file");
        /*string test_line = "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8";
        double[] split_line = test_line.split.map!(x => to!double(x)).array;
        writeln(split_line[2]);*/
        auto NACA0012 = load_c81_file("vsppropaf[1561].inp");
        double[] aoa_array = [-90.0,-60.0,-30.0,0.0,1.0,1.5,2.0,2.5,3.0,3.5,4.5,5.0,10.0,20.0,40.0,60.0,90.0,120.0,180.0];
        double[] mach_array = [0.4,0.45,0.5];
        double desired_aoa = 5;
        double desired_mach = 0.45;
        double[][] Cl_array = new double[][](aoa_array.length,mach_array.length);
        writeln("lenght of AOA array is : ",aoa_array.length);
        writeln("length od Cl arrat is : ",Cl_array.length);
        double[][] Cd_array = new double[][](aoa_array.length,mach_array.length);
        double[][] Cm_array = new double[][](aoa_array.length,mach_array.length);
        /*writeln("Getting cl");
        double desired_CL = NACA0012.get_Cl(desired_aoa,desired_mach);
        writeln("Getting cd");
        double desired_CD = NACA0012.get_Cd(desired_aoa,desired_mach);
        writeln("Getting cml");
        double desired_CM = NACA0012.get_Cm(desired_aoa,desired_mach);

        writeln("Cl =", desired_CL);
        writeln("Cd =", desired_CD);
        writeln("Cm =", desired_CM);*/

        for(int i=0; i<mach_array.length;i++){
            for(int j=0; j<aoa_array.length;j++){
                Cl_array[j][i] = NACA0012.get_Cl(aoa_array[j],mach_array[i]);
                Cd_array[j][i] = NACA0012.get_Cd(aoa_array[j],mach_array[i]);
                Cm_array[j][i] = NACA0012.get_Cm(aoa_array[j],mach_array[i]);
            }
        }

        /*auto coefficient = File("airfoil_data.txt","w");
        coefficient.write(to!(string)(aoa_array),"\n",to!(string)(mach_array));
        for(int i=0;i<Cl_array.length;i++){
            //coefficient.writeln((split(to!(string)(Cl_array[0..$][i])),","));
        }*/
        
        double[] Cl1 = Cl_array[0..$][0];
        writeln(Cl1.length);
        double[] Cl2 = Cl_array[2][0..$];
        writeln(Cl2.length);
        double[] Cl3 = Cl_array[1][0..$];
        plt.plot(aoa_array,Cl_array,["label": "$Cl$"]);
        plt.plot(aoa_array,Cd_array,["label": "$Cd$"]);
        plt.plot(aoa_array,Cm_array,["label": "$Cm$"]);
        plt.xlabel("Angle of attack (degrees)");
        plt.ylabel("Coefficient of Lift");
        plt.legend();
        plt.show();
        plt.grid();
        plt.savefig("Cl.png");
    }
    else{
        
        //auto afdata2D = create_aerodas_from_xfoil_polar("NACA23012_50000_polar.dat", 0.12);
        //auto afdata2D = create_aerodas_from_xfoil_polar("rc410_50000_polar.dat", 0.10);
        auto afdata2D = create_aerodas_from_xfoil_polar("sm701_Re250k.dat", 0.16);

        //writeln(afdata2D.CL);
        //double[110] aoa_array;
        //auto aoa_array = iota(-110*PI/180.0, 110*PI/180.0, 1*PI/180.0).array;
        auto aoa_array = iota(-20*PI/180.0, 110*PI/180.0, 1*PI/180.0).array;

        writeln("aoa_array.length: ", aoa_array.length);

        double mach = 0.5;
        //for (int i=0;i<110;i++){aoa_array[i]=(i-10.0)*(PI/180.0);}
        double[] CL_array = new double[aoa_array.length];
        double[] CD_array = new double[aoa_array.length];
        for(int i=0;i<aoa_array.length;i++){
            CL_array[i] = afdata2D.get_Cl(aoa_array[i],mach);
            CD_array[i] = afdata2D.get_Cd(aoa_array[i],mach);
        }

        plt.figure();
        plt.plot(aoa_array.map!(a => a*180.0/PI),CL_array,["label": "$AERODAS MODEL$"]);
        plt.plot(afdata2D.alpha,afdata2D.CL,["label": "XFOIL data"]);
        plt.xlabel("Angle of attack (degrees)");
        plt.ylabel("Coefficient of lift");
        plt.legend();
        plt.grid();
        plt.savefig("Cl.png", ["dpi": 500]);

        plt.clear();
        plt.figure();
        plt.plot(aoa_array.map!(a => a*180.0/PI),CD_array,["label": "$AERODAS MODEL$"]);
        plt.plot(afdata2D.alpha,afdata2D.CD,["label": "XFOIL data"]);
        plt.xlabel("Angle of attack (degrees)");
        plt.ylabel("Coefficient of drag");
        plt.legend();
        plt.grid();
        plt.savefig("Cd.png", ["dpi": 500]);
        /*double aoa = 30.0;
        double CL = afdata2D.get_Cl(aoa);
        double CD = afdata2D.get_Cd(aoa);
        writeln("CL = ", CL);
        writeln("Cd = ", CD);*/

    }
    


}
