import std.stdio;
import std.file;


void main(){
    File f = File("data", "r");
    
    auto Rotor_CT = f.rawRead(new double[48][4]);
    auto Wing_CL = f.rawRead(new double[24][2]);

    writeln("wing_CL = ", Wing_CL.size);
    writeln("rotor_CT = ", Rotor_CT.size);


}


