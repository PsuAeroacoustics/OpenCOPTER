module opencopter.liftmodel_interface;
import opencopter.memory;

interface Coefficent{
    double[] getCl(Chunk alphaQuery,Chunk machQuery);
    double[] getCd(Chunk alphaQuery,Chunk machQuery);
}