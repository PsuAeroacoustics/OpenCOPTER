module opencopter.liftmodel_interface;

interface Coefficent{
    double getCl(double alphaQuery,double machQuery);
    double getCd(double alphaQuery,double machQuery);
}