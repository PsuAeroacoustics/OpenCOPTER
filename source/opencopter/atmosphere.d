module opencopter.atmosphere;

extern (C++) struct Atmosphere {
	this(double density, double dynamic_viscosity, double speed_of_sound) {
		this.speed_of_sound = speed_of_sound;
		this.density = density;
		this.dynamic_viscosity = dynamic_viscosity;
		kinematic_viscosity = dynamic_viscosity/density;
	}
	
	immutable double density;
	immutable double dynamic_viscosity;
	immutable double kinematic_viscosity;
	immutable double speed_of_sound;
}
