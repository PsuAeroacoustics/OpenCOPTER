module opencopter.inflow.huangpeters;

import opencopter.aircraft;
import opencopter.config;
import opencopter.inflow;
import opencopter.math;
import opencopter.math.blas;
import opencopter.math.lapacke;
import opencopter.memory;
import opencopter.wake;

import numd.calculus.integration.forwardeuler;
import numd.calculus.integration.rk4;
import numd.utility : linspace;

import std.array;
import std.algorithm;
import std.complex;
import std.conv;
import std.math;
import std.range;
import std.stdio;
import std.typecons;

@nogc private double H(long m, long n) {
	return double_factorial(n + m - 1)*double_factorial(n - m - 1)/(double_factorial(n + m)*double_factorial(n - m));
}

@nogc private double K(long m, long n) {
	return ((PI/2.0)^^((-1.0)^^(n.to!double + m.to!double)))*H(m, n);
}

private void zero_matrix(ref double[][] M_c) {
	foreach(ref _M; M_c) {
		_M[] = 0;
	}
}

unittest {
	writeln("Hello hp");
	import numd.utility : linspace;

	double[] generate_radius_points(size_t n_sections) {
		auto i = linspace(0., (n_sections-1).to!double, n_sections);
		auto r_array = new double[n_sections];
		r_array[] = (1./n_sections)*(i[] + 0.5);
		return r_array;
	}

	//import plt = matplotlibd.pyplot;

	immutable double AR = 9.2204;
	immutable C_T = 0.008;
	immutable size_t elements = 48;
	immutable size_t num_blades = 2;
	immutable double R = 1;
	immutable double d_azimuth = (2.0*PI)/num_blades.to!double;

	immutable double theta_75 = 3*(PI/180.0);

	auto r = generate_radius_points(elements);
	double[] C_l_alpha = new double[elements];
	C_l_alpha[] = 2.0*PI;
	double[] c = new double[elements];
	c[] = 1/AR;
	double[] twist = new double[elements];
	double[] sweep = new double[elements];
	sweep[] = 0;

	immutable theta_tw_1 = -10.0*(PI/180.0);
	twist[] = theta_75 + (r[] - 0.75)*theta_tw_1;
	
	double[] alpha_0 = new double[elements];
	alpha_0[] = 0;
	

	auto rotor = RotorGeometry(
			num_blades,
			Vec3(0, 0, 0),
			R,
			0
		);

	foreach(b_idx, ref blade; rotor.blades) {
		blade.chunks = new BladeGeometryChunk[elements/chunk_size];
		blade.set_geometry_array!"r"(r);
		blade.set_geometry_array!"twist"(twist);
		blade.set_geometry_array!"sweep"(sweep);
		blade.set_geometry_array!"alpha_0"(alpha_0);
		blade.set_geometry_array!"chord"(c);
		blade.set_geometry_array!"C_l_alpha"(C_l_alpha);
		blade.azimuth_offset = b_idx*d_azimuth;
		import std.algorithm : sum;
		blade.average_chord = R*c.sum/c.length.to!double;

	}

	void simple_harmonic_30_ondisk() {
		double mu = 0.03651;
		double mu_z = 0;
		auto huang_peters = new HuangPetersInflow(12, 8, &rotor, mu);
		//auto huang_peters = new HuangPetersInflowT(6, 4, elements, mu);

		//writeln("tot states: ", huang_peters.total_states);

		huang_peters.tau_c[] = 0;
		huang_peters.tau_c[0] = abs(mu);//.03651;
		huang_peters.average_inflow = sqrt(C_T/2.0);
		huang_peters.simple_harmonic_solution(mu, mu_z);

		//writeln(huang_peters.alpha);
		//writeln;
		//writeln(huang_peters.a);

		auto x = linspace(-2.1, 2.1, 128);
		double[] v = new double[x.length];
		foreach(idx, ref x_chunk; x.chunks(chunk_size).enumerate) {
			immutable Chunk _x = x_chunk[].map!(a => a.to!double).staticArray!Chunk;
			immutable Chunk y = 0;
			immutable Chunk x_e = 0;
			immutable Chunk z = 0;

			immutable inflow = huang_peters.inflow_at(_x, y, z, x_e, 0.0);
			v[idx*chunk_size..idx*chunk_size + chunk_size] = inflow[];
		}

		double[] huang_x = [-1.9733103745157128, -1.9113215669393029, -1.8415841584158419, -1.7907877744296172, -1.7270770555316404, -1.6771416272061992, -1.6117089969866554, -1.5531640120533796, -1.4946190271201036, -1.442100731812312, -1.3732242789496345, -1.3241498062849768, -1.2569952647438662, -1.196728368489023, -1.139905294877314, -1.095996556177357, -1.0710288420146366, -1.0417563495479985, -1.0348687042617306, -1.0245372363323288, -1.0142057684029278, -1.000430477830392, -0.9900990099009903, -0.9789065863108053, -0.969436074042187, -0.958243650452002, -0.9444683598794665, -0.9306930693069309, -0.9040034438226432, -0.8859233749461906, -0.8592337494619029, -0.8368489022815326, -0.8153250107619461, -0.771416272061989, -0.737839001291434, -0.7059836418424454, -0.6482996125699527, -0.5931984502798109, -0.5277658200602666, -0.4769694360740424, -0.40378820490744816, -0.35643564356435675, -0.2849763237193277, -0.23504089539388806, -0.1661644425312101, -0.115368058544985, -0.0490744726646577, 0.005165733964700436, 0.06887645286267752, 0.12828239345673698, 0.1902712010331471, 0.22126560482135194, 0.25828669823504047, 0.3021954369349973, 0.34179939733103737, 0.373654756780025, 0.4175634954799823, 0.4571674558760228, 0.49591046061127786, 0.5389582436504519, 0.578562204046492, 0.619027120103314, 0.6586310804993536, 0.6896254842875589, 0.7188979767541963, 0.74214377959535, 0.7722772277227725, 0.7998278088678425, 0.8291003013344813, 0.8514851485148514, 0.8695652173913029, 0.889367197589324, 0.911752044769695, 0.92811020232458, 0.9393026259147654, 0.952216960826517, 0.9642703400774857, 0.9746018080068866, 0.9840723202755051, 0.9909599655617738, 0.9995695221696086, 1.0279810589754628, 1.050365906155832, 1.083082221265605, 1.1459319845027967, 1.1993112354713729, 1.2613000430477825, 1.3198450279810587, 1.391304347826086, 1.4429616874730948, 1.5187257856220402, 1.560912613000431, 1.6375376668101587, 1.6934997847610833, 1.7503228583727934, 1.8105897546276357, 1.876883340507964, 1.9207920792079203, 1.995695221696082];
		double[] huang_y = [0.15108795586883206, 0.1621207477781179, 0.1768311369904978, 0.1878639288997852, 0.20349371743794054, 0.22004290530186954, 0.24026969046889324, 0.2623352742874654, 0.28715905608335834, 0.31290223720502564, 0.35243640821330025, 0.38461538461538436, 0.4351823475329446, 0.4903463070793743, 0.555623659209316, 0.6218204106650316, 0.6733067729083662, 0.7505363162733676, 0.7873122893043205, 0.8553478394115837, 0.9233833895188461, 1.0153233220962299, 1.0851976708550408, 1.1513944223107564, 1.2185105730922459, 1.2856267238737353, 1.341710082745939, 1.393196444989273, 1.475022984983144, 1.5191541526202879, 1.5743181121667174, 1.6166104811523134, 1.6487894575543973, 1.7076310144039222, 1.7471651854121966, 1.7784247624885072, 1.8225559301256509, 1.8519767085504133, 1.8722034937174374, 1.8841556849524972, 1.8878332822555923, 1.885994483604045, 1.874961691694759, 1.8602513024823777, 1.8299111247318414, 1.804167943610174, 1.762794973950352, 1.724180202267851, 1.675452038001838, 1.6239656757585035, 1.5651241189089786, 1.5320257431811208, 1.4906527735212987, 1.4364082133006426, 1.3858412503830821, 1.3426294820717122, 1.279190928593318, 1.2203493717437932, 1.1587496169169467, 1.0851976708550402, 1.0153233220962294, 0.9371743794054541, 0.8562672387373573, 0.7873122893043201, 0.7192767391970569, 0.6604351823475318, 0.5822862396567565, 0.5078148942690763, 0.4213913576463364, 0.35335580753907303, 0.2926754520380006, 0.21912350597609498, 0.13086117070180725, 0.05547042598835317, -0.005209929512719036, -0.07508427827153019, -0.15047502298498383, -0.22954336500153305, -0.3030953110634398, -0.36193686791296464, -0.4327306159975497, -0.36377566656451243, -0.320563898253142, -0.28102972724486763, -0.220349371743795, -0.18265399938706883, -0.15139442231075995, -0.12749003984063867, -0.1091020533251621, -0.09806926141587535, -0.08519767085504348, -0.07784247624885365, -0.07048728164266116, -0.06589028501379213, -0.0603738890591492, -0.05669629175605362, -0.051179895801411135, -0.050260496475636796, -0.046582899172541214];

		// plt.figure;
		// plt.plot(x, v, "-", huang_x.map!(a => -a), huang_y, ".");
		// plt.legend(["opencopter", "Huang Dissertation"]);
		// // wake skew 30
		// plt.xlim(-2, 2);
		// plt.ylim(-1, 2);
		// plt.xlabel("$x$");
		// plt.ylabel("$v_z$");
		// plt.title("Simple harmonic solution\\nOn disk\\nWake skew = $"~(std.math.abs(huang_peters.wake_skew)*(180.0/PI)).to!string~" ^\\circ$");
		//plt.show;
	}

	void simple_harmonic_60_ondisk() {
		double mu = -4.*0.02738;
		double mu_z = 0;
		auto huang_peters = new HuangPetersInflow(12, 8, &rotor, mu);

		huang_peters.tau_c[] = 0;
		huang_peters.tau_c[0] = abs(mu);//.11;
		huang_peters.average_inflow = sqrt(C_T/2.0);
		huang_peters.simple_harmonic_solution(mu, mu_z);

		auto x = linspace(-2., 2., 256);
		double[] v = new double[x.length];
		foreach(idx, ref x_chunk; x.chunks(chunk_size).enumerate) {
			immutable Chunk _x = x_chunk[].map!(a => a.to!double).staticArray!Chunk;
			immutable Chunk y = 0;
			immutable Chunk x_e = 0;
			immutable Chunk z = 0;

			immutable inflow = huang_peters.inflow_at(_x, y, z, x_e, 0.0);
			v[idx*chunk_size..idx*chunk_size + chunk_size] = inflow[];
		}

		double[] huang_x = [-1.9965285311347365, -1.9505315686699938, -1.9106096767194622, -1.8437839010631376, -1.7960512041657628, -1.723150357995227, -1.684964200477327, -1.612931221523107, -1.5669342590583644, -1.5105228899978305, -1.461054458667824, -1.4185289650683446, -1.3681926665220228, -1.327402907355175, -1.2970275547841181, -1.261444998915166, -1.2180516380993711, -1.164243870687785, -1.1225862443046217, -1.0757214146235627, -1.040138858754611, -1.0288565849425042, -1.0149707094814502, -1.0054241701019746, -0.9993490995877632, -0.9785202863961817, -0.9594272076372317, -0.9221089173356478, -0.8735083532219574, -0.8361900629203733, -0.7650249511824694, -0.7164243870687788, -0.6400520720329794, -0.596658711217184, -0.5254935994792802, -0.4768930353655896, -0.403992189195052, -0.35452375786504664, -0.2894337166413541, -0.2512475591234544, -0.20091126057713193, -0.15838576697765205, -0.11325667172922582, -0.053373833803428194, -0.02299848123237158, 0.011716207420264357, 0.06465610761553453, 0.10110653070080255, 0.12974614883922708, 0.1818181818181812, 0.21913647211976572, 0.2538511607724012, 0.3050553265350402, 0.3449772184855715, 0.37448470384031163, 0.4222174007376869, 0.4586678238229549, 0.49251464525927524, 0.5428509438055964, 0.5784334996745493, 0.6148839227598168, 0.6643523540898233, 0.7060099804729871, 0.7372532002603589, 0.7797786938598397, 0.8170969841614228, 0.8596224777609018, 0.8978086352788002, 0.9229767845519632, 0.94901280104144, 0.9646344109351266, 0.9785202863961806, 0.9993490995877621, 1.018442178346711, 1.0297244521588187, 1.040138858754609, 1.0722499457582972, 1.0965502278151442, 1.15035799522673, 1.1954870904751567, 1.2692558038620088, 1.330006509004122, 1.3942286830114985, 1.4515079192883484, 1.5244087654588845, 1.588630939466261, 1.6615317856367966, 1.7795617270557598, 1.897591668474723, 1.9947927967021037];
		double[] huang_y = [0.6336939721792891, 0.6692426584234936, 0.6986089644513136, 0.7527047913446685, 0.7944358578052553, 0.8655332302936634, 0.9057187017001547, 0.990726429675425, 1.0510046367851622, 1.1313755795981457, 1.2086553323029374, 1.2782071097372492, 1.366306027820711, 1.4435857805255028, 1.502318392581144, 1.5718701700154565, 1.6568778979907268, 1.7604327666151474, 1.8346213292117468, 1.9165378670788258, 1.9814528593508502, 2.0803709428129835, 2.2040185471406524, 2.2859350850077282, 2.335394126738795, 2.398763523956724, 2.4544049459041735, 2.5069551777434316, 2.5347758887171565, 2.54095826893354, 2.5332302936630606, 2.5162287480680066, 2.4744976816074193, 2.4451313755795985, 2.387944358578053, 2.343122102009274, 2.268933539412674, 2.21483771251932, 2.1391035548686252, 2.091190108191654, 2.030911901081917, 1.9675425038639884, 1.9010819165378676, 1.8129829984544057, 1.7666151468315305, 1.7125193199381765, 1.6259659969088103, 1.5641421947449772, 1.519319938176198, 1.4281298299845444, 1.3601236476043281, 1.2967542503863994, 1.1978361669242656, 1.1221020092735707, 1.061823802163833, 0.9613601236476041, 0.8840803709428129, 0.8098918083462134, 0.6986089644513136, 0.6151468315301392, 0.5255023183925811, 0.40340030911901126, 0.2952086553323032, 0.2055641421947456, 0.08500772797527079, -0.02472952086553315, -0.16074188562596525, -0.29366306027820777, -0.3925811437403395, -0.5100463678516225, -0.5842349304482224, -0.6769706336939723, -0.7727975270479135, -0.6707882534775886, -0.599690880989181, -0.5162287480680061, -0.4574961360123644, -0.4142194744976835, -0.32766615146831546, -0.2673879443585787, -0.19629057187016974, -0.1545595054095834, -0.1236476043276662, -0.1035548686244212, -0.08500772797526901, -0.07573415765069447, -0.06336939721792811, -0.05255023183925722, -0.04482225656877992, -0.038639876352396296];

		// plt.figure;
		// plt.plot(x, v, "-", huang_x, huang_y, ".");
		// plt.legend(["opencopter", "Huang Dissertation"]);
		// plt.xlim(-2, 2);
		// plt.ylim(-2, 3);
		// plt.xlabel("$x$");
		// plt.ylabel("$v_z$");
		// //plt.title("$\\chi$ = "~(std.math.abs(huang_peters.wake_skew)*(180.0/PI)).to!string);
		// plt.title("Simple harmonic solution\\nOn disk\\nWake skew = $"~(std.math.abs(huang_peters.wake_skew)*(180.0/PI)).to!string~" ^\\circ$");
		//plt.show;
	}

	void simple_harmonic_85_above_disk() {
		double mu = 26.35*0.0274;
		double mu_z = 0.0;
		//double mu = -45*0.0274;
		//auto huang_peters = new HuangPetersInflow(12, 8, &rotor, mu);
		auto huang_peters = new HuangPetersInflow(6, 4, &rotor, mu);

		//writeln("tot states: ", huang_peters.total_states);
		huang_peters.tau_c[] = 0;
		huang_peters.tau_c[0] = abs(mu);//.724;
		huang_peters.average_inflow = sqrt(C_T/2.0);

		huang_peters.times = linspace(0., 50, huang_peters.time_history);
		huang_peters.curr_state = 0;//huang_peters.time_history - 10;

		huang_peters.simple_harmonic_solution(mu, mu_z);

		foreach(idx; 0..huang_peters.time_history - 1) {
			huang_peters.state_history[idx + 1][] = huang_peters.state_history[0][];
		}

		huang_peters.curr_state = huang_peters.time_history - 10;

		/+huang_peters.simple_harmonic_solution(mu, mu_z);

		foreach(idx; 0..huang_peters.time_history - 1) {
			huang_peters.state_history[idx + 1][] = huang_peters.state_history[0][];
		}+/

		

		auto x = linspace(-4.1, 4.1, 1024);
		double[] v = new double[x.length];
		foreach(idx, ref x_chunk; x.chunks(chunk_size).enumerate) {
			immutable Chunk _x = x_chunk[].map!(a => a.to!double).staticArray!Chunk;
			immutable Chunk y = 0;
			immutable Chunk x_e = 0;
			immutable Chunk z = -0.4;

			immutable inflow = huang_peters.inflow_at(_x, y, z, x_e, 0.0);
			v[idx*chunk_size..chunk_size*(idx + 1)] = inflow[];
		}

		double[] huang_x = [-1.9912910951447855, -1.943392118441106, -1.8928804702808621, -1.847594165033747, -1.8040496407576745, -1.7535379925974308, -1.70476812540823, -1.655998258219029, -1.597648595689092, -1.5532331809274984, -1.4931417374265186, -1.4417591987807534, -1.3755715218811237, -1.3146091878946224, -1.247550620509471, -1.1944263008926628, -1.1230132810799043, -1.0812105377748749, -0.9967341606792948, -0.9523187459177009, -0.925321140866536, -0.8774221641628566, -0.8042673633790554, -0.7554974961898544, -0.6945351622033531, -0.6405399521010235, -0.5830611800566079, -0.5360330938384501, -0.48116699325059864, -0.4263008926627474, -0.3670803396472895, -0.308730677117353, -0.27389505769649514, -0.23383409536250843, -0.18506422817330792, -0.14413237535379997, -0.10668408447637656, -0.06139777922926104, -0.01872414543871148, 0.014369693011103468, 0.05965599825821899, 0.09971696059220525, 0.13803614195514857, 0.18593511865882828, 0.22164162856520875, 0.2599608099281503, 0.3078597866318309, 0.34356629653821, 0.3810145874156321, 0.4323971260613986, 0.46984541693881976, 0.5046810363596772, 0.5543217940343999, 0.5969954278249507, 0.6266057043326798, 0.676246462007402, 0.713694752884825, 0.7459177008491169, 0.799912910951448, 0.8356194208578267, 0.8678423688221195, 0.914870455040278, 0.9540605268887434, 1.0341824515567177, 1.09688656651426, 1.1796211626387976, 1.244067058567385, 1.3119965164380578, 1.3799259743087302, 1.4704985848029604, 1.5384280426736332, 1.6028739386022206, 1.6690616155018496, 1.7247986065752219, 1.7857609405617234, 1.8536903984323958, 1.9233616372741125, 1.9947746570868707];
		double[] huang_y = [1.3589147286821703, 1.383720930232558, 1.4100775193798447, 1.4341085271317826, 1.457364341085271, 1.4837209302325578, 1.5116279069767438, 1.5372093023255813, 1.5689922480620155, 1.5914728682170538, 1.6217054263565889, 1.6472868217054262, 1.6759689922480618, 1.69922480620155, 1.721705426356589, 1.7348837209302321, 1.7457364341085269, 1.748062015503876, 1.7403100775193796, 1.7294573643410849, 1.71937984496124, 1.7093023255813948, 1.6852713178294572, 1.6635658914728682, 1.629457364341085, 1.595348837209302, 1.5511627906976742, 1.511627906976744, 1.4620155038759688, 1.409302325581395, 1.3472868217054264, 1.2821705426356584, 1.242635658914728, 1.1968992248062014, 1.1379844961240306, 1.0891472868217051, 1.042635658914728, 0.9860465116279065, 0.9317829457364337, 0.8883720930232555, 0.8286821705426353, 0.7775193798449611, 0.7279069767441857, 0.6651162790697671, 0.6186046511627903, 0.5689922480620149, 0.5077519379844959, 0.46201550387596857, 0.4147286821705425, 0.3519379844961237, 0.30697674418604626, 0.265116279069767, 0.20930232558139505, 0.16279069767441823, 0.1325581395348836, 0.08372093023255767, 0.04961240310077475, 0.023255813953488857, -0.01705426356589168, -0.04031007751938054, -0.058139534883721034, -0.081395348837209, -0.09457364341085306, -0.11317829457364459, -0.11627906976744251, -0.11240310077519444, -0.10387596899224816, -0.0922480620155044, -0.07984496124031093, -0.06511627906976791, -0.05581395348837237, -0.048062015503876676, -0.04186046511627861, -0.03720930232558217, -0.032558139534884845, -0.02790697674418663, -0.024806201550388263, -0.022480620155039155];

		// plt.figure;
		// plt.plot(x, v, "-", huang_x, huang_y, ".");
		// plt.legend(["opencopter", "Huang Dissertation"]);
		// plt.xlim(-2, 2);
		// plt.ylim(-0.5, 5);
		// plt.xlabel("$x$");
		// plt.ylabel("$v_z$");
		// plt.title("Simple harmonic solution\\n0.4 radii above disk\\nWake skew = $"~(std.math.abs(huang_peters.wake_skew)*(180.0/PI)).to!string~" ^\\circ$");
		//plt.show;
	}

	void simple_harmonic_85_below_disk() {
		double mu = 26.35*0.0274;
		double mu_z = 0.0;
		//double mu = 0.0000;
		//double mu = -45*0.0274;
		//auto huang_peters = new HuangPetersInflow(12, 8, &rotor, mu);
		//auto huang_peters = new HuangPetersInflow(10, 6, &rotor, mu);
		auto huang_peters = new HuangPetersInflow(6, 4, &rotor, mu);

		//auto z_idx = huang_peters.find_z_bracket(10);
		//writeln(z_idx);
		//writeln(huang_peters.contraction_z_array_alias);
		huang_peters.tau_c[] = 0;
		huang_peters.tau_c[0] = abs(mu);//.71;
		huang_peters.average_inflow = sqrt(C_T/2.0);
		debug writeln("huang_peters.average_inflow: ", huang_peters.average_inflow);
		huang_peters.simple_harmonic_solution(mu, mu_z);
		debug writeln("huang_peters.average_inflow: ", huang_peters.average_inflow);

		debug writeln(huang_peters.state_history[0][]);
		foreach(idx; 0..huang_peters.time_history - 1) {
			huang_peters.state_history[idx + 1][] = huang_peters.state_history[0][];
		}

		huang_peters.times = linspace(0., 50, huang_peters.time_history);
		huang_peters.curr_state = huang_peters.time_history - 10;

		//auto x = linspace(-2.1, 2.1, 256);
		auto x = linspace(-6.1, 2.1, 4*4096);
		double[] v = new double[x.length];

		debug writeln("huang_peters.average_inflow: ", huang_peters.average_inflow);
		foreach(idx, ref x_chunk; x.chunks(chunk_size).enumerate) {
			immutable Chunk z = 0.4;
			immutable Chunk x_offset = z[]*tan(huang_peters.chi);
			immutable Chunk _x = x_chunk[].map!(a => a.to!double).staticArray!Chunk[] - x_offset[];
			immutable Chunk y = 0;
			immutable Chunk x_e = 0;
			

			immutable inflow = huang_peters.inflow_at(_x, y, z, x_e, 0.0);
			v[idx*chunk_size..idx*chunk_size + chunk_size] = inflow[];
		}

		double[] huang_x = [-1.996155694377703, -1.917347429120615, -1.8443056222969725, -1.7520422873618453, -1.6751561749159056, -1.5828928399807785, -1.4560307544449782, -1.4022104757328206, -1.3176357520422872, -1.2369053339740508, -1.1869293608841902, -1.1331090821720324, -1.0831331090821719, -1.0254685247477173, -0.9985583853916387, -0.9639596347909658, -0.8870735223450263, -0.8294089380105716, -0.7602114368092263, -0.7179240749639595, -0.64103796251802, -0.5987506006727539, -0.5007208073041804, -0.37962518020182534, -0.26045170591061995, -0.13743392599711646, -0.033637674195098555, 0.09899086977414528, 0.22200864968765144, 0.3411821239788573, 0.46227775108121083, 0.5776069197501204, 0.7006246996636247, 0.8217203267659783, 0.9389716482460355, 0.9870254685247479, 1.0158577606919752, 1.061989428159539, 1.1081210956271024, 1.1696299855838541, 1.2465160980297938, 1.315713599231139, 1.3887554060547815, 1.4348870735223453, 1.5194617972128786, 1.5905814512253724, 1.665545410860163, 1.7616530514175879, 1.830850552618934, 1.892359442575685, 1.957712638154733, 1.994233541566555];

		double[] huang_y = [2.3801369863013706, 2.4568493150684936, 2.5335616438356174, 2.629452054794521, 2.7157534246575343, 2.821232876712329, 2.9650684931506857, 3.0178082191780824, 3.089726027397261, 3.1280821917808224, 3.1280821917808224, 3.1089041095890413, 3.07054794520548, 3.0417808219178086, 3.051369863013699, 3.1089041095890413, 3.1952054794520555, 3.238356164383562, 3.2863013698630144, 3.305479452054795, 3.3438356164383567, 3.3630136986301378, 3.396575342465754, 3.425342465753425, 3.449315068493151, 3.458904109589042, 3.468493150684932, 3.463698630136987, 3.454109589041096, 3.4349315068493156, 3.406164383561644, 3.3726027397260276, 3.3198630136986305, 3.2479452054794526, 3.147260273972603, 3.065753424657535, 3.0369863013698635, 3.0561643835616445, 3.099315068493151, 3.1328767123287675, 3.1280821917808224, 3.0945205479452063, 3.0369863013698635, 2.9890410958904114, 2.902739726027398, 2.821232876712329, 2.7349315068493154, 2.629452054794521, 2.5527397260273976, 2.49041095890411, 2.4280821917808226, 2.389726027397261];

		// plt.figure;
		// plt.plot(x, v, "-", huang_x.map!(a => -a), huang_y, ".");
		// plt.legend(["opencopter", "Huang Dissertation"]);
		// //plt.xlim(-2, 2);
		// plt.xlim(-6, 2);
		// plt.ylim(-1, 6);
		// plt.xlabel("$x$");
		// plt.ylabel("$v_z$");
		// plt.title("Simple harmonic solution\\n0.4 radii below disk\\nWake skew = $"~(std.math.abs(huang_peters.wake_skew)*(180.0/PI)).to!string~" ^\\circ$");
		//plt.show;
	}

	//simple_harmonic_30_ondisk;
	//simple_harmonic_60_ondisk;
	//simple_harmonic_85_above_disk;
	//simple_harmonic_85_below_disk;
	//plt.show;
}

@nogc private Chunk pow(immutable Chunk x, long power) {
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	Chunk res = 1;
	if(power == 0) {
		res[] = 1;
	} else {
		foreach(idx; 0..power) {
			res[] *= x[];
		}
	}
	return res;
}

@nogc private Chunk associated_legendre_polynomial(bool reduce_order = false)(long m, long n, Chunk x, double[] c) {
	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	//x[] = 2.0*x[] - 1.0;

	Chunk p;
	p[] = 0;
	foreach(k; m..n + 1) {
		p[] += c[k - m + 1]*pow(x[], k - m)[];
	}

	immutable Chunk one_m_x2 = 1.0 - x[]*x[];
	Chunk radicand = pow(one_m_x2, m/2)[];
	if((m & 0x1) != 0) {
		// m is odd so we have an integer power plus a square root
		radicand[] *= sqrt(one_m_x2)[];
	}
	p[] *= c[0]*radicand[];
	return p;
}

@nogc private Chunk associated_legendre_polynomial_nh(bool reduce_order = false)(long r, long j, immutable Chunk x, double[] c) {

	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	//immutable Chunk x = 2.0*_x[] - 1.0;
	Chunk p = 0;
	Chunk r_bar = 1.0 - x[]*x[];
	r_bar = sqrt(r_bar);
	foreach(qidx, q; iota(r, j, 2).enumerate) {
		p[] += c[q - r + 1]*pow(r_bar[], q)[];
	}
	p[] *= c[0];
	return p;
}

@nogc private void associated_legendre_function(Chunk x, ref Chunk[][] Qmn_bar, double[][] K_table) {

	version(LDC) pragma(inline, true);
	version(GNU) pragma(inline, true);

	immutable Chunk xs = x[]*x[];

	immutable Chunk one_plus_xs = 1.0 + xs[];
	immutable Chunk sqrt_opxs = sqrt(one_plus_xs);
	immutable Chunk one_over_sqrt_opxs = 1.0/sqrt_opxs[];

	Qmn_bar[0][0][] = (2.0/PI)*((PI/2.0) - atan(x)[])[];
	Qmn_bar[0][1][] = 1.0 - x[]*PI/2.0*Qmn_bar[0][0][];

	foreach(n; 1..Qmn_bar[0].length - 1) {
		immutable long m = 0;
		immutable Kc1 = (2.0*n.to!double + 1.0)*K_table[m][n];
		
		Qmn_bar[0][n + 1][] = Qmn_bar[0][n - 1][] - Kc1*x[]*Qmn_bar[0][n][];
	}

	foreach(m; 0..Qmn_bar.length - 1) {
		foreach(n; 1..Qmn_bar[m + 1].length) {

			immutable Knmm = (n.to!double - m.to!double)*K_table[m][n];
			immutable Chunk Qtmp = Knmm*x[]*Qmn_bar[m][n][];

			immutable Chunk Qmn_delta = Qmn_bar[m][n - 1][] - Qtmp[];

			Qmn_bar[m + 1][n][] = Qmn_delta[];
			//Qmn_bar[m + 1][n][] = Qmn_delta[]*one_over_sqrt_opxs[];
		}
	}
	foreach(m; 0..Qmn_bar.length - 1) {
		foreach(n; 1..Qmn_bar[m + 1].length) {
			foreach(_; 1..m) {
				Qmn_bar[m + 1][n][] *= one_over_sqrt_opxs[];
			}
		}
	}

}

@nogc Chunk sign(Chunk x) {
	Chunk res;
	foreach(idx, ref _x; x) {
		if(_x > 0.0) {
			res[idx] = 1.0;
		} else {
			res[idx] = -1.0;
		}
	}
	return res;
}

@nogc double sign(double x) {
	if(x > 0.0) {
		return 1.0;
	} else {
		return -1.0;
	}
}

private struct ElipticalCoords {
	Chunk nu;
	Chunk eta;
	Chunk psi;
}

private struct CartisianCoords {
	Chunk x;
	Chunk y;
	Chunk z;
}

@nogc private auto to_eliptical(immutable CartisianCoords coords) {
	immutable Chunk S = coords.x[]*coords.x[] + coords.y[]*coords.y[] + coords.z[]*coords.z[];
	immutable Chunk one_m_s = (1.0 - S[]);
	immutable Chunk tmp1 = one_m_s[]*one_m_s[] + 4.0*coords.z[]*coords.z[];
	immutable Chunk tmp2 = opencopter.math.sqrt(tmp1);
	immutable Chunk nu_tmp = one_m_s[] + tmp2[];
	immutable Chunk eta_tmp = -one_m_s[] + tmp2[];
	//immutable Chunk nu_tmp = 1.0 - S[] + tmp2[];
	//immutable Chunk eta_tmp = S[] - 1.0 + tmp2[];

	Chunk nu = (-sign(coords.z)[]/sqrt(2.0))*sqrt(nu_tmp)[];
	//Chunk nu = (-sgn(coords.z)[]/sqrt(2.0))*sqrt(nu_tmp)[];
	nu = nu[].map!(a => a > 1.0 ? 1.0 : a).map!(a => a < -1.0 ? -1.0 : a).staticArray!Chunk[];
	immutable Chunk eta = (1.0/sqrt(2.0))*sqrt(eta_tmp)[];

	immutable Chunk neg_y = -coords.y[];
	immutable Chunk psi = atan2(neg_y, coords.x);

	return ElipticalCoords(nu, eta, psi);
}

private auto compute_velocities_nh(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients, T[] coeffiecients_sin, immutable ElipticalCoords coords) {

	Chunk V = 0;

	size_t sin_idx = 0;
	auto _idx = iterate_odds!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Pmn = infl.P_nh_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			if(m == 0) {
				immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
				V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			} else {
				immutable Chunk sin_mpsi = infl.sin_mpsi_buff[m];
				immutable Chunk beta = coeffiecients_sin[sin_idx];
				immutable Chunk Qmn_bar = infl.Qmn_bar[m][m+1][];
				sin_idx++;
				immutable Chunk a_cos = alpha[]*cos_mpsi[];
				immutable Chunk b_sin = beta[]*sin_mpsi[];

				V[] += Pmn[]*Qmn_bar[]*(a_cos[] + b_sin[]);
			}
		}
	)(infl.Mo, 0);

	iterate_evens!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			//immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			if(m == 0) {
				immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
				V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			} else {
				immutable Chunk sin_mpsi = infl.sin_mpsi_buff[m];
				immutable Chunk beta = coeffiecients_sin[sin_idx];
				immutable Chunk Qmn_bar = infl.Qmn_bar[m][m+1][];
				sin_idx++;
				immutable Chunk a_cos = alpha[]*cos_mpsi[];
				immutable Chunk b_sin = beta[]*sin_mpsi[];

				V[] += Pmn[]*Qmn_bar[]*(a_cos[] + b_sin[]);
				//V[] += Pmn[]*Qmn_bar[]*(alpha[]*cos_mpsi[] + beta[]*sin_mpsi[]);
			}
			//V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
		}
	)(infl.Me, _idx);

	return V;
}

private auto compute_velocities_md(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients, T[] coeffiecients_sin, immutable ElipticalCoords coords) {

	Chunk V = 0;

	size_t sin_idx = 0;
	auto _idx = iterate_odds!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			//immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][n][]*cos_mpsi[];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			//V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			if(m == 0) {
				immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
				V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			} else {
				immutable Chunk sin_mpsi = infl.sin_mpsi_buff[m];
				immutable Chunk beta = coeffiecients_sin[sin_idx];
				immutable Chunk Qmn_bar = infl.Qmn_bar[m][m+1][];
				sin_idx++;
				V[] += Pmn[]*Qmn_bar[]*(alpha[]*cos_mpsi[] + beta[]*sin_mpsi[]);
			}

		}
	)(infl.Mo, 0);

	iterate_evens!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			//immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][n][]*cos_mpsi[];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			//V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			if(m == 0) {
				immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
				V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
			} else {
				immutable Chunk sin_mpsi = infl.sin_mpsi_buff[m];
				immutable Chunk beta = coeffiecients_sin[sin_idx];
				immutable Chunk Qmn_bar = infl.Qmn_bar[m][m+1][];
				sin_idx++;
				V[] += Pmn[]*Qmn_bar[]*(alpha[]*cos_mpsi[] + beta[]*sin_mpsi[]);
			}
		}
	)(infl.Me, _idx);

	return V;
}
/+
@nogc private auto compute_velocities_nh(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients, T[] coeffiecients_sin, immutable ElipticalCoords coords) {

	Chunk V = 0;

	auto _idx = iterate_odds!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
			immutable Chunk Pmn = infl.P_nh_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
		}
	)(infl.Mo, 0);

	iterate_evens!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][m+1][]*cos_mpsi[];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
		}
	)(infl.Me, _idx);

	return V;
}

@nogc private auto compute_velocities_md(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients, T[] coeffiecients_sin, immutable ElipticalCoords coords) {

	Chunk V = 0;

	auto _idx = iterate_odds!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][n][]*cos_mpsi[];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];

		}
	)(infl.Mo, 0);

	iterate_evens!(
		(m, n, idx) {
			immutable Chunk cos_mpsi = infl.mpsi_buff[m];
			immutable Chunk Qmn_cos_mpsi = infl.Qmn_bar[m][n][]*cos_mpsi[];
			immutable Chunk Pmn = infl.P_md_mn_bar[m][n];

			immutable Chunk alpha = coefficients[idx];
			V[] += alpha[]*Pmn[]*Qmn_cos_mpsi[];
		}
	)(infl.Me, _idx);

	return V;
}
+/
private auto compute_velocities_bl(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients_md, T[] coefficients_nh, T[] coefficients_md_sin, T[] coefficients_nh_sin, immutable CartisianCoords ccoords, immutable ElipticalCoords coords) {

	// Recomended buffer zone in Huang's dissertation
	immutable double eps = 0.01;

	immutable Chunk h =  zip(coords.eta[], eps.repeat).map!(a => a[0] < a[1] ? 0 : a[0] - a[1]).staticArray!Chunk;
	immutable abs_y = abs(ccoords.y);
	
	immutable Chunk b = zip(ccoords.x[], abs_y[], coords.eta[], infl.sin_chi.repeat).map!((xye) {
		if(xye[0].isNaN || xye[1].isNaN || xye[2].isNaN || xye[3].isNaN) {
			return 0;
		}
		if(xye[0] <= 0.0 && xye[1] <= 1.0) {
			return 20.0*(1.0 - xye[1]*xye[1]*xye[3]/(1.0 + xye[2]*xye[2]));
		} else if(xye[0] <= 0.0 && xye[1] > 1.0) {
			return 20.0*(1.0 - xye[1]*xye[1]*xye[3]/((1.0 + xye[2]*xye[2]) + 0.615*(xye[1]*xye[1] - 1.0)));
		} else if(xye[0] > 0.0 && xye[1] <= 1.0) {
			return 20.0*(1.0 - (xye[0]*xye[0] + xye[1]*xye[1])*xye[3]/(1.0 + xye[2]*xye[2]));
		} else if(xye[0] > 0.0 && xye[1] > 1.0) {
			return 20.0*(1.0 - (xye[0]*xye[0] + xye[1]*xye[1])*xye[3]/((1.0 + xye[2]*xye[2]) + 0.615*(xye[1]*xye[1] - 1.0)));
		}
		assert(false);
	}).staticArray!Chunk;

	immutable Chunk blend_denom = 1.0 + b[]*h[];

	immutable Chunk blend_nh = 1.0/blend_denom[];
	immutable Chunk blend_md = b[]*h[]/blend_denom[];
	
	immutable bool compute_nh = blend_nh[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);
	immutable bool compute_md = blend_md[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);

	auto _idx = iterate_odds!(
			(m, n, idx) {
				if(compute_md) {
					immutable Chunk Pmn_md = associated_legendre_polynomial(m, n, coords.nu, infl.P_coefficients[idx]);
					infl.P_md_mn_bar[m][n][] = Pmn_md[];
				}

				if(compute_nh) {
					immutable Chunk Pmn_nh = associated_legendre_polynomial_nh(m, n, coords.nu, infl.P_coefficients_nh[idx]);
					infl.P_nh_mn_bar[m][n][] = Pmn_nh[];
				}
			}
	)(infl.Mo, 0);

	iterate_evens!(
		(m, n, idx) {
			immutable Chunk Pmn_md = associated_legendre_polynomial(m, n, coords.nu, infl.P_coefficients[idx]);
			infl.P_md_mn_bar[m][n][] = Pmn_md[];
		}
	)(infl.Me, _idx);

	associated_legendre_function(coords.eta, infl.Qmn_bar, infl.K_table);

	foreach(m; 0..max(infl.Mo, infl.Me) + 1) {
		if(m == 0) {
			infl.mpsi_buff[m] = 1.0;
			infl.sin_mpsi_buff[m] = 0.0;
		} else {
			immutable Chunk mpsi = m.to!double*coords.psi[];
			immutable Chunk[2] sin_cos = sincos(mpsi);
			infl.mpsi_buff[m] = sin_cos[1];
			infl.sin_mpsi_buff[m] = sin_cos[0];
		}
	}

	immutable Chunk V_md = compute_md ? compute_velocities_md(infl, coefficients_md, coefficients_md_sin, coords) : Z;
	immutable Chunk V_nh = compute_nh ? compute_velocities_nh(infl, coefficients_nh, coefficients_nh_sin, coords) : Z;
	
	immutable Chunk V = blend_nh[]*V_nh[] + blend_md[]*V_md[];

	return V;
}

@nogc private auto final_blend(double cos_chi, double sin_chi, immutable Chunk y, immutable Chunk sigma, immutable Chunk z, immutable Chunk x, immutable Chunk s_0) {

	immutable g = 1.84*sqrt(cos_chi) - 4.06*cos_chi + 11.84*cos_chi^^(1.5);
	immutable double sin_xi_2 = sin_chi*sin_chi;

	immutable Chunk f =
		/+zip(sigma[], sin_xi_2.repeat, g.repeat, y[], z[], x[], s_0[])
		.map!(
			(sxgy) {
				if((sxgy[0] < 0.0) || (sxgy[1] <= 1.0e-14) || ((sxgy[0] < 0.0) && (abs(sxgy[4]) <= 1.0e-14)) || (sxgy[5] > -sxgy[6]) /+|| (sxgy[4] > 0.0)+/) {
					return 0.0;
				} else {
					if(abs(sxgy[3]) <= 1.0) {
						return sxgy[1]/(sxgy[1] + sxgy[0]*sxgy[2]);
					} else {
						return sxgy[1]/(sxgy[1] + (sxgy[0] + 1.5*sqrt(sxgy[3]*sxgy[3] - 1.0))*sxgy[2]);
					}
				}
				assert(false);
			}
		)
		.staticArray!Chunk;+/
		zip(sigma[], sin_xi_2.repeat, g.repeat, y[], z[])
		.map!(
			(sxgy) {
				if((sxgy[0] < 0.0) || (sxgy[1] <= 1.0e-14) || ((sxgy[0] < 0.0) && (abs(sxgy[4]) <= 1.0e-14))) {
					return 0.0;
				} else {
					if(abs(sxgy[3]) <= 1.0) {
						return sxgy[1]/(sxgy[1] + sxgy[0]*sxgy[2]);
					} else {
						return sxgy[1]/(sxgy[1] + (sxgy[0] + 1.5*sqrt(sxgy[3]*sxgy[3] - 1.0))*sxgy[2]);
					}
				}
				assert(false);
			}
		)
		.staticArray!Chunk;
	return f;
}

private auto compute_velocities_final(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] a, T[] alpha, T[] delta, T[] lambda, T[] b, T[] beta, T[] delta_s, T[] lambda_s, immutable CartisianCoords ccoords, Chunk t, immutable Chunk x_0, immutable Chunk z_0, bool above_disk) {

	immutable auto coords = to_eliptical(ccoords);

	// if(infl.debug_coords) {
	// 	writeln("cart coords: ", ccoords);
	// 	writeln("elip coords: ", coords);
	// }

	immutable Chunk rho_axial = 1;

	immutable Chunk y2z2 = ccoords.y[]*ccoords.y[] + ccoords.z[]*ccoords.z[];
	//immutable Chunk y2z2 = ccoords.y[]*ccoords.y[] + z_0[]*z_0[];
	//immutable Chunk y2z2 = ccoords.y[]*ccoords.y[];
	immutable Chunk x2y2z2 = ccoords.x[]*ccoords.x[] + ccoords.y[]*ccoords.y[] + ccoords.z[]*ccoords.z[];
	immutable Chunk rho2 = rho_axial[]*rho_axial[];
	immutable Chunk s_0 = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;
	// Chunk s_0 = 0;
	
	// if(above_disk) {
	// 	//s_0[] = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;
	// 	s_0[] = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;//[] - z_0[]*infl.tan_chi;
	// 	//s_0[] = zip(y2z2[], rho2[], ccoords.y[]).map!(a => ((a[0] < a[1]) && (abs(a[2]) <= a[1])) ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[];
	// } else {
	// 	s_0[] = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[] - z_0[]*infl.tan_chi;
	// 	//s_0[] = zip(y2z2[], rho2[], ccoords.y[]).map!(a => ((a[0] < a[1]) && (abs(a[2]) <= a[1])) ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[] + z_0[]*infl.tan_chi;
	// }

	immutable Chunk neg_s_0 = -s_0[];
	immutable Chunk neg_y = -ccoords.y[];
	immutable Chunk sigma = -ccoords.x[] - s_0[];
	immutable Chunk sigma_p_s = sigma[] + s_0[];

	//immutable bool compute_ds = zip(x2y2z2[], ccoords.x[]).map!(a => (a[1] <= 0.0) && (a[0] >= rho2[0])).fold!((res, a) => res |= a)(false);
	
	//Chunk f = 0.0;
	//Chunk f = final_blend(infl.cos_chi, infl.sin_chi, ccoords.y, sigma, z_0, x_0, s_0);
	Chunk f = final_blend(infl.cos_chi, infl.sin_chi, ccoords.y, sigma, ccoords.z, ccoords.x, s_0);
	immutable bool compute_ds = f[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);
	
	immutable Chunk V_bl = compute_velocities_bl(infl, a, alpha, b, beta, ccoords, coords);
	Chunk V_ds = 0.0;
	if(compute_ds) {

		//f = final_blend(infl.cos_chi, infl.sin_chi, ccoords.y, sigma, ccoords.z);
		double devisor;
		devisor = infl.advance_ratio;

		Chunk t_delay = t[] - sigma[]*infl.sin_chi;///devisor;

		foreach(c_idx, ref t_d; t_delay) {
			if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
				while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
					if(abs(f[c_idx]) > 1.0e-14) {
						writeln("Dealing with the future: ", ccoords.x[c_idx], ", ", ccoords.y[c_idx], ", ", ccoords.z[c_idx], ", ", ccoords.x[c_idx]^^2.0+ccoords.y[c_idx]^^2.0 + ccoords.z[c_idx]^^2.0, ", f: ", f[c_idx]);
					}
					
					t_d -= 2.0*PI;
				}
			}
		}

		//t_delay[] = t_delay[]/devisor;

		infl.interpolate_to_time(t_delay, infl.time_delay_alpha_2, infl.time_delay_a_2, infl.time_delay_beta_2, infl.time_delay_b_2);

		immutable CartisianCoords[3] ccoords_ds = [
			CartisianCoords(neg_s_0, ccoords.y, ccoords.z),
			CartisianCoords(s_0, neg_y, ccoords.z),
			CartisianCoords(sigma_p_s, neg_y, ccoords.z)
		];

		immutable ElipticalCoords[3] coords_ds = [
			to_eliptical(ccoords_ds[0]),
			to_eliptical(ccoords_ds[1]),
			to_eliptical(ccoords_ds[2])
		];

		immutable Chunk v_ds_bl = compute_velocities_bl(infl,
			infl.time_delay_a_2[0..infl.total_states],
			infl.time_delay_alpha_2[0..infl.total_states],
			
			//infl.c_zero,
			//infl.c_zero,
			infl.time_delay_b_2[0..infl.total_sin_states],
			infl.time_delay_beta_2[0..infl.total_sin_states],
			
			ccoords_ds[0],
			coords_ds[0]
		);
		immutable Chunk v_ds_bl_a_1 = compute_velocities_bl(infl,
			infl.time_delay_a_2[infl.total_states..$],
			infl.time_delay_alpha_2[infl.total_states..$],
			
			//infl.c_zero,
			//infl.c_zero,
			infl.time_delay_b_2[infl.total_sin_states..$],
			infl.time_delay_beta_2[infl.total_sin_states..$],
			
			ccoords_ds[1],
			coords_ds[1]
		);
		immutable Chunk v_ds_bl_a_2 = compute_velocities_bl(infl,
			delta,
			lambda,
			delta_s,
			lambda_s,
			ccoords_ds[2],
			coords_ds[2]
		);

		V_ds[] = v_ds_bl[] + v_ds_bl_a_1[] - v_ds_bl_a_2[];
	}
	Chunk one_m_f = 1.0 - f[];

	immutable Chunk V_ds_f = V_ds[]*f[];
	immutable Chunk V = V_bl[]*one_m_f[] + V_ds_f[];

	return V;
}

private auto compute_velocities_final_adjoint(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] a, T[] alpha, T[] delta, T[] lambda, T[] b, T[] beta, T[] delta_s, T[] lambda_s, immutable CartisianCoords ccoords, Chunk t, immutable Chunk x_0, immutable Chunk z_0) {

	immutable auto coords = to_eliptical(ccoords);

	immutable Chunk rho_axial = 1;

	immutable Chunk y2z2 = ccoords.y[]*ccoords.y[] + ccoords.z[]*ccoords.z[];
	//immutable Chunk y2z2 = ccoords.y[]*ccoords.y[] + z_0[]*z_0[];
	//immutable Chunk y2z2 = ccoords.y[]*ccoords.y[];
	immutable Chunk x2y2z2 = ccoords.x[]*ccoords.x[] + ccoords.y[]*ccoords.y[] + ccoords.z[]*ccoords.z[];
	immutable Chunk rho2 = rho_axial[]*rho_axial[];
	immutable Chunk s_0 = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;
	//immutable Chunk s_0 = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[] - z_0[]*infl.tan_chi;
	//immutable Chunk s_0 = zip(y2z2[], rho2[], ccoords.y[]).map!(a => ((a[0] < a[1]) && (abs(a[2]) <= a[1])) ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[] + z_0[]*infl.tan_chi;

	// Chunk s_0 = 0;
	
	// if(above_disk) {
	// 	//s_0[] = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;
	// 	s_0[] = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;//[] - z_0[]*infl.tan_chi;
	// 	//s_0[] = zip(y2z2[], rho2[], ccoords.y[]).map!(a => ((a[0] < a[1]) && (abs(a[2]) <= a[1])) ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[];
	// } else {
	// 	s_0[] = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[] - z_0[]*infl.tan_chi;
	// 	//s_0[] = zip(y2z2[], rho2[], ccoords.y[]).map!(a => ((a[0] < a[1]) && (abs(a[2]) <= a[1])) ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[] + z_0[]*infl.tan_chi;
	// }


	immutable Chunk neg_s_0 = -s_0[];
	immutable Chunk neg_y = -ccoords.y[];
	immutable Chunk sigma = -ccoords.x[] - s_0[];
	immutable Chunk sigma_p_s = sigma[] + s_0[];

	//Chunk f = 0.0;
	//Chunk f = final_blend(infl.cos_chi, infl.sin_chi, ccoords.y, sigma, z_0, x_0, s_0);
	Chunk f = final_blend(infl.cos_chi, infl.sin_chi, ccoords.y, sigma, ccoords.z, ccoords.x, s_0);

	//immutable bool compute_ds = zip(x2y2z2[], ccoords.x[]).map!(a => (a[1] <= 0.0) && (a[0] >= rho2[0])).fold!((res, a) => res |= a)(false);
	immutable bool compute_ds = f[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);

	immutable Chunk V_bl = compute_velocities_bl(infl, delta, lambda, delta_s, lambda_s, ccoords, coords);
	Chunk V_ds = 0.0;

	if(compute_ds) {

		//f = final_blend(infl.cos_chi, infl.sin_chi, ccoords.y, sigma, ccoords.z);

		double devisor;
		devisor = infl.advance_ratio;

		Chunk t_delay = t[] + sigma[]*infl.sin_chi;///devisor;

		foreach(c_idx, ref t_d; t_delay) {
			if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
				while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
					t_d -= 2.0*PI;
				}
			}
		}

		//t_delay[] = t_delay[]/devisor;

		infl.interpolate_to_time(t_delay, infl.time_delay_alpha_2, infl.time_delay_a_2, infl.time_delay_beta_2, infl.time_delay_b_2);

		immutable CartisianCoords[3] ccoords_ds = [
			CartisianCoords(neg_s_0, ccoords.y, ccoords.z),
			CartisianCoords(s_0, neg_y, ccoords.z),
			CartisianCoords(sigma_p_s, neg_y, ccoords.z)
		];

		immutable ElipticalCoords[3] coords_ds = [
			to_eliptical(ccoords_ds[0]),
			to_eliptical(ccoords_ds[1]),
			to_eliptical(ccoords_ds[2])
		];

		immutable Chunk v_ds_bl_a = compute_velocities_bl(infl,
			infl.time_delay_a_2[infl.total_states..$],
			infl.time_delay_alpha_2[infl.total_states..$],
			
			//infl.c_zero,
			//infl.c_zero,
			infl.time_delay_b_2[infl.total_sin_states..$],
			infl.time_delay_beta_2[infl.total_sin_states..$],

			ccoords_ds[0],
			coords_ds[0]
		);

		immutable Chunk v_ds_bl_1 = compute_velocities_bl(infl,
			infl.time_delay_a_2[0..infl.total_states],
			infl.time_delay_alpha_2[0..infl.total_states],
			
			//infl.c_zero,
			//infl.c_zero,
			infl.time_delay_b_2[0..infl.total_sin_states],
			infl.time_delay_beta_2[0..infl.total_sin_states],

			ccoords_ds[1],
			coords_ds[1]
		);
		immutable Chunk v_ds_bl_2 = compute_velocities_bl(infl,
			a,
			alpha,
			b,
			beta,
			ccoords_ds[2],
			coords_ds[2]
		);

		V_ds[] = v_ds_bl_a[] + v_ds_bl_1[] - v_ds_bl_2[];
	}
	
	Chunk one_m_f = 1.0 - f[];

	immutable Chunk V_ds_f = V_ds[]*f[];
	immutable Chunk V = V_bl[]*one_m_f[] + V_ds_f[];

	return V;
}

private auto under_disk_solution(ArrayContainer AC, T, C)(HuangPetersInflowT!AC infl, T[] lambda, T[] delta, T[] lambda_s, T[] delta_s, auto ref C x, auto ref C y, auto ref C z, Chunk t) {
	immutable tan_chi = tan(infl.chi);
	immutable Chunk x_offset = z[]*tan_chi;
	immutable Chunk x0_1 = x[] + x_offset[];
	immutable Chunk y0_1 = y[];
	immutable Chunk x0_2 = -x0_1[];
	immutable Chunk y0_2 = -y[];
	immutable Chunk zero = 0;

	Chunk s = z[]*z[] + x_offset[]*x_offset[];
	s = sqrt(s);
	//immutable Chunk t = infl.times[infl.get_circular_index(infl.curr_state)];

	double devisor;
	if(abs(infl.chi) < 1.0e-12) {
		devisor = ((infl.v_0 + infl.axial_advance_ratio*infl.cos_chi));
	} else {
		devisor = ((infl.v_0 + infl.advance_ratio/infl.sin_chi + infl.axial_advance_ratio*infl.cos_chi));
	}
	Chunk t_minus = s[];
	Chunk t_delay = t[] - t_minus[];///devisor;

	foreach(c_idx, ref t_d; t_delay) {
		if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {

			while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
				t_d -= 2.0*PI;
			}
		}
	}

	//t_delay[] = t_delay[]/devisor;
	infl.interpolate_to_time(t_delay, infl.time_delay_alpha, infl.time_delay_a, infl.time_delay_beta, infl.time_delay_b);

	CartisianCoords ccoords1 = CartisianCoords(x0_1, y0_1, zero);
	CartisianCoords ccoords2 = CartisianCoords(x0_2, y0_2, zero);
	
	immutable Chunk neg_z = -z[];
	immutable Chunk x3 = -x[];
	CartisianCoords ccoords3 = CartisianCoords(x3, y0_2, neg_z);

	ElipticalCoords coords1 = to_eliptical(ccoords1);
	ElipticalCoords coords2 = to_eliptical(ccoords2);
	ElipticalCoords coords3 = to_eliptical(ccoords3);

	//private auto compute_velocities_bl(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients_md, T[] coefficients_nh, T[] coefficients_md_sin, T[] coefficients_nh_sin, immutable CartisianCoords ccoords, immutable ElipticalCoords coords) {
	immutable Chunk v_f = compute_velocities_bl(infl,
		infl.time_delay_a[0..infl.total_states],
		infl.time_delay_alpha[0..infl.total_states],
		//infl.time_delay_a[infl.total_states..$],
		//infl.time_delay_alpha[infl.total_states..$],
		infl.time_delay_b[0..infl.total_sin_states],
		infl.time_delay_beta[0..infl.total_sin_states],
		
		//infl.c_zero,
		//infl.c_zero,
		//infl.time_delay_b[infl.total_sin_states..$],
		//infl.time_delay_beta[infl.total_sin_states..$],

		ccoords1,
		coords1
		//t_delay,
		//x,
		//z,
		//false
	);

	immutable Chunk v_f_a_1 = compute_velocities_bl(infl,
		//infl.time_delay_a[0..infl.total_states],
		//infl.time_delay_alpha[0..infl.total_states],
		infl.time_delay_a[infl.total_states..$],
		infl.time_delay_alpha[infl.total_states..$],
		
		//infl.c_zero,
		//infl.c_zero,
		//infl.time_delay_b[0..infl.total_sin_states],
		//infl.time_delay_beta[0..infl.total_sin_states],

		infl.time_delay_b[infl.total_sin_states..$],
		infl.time_delay_beta[infl.total_sin_states..$],
		ccoords2,
		coords2
		//coords2,
		//t_delay,
		//x,
		//z
	);

	

	immutable Chunk v_f_a_2 = compute_velocities_bl(infl,
		//infl.a,
		//infl.alpha,
		delta,
		lambda,
		
		//infl.zero,
		//infl.zero,
		//infl.b,
		//infl.beta,

		delta_s,
		lambda_s,

		ccoords3,
		coords3

		//coords3,
		//t,
		//x,
		//z
	);

	immutable Chunk V_below = v_f[] + v_f_a_1[] - v_f_a_2[];
	return V_below;
}

private auto under_disk_solution_adjoint(ArrayContainer AC, T, C)(HuangPetersInflowT!AC infl, T[] alpha, T[] a, T[] beta, T[] b, auto ref C x, auto ref C y, auto ref C z, Chunk t) {
	immutable tan_chi = tan(infl.chi);
	immutable Chunk x_offset = z[]*tan_chi;
	immutable Chunk x0_1 = x[] + x_offset[];
	immutable Chunk y0_1 = y[];
	immutable Chunk x0_2 = -x0_1[];
	immutable Chunk y0_2 = -y[];
	immutable Chunk zero = 0;

	Chunk s = z[]*z[] + x_offset[]*x_offset[];
	s = sqrt(s);

	/+double devisor;
	if(abs(infl.chi) < 1.0e-12) {
		devisor = ((infl.v_0 + infl.axial_advance_ratio*infl.cos_chi));
	} else {
		devisor = ((infl.v_0 + infl.advance_ratio/infl.sin_chi + infl.axial_advance_ratio*infl.cos_chi));
	}+/
	Chunk t_minus = s[];
	Chunk t_delay = t[] + t_minus[];///devisor;

	foreach(c_idx, ref t_d; t_delay) {
		if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {

			while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
				//t_d -= 2.0*PI;
				t_d -= 1.0;
			}
		}
	}

	//t_delay[] = t_delay[]/devisor;
	infl.interpolate_to_time(t_delay, infl.time_delay_alpha, infl.time_delay_a, infl.time_delay_beta, infl.time_delay_b);

	CartisianCoords ccoords1 = CartisianCoords(x0_1, y0_1, zero);
	CartisianCoords ccoords2 = CartisianCoords(x0_2, y0_2, zero);
	
	immutable Chunk neg_z = -z[];
	immutable Chunk x3 = -x[];
	CartisianCoords ccoords3 = CartisianCoords(x3, y0_2, neg_z);

	ElipticalCoords coords1 = to_eliptical(ccoords1);
	ElipticalCoords coords2 = to_eliptical(ccoords2);
	ElipticalCoords coords3 = to_eliptical(ccoords3);

	//private auto compute_velocities_bl(ArrayContainer AC, T)(HuangPetersInflowT!AC infl, T[] coefficients_md, T[] coefficients_nh, T[] coefficients_md_sin, T[] coefficients_nh_sin, immutable CartisianCoords ccoords, immutable ElipticalCoords coords) {

	// a + r - r
	immutable Chunk v_f_a = compute_velocities_bl(infl,
		infl.time_delay_a[infl.total_states..$],
		infl.time_delay_alpha[infl.total_states..$],
		infl.time_delay_b[infl.total_sin_states..$],
		infl.time_delay_beta[infl.total_sin_states..$],
		ccoords1,
		coords1
	);

	immutable Chunk v_f_1 = compute_velocities_bl(infl,
		infl.time_delay_a[0..infl.total_states],
		infl.time_delay_alpha[0..infl.total_states],
		infl.time_delay_b[0..infl.total_sin_states],
		infl.time_delay_beta[0..infl.total_sin_states],
		ccoords2,
		coords2
	);

	immutable Chunk v_f_2 = compute_velocities_bl(infl,
		a,
		alpha,
		b,
		beta,
		ccoords3,
		coords3
	);

	immutable Chunk V_below = v_f_a[] + v_f_1[] - v_f_2[];
	return V_below;
}

private auto inflow_at_impl(ArrayContainer AC, C)(HuangPetersInflowT!AC infl, auto ref C x, auto ref C y, auto ref C z) {
	
	static import opencopter.math;

	bool all_nan = x[].map!(a => a.isNaN).fold!((res, a) => res && a)(true);
	if(all_nan) {
		Chunk V = 0;
		return V;
	}

	immutable all_below_disk = z[].map!(a => a > 0).fold!((res, a) => a && res)(true);
	immutable all_above_or_on_disk = z[].map!(a => (a <= 0)).fold!((res, a) => a && res)(true);
	immutable all_somewhere = !all_above_or_on_disk && !all_below_disk;// && !all_in_upper_hemi;

	import core.stdc.stdio : printf;

	CartisianCoords coords = CartisianCoords(x, y, z);

	Chunk V_above = 0;
	Chunk V_below = 0;

	if(all_above_or_on_disk || all_somewhere) {
		Chunk t = infl.times[infl.get_circular_index(infl.curr_state)];
		// if(infl.debug_coords) {
		// 	writeln("infl.b: ", infl.b);
		// 	writeln("infl.beta: ", infl.beta);
		// }
		V_above = compute_velocities_final(infl, infl.a, infl.alpha, infl.delta, infl.lambda, infl.b, infl.beta, infl.delta_s, infl.lambda_s, coords, t, x, z, true);
		if(all_above_or_on_disk && !all_somewhere) {
			return V_above;
		}
	}
	
	if(all_below_disk || all_somewhere) {
		/+
		immutable Chunk rho_axial = 1;
		//immutable Chunk y2z2 = y[]*y[];
		immutable Chunk y2z2 = y[]*y[];// + z[]*z[];
		//immutable Chunk y2z2 = y[]*y[] + z[]*z[];
		//immutable Chunk x2y2z2 = ccoords.x[]*ccoords.x[] + ccoords.y[]*ccoords.y[] + ccoords.z[]*ccoords.z[];
		immutable Chunk rho2 = rho_axial[]*rho_axial[];
		immutable Chunk s_0 = zip(y2z2[], rho2[]).map!(a => a[0] < a[1] ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk;
		//immutable Chunk s_0 = zip(y2z2[], rho2[], y[]).map!(a => ((a[0] < a[1]) && (abs(a[2]) <= a[1])) ? sqrt(a[1] - a[0]) : 0).staticArray!Chunk[] + z[]*infl.tan_chi;
	
		immutable all_downstream = zip(x[], s_0[]).map!(a => a[0] <= -a[1]).fold!((res, a) => a && res)(true);
		immutable all_elsewhere = zip(x[], s_0[]).map!(a => a[0] > -a[1]).fold!((res, a) => a && res)(true);

		immutable all_somewhere_under = !all_downstream && !all_elsewhere;

		//Chunk V_ds_below = 0;
		//Chunk V_below = 0;

		immutable Chunk t = infl.times[infl.get_circular_index(infl.curr_state)];

		V_below = under_disk_solution(infl, infl.lambda, infl.delta, infl.lambda_s, infl.delta_s, x, y, z, t);



		Chunk V_ds = 0;
		immutable Chunk neg_s_0 = -s_0[];
		immutable Chunk neg_y = -y[];
		immutable Chunk sigma = -x[] - s_0[];
		immutable Chunk sigma_p_s = sigma[] + s_0[];

		//immutable bool compute_ds = zip(x2y2z2[], ccoords.x[]).map!(a => (a[1] <= 0.0) && (a[0] >= rho2[0])).fold!((res, a) => res |= a)(false);
		
		//Chunk f = 0.0;
		Chunk f = final_blend(infl.cos_chi, infl.sin_chi, y, sigma, z, x, s_0);
		immutable bool compute_ds = f[].map!(a => !a.isClose(0.0)).fold!((res, a) => res |= a)(false);

		if((all_downstream || all_somewhere_under) && (abs(infl.chi) > 1.0e-14)) {

			Chunk t_delay = t[] - sigma[];//*infl.sin_chi;///devisor;
			//Chunk t_delay_2 = t[] + sigma[]*infl.sin_chi;///devisor;

			foreach(c_idx, ref t_d; t_delay) {
				if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
					while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
						if(abs(f[c_idx]) > 1.0e-14) {
							//writeln("Dealing with the future: ", ccoords.x[c_idx], ", ", ccoords.y[c_idx], ", ", ccoords.z[c_idx], ", ", ccoords.x[c_idx]^^2.0+ccoords.y[c_idx]^^2.0 + ccoords.z[c_idx]^^2.0, ", f: ", f[c_idx]);
						}
						
						t_d -= 2.0*PI;
					}
				}
			}

			/+foreach(c_idx, ref t_d; t_delay_2) {
				if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
					while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
						if(abs(f[c_idx]) > 1.0e-14) {
							//writeln("Dealing with the future: ", ccoords.x[c_idx], ", ", ccoords.y[c_idx], ", ", ccoords.z[c_idx], ", ", ccoords.x[c_idx]^^2.0+ccoords.y[c_idx]^^2.0 + ccoords.z[c_idx]^^2.0, ", f: ", f[c_idx]);
						}
						
						t_d -= 2.0*PI;
					}
				}
			}+/

			//t_delay[] = t_delay[]/devisor;

			infl.interpolate_to_time(t_delay, infl.time_delay_alpha_2, infl.time_delay_a_2, infl.time_delay_beta_2, infl.time_delay_b_2);

			immutable CartisianCoords[3] ccoords_ds = [
				CartisianCoords(neg_s_0, y, z),
				CartisianCoords(s_0, neg_y, z),
				CartisianCoords(sigma_p_s, neg_y, z)
			];

			/+immutable ElipticalCoords[3] coords_ds = [
				to_eliptical(ccoords_ds[0]),
				to_eliptical(ccoords_ds[1]),
				to_eliptical(ccoords_ds[2])
			];+/
			//private auto under_disk_solution(ArrayContainer AC, C)(HuangPetersInflowT!AC infl, C lambda, C delta, C lambda_s, C delta_s, auto ref C x, auto ref C y, auto ref C z, Chunk t)
			immutable Chunk v_ds_bl = under_disk_solution(infl,
				infl.time_delay_a_2[infl.total_states..$],
				infl.time_delay_alpha_2[infl.total_states..$],
				infl.time_delay_b_2[infl.total_sin_states..$],
				infl.time_delay_beta_2[infl.total_sin_states..$],
				ccoords_ds[0].x,
				ccoords_ds[0].y,
				ccoords_ds[0].z,
				t_delay
			);

			immutable Chunk v_ds_bl_a_1 = under_disk_solution_adjoint(infl,
				infl.time_delay_a_2[infl.total_states..$],
				infl.time_delay_alpha_2[infl.total_states..$],
				infl.time_delay_b_2[infl.total_sin_states..$],
				infl.time_delay_beta_2[infl.total_sin_states..$],
				ccoords_ds[1].x,
				ccoords_ds[1].y,
				ccoords_ds[1].z,
				t_delay
			);
			immutable Chunk v_ds_bl_a_2 = under_disk_solution_adjoint(infl,
				infl.a,
				infl.alpha,
				infl.b,
				infl.beta,
				ccoords_ds[2].x,
				ccoords_ds[2].y,
				ccoords_ds[2].z,
				t
			);

			V_ds[] = v_ds_bl[] + v_ds_bl_a_1[] - v_ds_bl_a_2[];
		}

		Chunk one_m_f = 1.0 - f[];

		immutable Chunk V_ds_f = V_ds[]*f[];
		V_below = V_below[]*one_m_f[] + V_ds_f[];
		+/
		
		immutable tan_chi = tan(infl.chi);
		immutable Chunk x_offset = z[]*tan_chi;
		immutable Chunk x0_1 = x[] + x_offset[];
		immutable Chunk y0_1 = y[];
		immutable Chunk x0_2 = -x0_1[];
		immutable Chunk y0_2 = -y[];
		immutable Chunk zero = 0;

		Chunk s = z[]*z[] + x_offset[]*x_offset[];
		s = sqrt(s);
		immutable Chunk t = infl.times[infl.get_circular_index(infl.curr_state)];

		double devisor;
		if(abs(infl.chi) < 1.0e-12) {
			devisor = ((infl.v_0 + infl.axial_advance_ratio*infl.cos_chi));
		} else {
			devisor = ((infl.v_0 + infl.advance_ratio/infl.sin_chi + infl.axial_advance_ratio*infl.cos_chi));
		}
		Chunk t_minus = s[];
		Chunk t_delay = t[] - t_minus[];///devisor;

		foreach(c_idx, ref t_d; t_delay) {
			if(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {

				while(t_d > infl.times[infl.get_circular_index(infl.curr_state)]) {
					t_d -= 2.0*PI;
				}
			}
		}

		//t_delay[] = t_delay[]/devisor;
		infl.interpolate_to_time(t_delay, infl.time_delay_alpha, infl.time_delay_a, infl.time_delay_beta, infl.time_delay_b);

		CartisianCoords coords1 = CartisianCoords(x0_1, y0_1, zero);
		CartisianCoords coords2 = CartisianCoords(x0_2, y0_2, zero);

		immutable Chunk v_f = compute_velocities_final(infl,
			infl.time_delay_a[0..infl.total_states],
			infl.time_delay_alpha[0..infl.total_states],
			infl.time_delay_a[infl.total_states..$],
			infl.time_delay_alpha[infl.total_states..$],
			infl.time_delay_b[0..infl.total_sin_states],
			infl.time_delay_beta[0..infl.total_sin_states],
			
			//infl.c_zero,
			//infl.c_zero,
			infl.time_delay_b[infl.total_sin_states..$],
			infl.time_delay_beta[infl.total_sin_states..$],

			coords1,
			t_delay,
			x,
			z,
			false
		);

		immutable Chunk v_f_a_1 = compute_velocities_final_adjoint(infl,
			infl.time_delay_a[0..infl.total_states],
			infl.time_delay_alpha[0..infl.total_states],
			infl.time_delay_a[infl.total_states..$],
			infl.time_delay_alpha[infl.total_states..$],
			
			//infl.c_zero,
			//infl.c_zero,
			infl.time_delay_b[0..infl.total_sin_states],
			infl.time_delay_beta[0..infl.total_sin_states],

			infl.time_delay_b[infl.total_sin_states..$],
			infl.time_delay_beta[infl.total_sin_states..$],
			coords2,
			t_delay,
			x,
			z
		);

		immutable Chunk neg_z = -z[];
		
		immutable Chunk x3 = -x[];

		CartisianCoords coords3 = CartisianCoords(x3, y0_2, neg_z);

		immutable Chunk v_f_a_2 = compute_velocities_final_adjoint(infl,
			infl.a,
			infl.alpha,
			infl.delta,
			infl.lambda,
			
			//infl.zero,
			//infl.zero,
			infl.b,
			infl.beta,

			infl.delta_s,
			infl.lambda_s,
			coords3,
			t,
			x,
			z
		);

		V_below = v_f[] + v_f_a_1[] - v_f_a_2[];
		
		if(all_below_disk && !all_somewhere) {
			return V_below;
		}
	}

	// Merge results
	if(all_somewhere) {

		Chunk V = 0;
		foreach(c_idx; 0..chunk_size) {
			if(z[c_idx] > 0) {
				V[c_idx] = V_below[c_idx];
			} else {
				V[c_idx] = V_above[c_idx];
			}
		}
		return V;
	}

	assert(false);
}

private immutable odd_states = [
	[[1]],  // 0th
	[[1],[2]],  // 1st
	[[1, 3], [2], [3]],  // 2nd
	[[1, 3], [2, 4], [3], [4]],  // 3rd
	[[1, 3, 5], [2, 4], [3, 5], [4], [5]],  // 4th
	[[1, 3, 5], [2, 4, 6], [3, 5], [4, 6], [5], [6]],  // 5th
	[[1, 3, 5, 7], [2, 4, 6], [3, 5, 7], [4, 6], [5, 7], [6], [7]],  // 6th
	[[1, 3, 5, 7], [2, 4, 6, 8], [3, 5, 7], [4, 6, 8], [5, 7], [6, 8], [7], [8]],  // 7th
	[[1, 3, 5, 7, 9], [2, 4, 6, 8], [3, 5, 7, 9], [4, 6, 8], [5, 7, 9], [6, 8], [7, 9], [8], [9]],  // 8th
	[[1, 3, 5, 7, 9], [2, 4, 6, 8, 10], [3, 5, 7, 9], [4, 6, 8, 10], [5, 7, 9], [6, 8, 10], [7, 9], [8, 10], [9], [10]],  // 9th
	[[1, 3, 5, 7, 9, 11], [2, 4, 6, 8, 10], [3, 5, 7, 9, 11], [4, 6, 8, 10], [5, 7, 9, 11], [6, 8, 10], [7, 9, 11], [8, 10], [9, 11], [10], [11]],  // 10th
	[[1, 3, 5, 7, 9, 11], [2, 4, 6, 8, 10, 12], [3, 5, 7, 9, 11], [4, 6, 8, 10, 12], [5, 7, 9, 11], [6, 8, 10, 12], [7, 9, 11], [8, 10, 12], [9, 11], [10, 12], [11], [12]],  // 11th
	[[1, 3, 5, 7, 9, 11, 13], [2, 4, 6, 8, 10, 12], [3, 5, 7, 9, 11, 13], [4, 6, 8, 10, 12], [5, 7, 9, 11, 13], [6, 8, 10, 12], [7, 9, 11, 13], [8, 10, 12], [9, 11, 13], [10, 12], [11, 13], [12], [13]]  // 12th
];

private immutable even_states = [
	[ // 1st
		[2]
	],
	[ // 2nd
		[2], [3]
	],
	[ // 3rd
		[2, 4], [3], [4]
	],
	[ // 4th
		[2, 4], [3, 5], [4], [5]
	],
	[ // 5th
		[2, 4, 6], [3, 5], [4, 6], [5], [6]
	],
	[ // 6th
		[2, 4, 6], [3, 5, 7], [4, 6], [5, 7], [6], [7]
	],
	[ // 7th
		[2, 4, 6, 8], [3, 5, 7], [4, 6, 8], [5, 7], [6, 8], [7], [8]
	],
	[ // 8th
		[2, 4, 6, 8], [3, 5, 7, 9], [4, 6, 8], [5, 7, 9], [6, 8], [7, 9], [8], [9]
	],
	[ // 9th
		[2, 4, 6, 8, 10], [3, 5, 7, 9], [4, 6, 8, 10], [5, 7, 9], [6, 8, 10], [7, 9], [8, 10], [9], [10]
	],
	[ // 10th
		[2, 4, 6, 8, 10], [3, 5, 7, 9, 11], [4, 6, 8, 10], [5, 7, 9, 11], [6, 8, 10], [7, 9, 11], [8, 10], [9, 11], [10], [11]
	],
	[ // 11th
		[2, 4, 6, 8, 10, 12], [3, 5, 7, 9, 11], [4, 6, 8, 10, 12], [5, 7, 9, 11], [6, 8, 10, 12], [7, 9, 11], [8, 10, 12], [9, 11], [10, 12], [11], [12]
	],
	[ // 12th
		[2, 4, 6, 8, 10, 12], [3, 5, 7, 9, 11, 13], [4, 6, 8, 10, 12], [5, 7, 9, 11, 13], [6, 8, 10, 12], [7, 9, 11, 13], [8, 10, 12], [9, 11, 13], [10, 12], [11, 13], [12], [13]
	],
	[ // 13th
		[2, 4, 6, 8, 10, 12, 14], [3, 5, 7, 9, 11, 13], [4, 6, 8, 10, 12, 14], [5, 7, 9, 11, 13], [6, 8, 10, 12, 14], [7, 9, 11, 13], [8, 10, 12, 14], [9, 11, 13], [10, 12, 14], [11, 13], [12, 14], [13], [14]
	]
];

private size_t iterate_odds_sin(alias f, Args...)(long Mo, size_t idx, auto ref Args args) {
	foreach(long _m; 1..Mo + 1) {
		foreach(long _n; odd_states[Mo][_m]) {
			f(_m.to!long, _n.to!long, idx, args);
			idx++;
		}
	}

	return idx;
}

private size_t iterate_evens_sin(alias f, Args...)(long Me, size_t idx, auto ref Args args) {
	foreach(long _m; 1..Me + 1) {
		foreach(long _n; even_states[Me][_m]) {
			f(_m.to!long, _n.to!long, idx, args);
			idx++;
		}
	}

	return idx;
}

size_t iterate_even_odd_sin(alias f, Args...)(long Mo, long Me, size_t start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_evens_sin!((r, j, row_idx) {
		iterate_odds_sin!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Me, start_row_idx);
}

size_t iterate_odd_even_sin(alias f, Args...)(long Mo, long Me, size_t start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_odds_sin!((r, j, row_idx) {
		iterate_evens_sin!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Me, start_row_idx);
}

size_t iterate_odd_odd_sin(alias f, Args...)(long Mo, long start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_odds_sin!((r, j, row_idx) {
		iterate_odds_sin!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Mo, start_row_idx);
}

size_t iterate_even_even_sin(alias f, Args...)(long Me, long start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_evens_sin!((r, j, row_idx) {
		iterate_evens_sin!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Me, start_col_idx);
	})(Me, start_row_idx);
}

void iterate_whole_matrix_sin(alias f_oo, alias f_oe, alias f_eo, alias f_ee, Args...)(long Mo, long Me, auto ref Args args) {

	auto _row_idx = iterate_odds_sin!((r, j, row_idx) {
		auto _col_idx = iterate_odds_sin!(
			(m, n, col_idx) => f_oo(r, j, m, n, row_idx, col_idx, args)
		)(Mo, 0);

		iterate_evens_sin!(
			(m, n, col_idx) => f_oe(r, j, m, n, row_idx, col_idx, args)
		)(Me, _col_idx);
	})(Mo, 0);

	iterate_evens_sin!((r, j, row_idx) {
		auto _col_idx = iterate_odds_sin!(
			(m, n, col_idx) => f_eo(r, j, m, n, row_idx, col_idx, args)
		)(Mo, 0);

		iterate_evens_sin!(
			(m, n, col_idx) => f_ee(r, j, m, n, row_idx, col_idx, args)
		)(Me, _col_idx);
	})(Me, _row_idx);
}

private size_t iterate_odds(alias f, Args...)(long Mo, size_t idx, auto ref Args args) {
	foreach(long _m; 0..Mo + 1) {
		foreach(long _n; odd_states[Mo][_m]) {
			f(_m.to!long, _n.to!long, idx, args);
			idx++;
		}
	}

	return idx;
}

private size_t iterate_evens(alias f, Args...)(long Me, size_t idx, auto ref Args args) {
	foreach(long _m; 0..Me + 1) {
		foreach(long _n; even_states[Me + 0][_m]) {
			f(_m.to!long, _n.to!long, idx, args);
			idx++;
		}
	}

	return idx;
}

size_t iterate_even_odd(alias f, Args...)(long Mo, long Me, size_t start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_evens!((r, j, row_idx) {
		iterate_odds!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Me, start_row_idx);
}

size_t iterate_odd_even(alias f, Args...)(long Mo, long Me, size_t start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_odds!((r, j, row_idx) {
		iterate_evens!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Me, start_row_idx);
}

size_t iterate_odd_odd(alias f, Args...)(long Mo, long start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_odds!((r, j, row_idx) {
		iterate_odds!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Mo, start_col_idx);
	})(Mo, start_row_idx);
}

size_t iterate_even_even(alias f, Args...)(long Me, long start_row_idx, size_t start_col_idx, auto ref Args args) {
	return iterate_evens!((r, j, row_idx) {
		iterate_evens!(
			(m, n, col_idx) => f(r, j, m, n, row_idx, col_idx, args)
		)(Me, start_col_idx);
	})(Me, start_row_idx);
}

void iterate_whole_matrix(alias f_oo, alias f_oe, alias f_eo, alias f_ee, Args...)(long Mo, long Me, auto ref Args args) {

	auto _row_idx = iterate_odds!((r, j, row_idx) {
		auto _col_idx = iterate_odds!(
			(m, n, col_idx) => f_oo(r, j, m, n, row_idx, col_idx, args)
		)(Mo, 0);

		iterate_evens!(
			(m, n, col_idx) => f_oe(r, j, m, n, row_idx, col_idx, args)
		)(Me, _col_idx);
	})(Mo, 0);

	iterate_evens!((r, j, row_idx) {
		auto _col_idx = iterate_odds!(
			(m, n, col_idx) => f_eo(r, j, m, n, row_idx, col_idx, args)
		)(Mo, 0);

		iterate_evens!(
			(m, n, col_idx) => f_ee(r, j, m, n, row_idx, col_idx, args)
		)(Me, _col_idx);
	})(Me, _row_idx);
}

void build_vlm_matrix(ArrayContainer AC)(HuangPetersInflowT!AC infl, double advance_ratio, double axial_advance_ratio) {

	double[] alpha = infl.state_history[infl.get_circular_index(infl.curr_state)][0..infl.total_states];
	
	if(infl.average_inflow != 0) {
		infl.chi = atan2(advance_ratio, (infl.average_inflow + axial_advance_ratio));
		//infl.chi = atan2(advance_ratio, (infl.average_inflow));
	}

	infl.sin_chi = sin(infl.chi);
	infl.cos_chi = cos(infl.chi);
	infl.tan_chi = tan(infl.chi);
	//infl.cot_chi = cot(infl.chi);

	immutable X = tan(-infl.chi/2.0);

	//writeln("X: ", X);

	iterate_whole_matrix!(
		// odd-odd
		(r, j, m, n, row_idx, col_idx) {
			immutable l = min(r, m);
			immutable exp1 = abs(r - m);
			immutable exp2 = abs(r + m);

			immutable l_s = min(r + 1, m + 1);
			immutable exp1_s = abs((r + 1) - (m + 1));
			immutable exp2_s = abs((r + 1) + (m + 1));
			if(r == 0) {
				infl.L_c[row_idx][col_idx] = X^^(m.to!double)*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = (X^^exp1_s - ((-1.0)^^l_s)*X^^exp2_s)*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = 0;
			} else {
				//immutable double l = min(r, m);
				//immutable double exp1 = abs(r - m);
				//immutable double exp2 = abs(r + m);

				infl.L_c[row_idx][col_idx] = (X^^exp1 + ((-1.0)^^l)*X^^exp2)*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = (X^^exp1_s - ((-1.0)^^l_s)*X^^exp2_s)*infl.Gamma[row_idx][col_idx];
			}
		},
		// odd-even
		(r, j, m, n, row_idx, col_idx) {
			immutable l = min(r, m);
			immutable exp1 = abs(r - m);
			immutable exp2 = abs(r + m);

			immutable l_s = min(r + 1, m + 1);
			immutable exp1_s = abs((r + 1) - (m + 1));
			immutable exp2_s = abs((r + 1) + (m + 1));
			if(r == 0) {
				infl.L_c[row_idx][col_idx] = X^^(m.to!double)*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = 0;
				//infl.L_s[row_idx][col_idx] = (X^^exp1_s - ((-1.0)^^l_s)*X^^exp2_s)*infl.Gamma[row_idx][col_idx];
			} else {
				//immutable double l = min(r, m);
				//immutable double exp1 = abs(r - m);
				//immutable double exp2 = abs(r + m);

				infl.L_c[row_idx][col_idx] = (X^^exp1 + ((-1.0)^^l)*X^^exp2)*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = (X^^exp1_s - ((-1.0)^^l_s)*X^^exp2_s)*infl.Gamma[row_idx][col_idx];
			}
		},
		// even-odd
		(r, j, m, n, row_idx, col_idx) {
			immutable l = min(r, m);
			immutable exp1 = abs(r - m);
			immutable exp2 = abs(r + m);

			immutable l_s = min(r + 1, m + 1);
			immutable exp1_s = abs((r + 1) - (m + 1));
			immutable exp2_s = abs((r + 1) + (m + 1));
			if(r == 0) {
				infl.L_c[row_idx][col_idx] = X^^m*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = 0;
				//infl.L_s[row_idx][col_idx] = (X^^exp1_s - ((-1.0)^^l_s)*X^^exp2_s)*infl.Gamma[row_idx][col_idx];
			} else {
				//immutable l = min(r, m);
				//immutable exp1 = abs(r - m);
				//immutable exp2 = abs(r + m);

				infl.L_c[row_idx][col_idx] = (X^^exp1 + ((-1.0)^^l)*X^^exp2)*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = (X^^exp1_s - ((-1.0)^^l_s)*X^^exp2_s)*infl.Gamma[row_idx][col_idx];
			}
		},
		// even-even
		(r, j, m, n, row_idx, col_idx) {
			immutable l = min(r, m);
			immutable exp1 = abs(r - m);
			immutable exp2 = abs(r + m);

			immutable l_s = min(r + 1, m + 1);
			immutable exp1_s = abs((r + 1) - (m + 1));
			immutable exp2_s = abs((r + 1) + (m + 1));
			if(r == 0) {
				infl.L_c[row_idx][col_idx] = X^^m*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = 0;
				//infl.L_s[row_idx][col_idx] = (X^^exp1_s - ((-1.0)^^l_s)*X^^exp2_s)*infl.Gamma[row_idx][col_idx];
			} else {
				

				infl.L_c[row_idx][col_idx] = (X^^exp1 + ((-1.0)^^l)*X^^exp2)*infl.Gamma[row_idx][col_idx];
				//infl.L_s[row_idx][col_idx] = (X^^exp1_s - ((-1.0)^^l_s)*X^^exp2_s)*infl.Gamma[row_idx][col_idx];
			}
		}
	)(infl.Mo, infl.Me);


	iterate_whole_matrix_sin!(
		// odd-odd
		(r, j, m, n, row_idx, col_idx) {

			//writeln("oo r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
			//if(row_idx == col_idx) {
			//	infl.L_s[row_idx][col_idx] = 1;
			//} else {
			//	infl.L_s[row_idx][col_idx] = 0;
			//}
			
			// if(r == 0) {
			// 	M_s[row_idx][col_idx] = 0;
			// } else {
			// 	immutable double l = min(r, m);
			// 	immutable double exp1 = abs(r - m);
			// 	immutable double exp2 = abs(r + m);

			// 	M_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[row_idx][col_idx];
			// }

			immutable double l = min(r, m);
			immutable double exp1 = abs(r - m);
			immutable double exp2 = abs(r + m);

			immutable gamma_col_idx = col_idx + odd_states[infl.Mo][0].length;
			immutable gamma_row_idx = row_idx + odd_states[infl.Mo][0].length;

			infl.L_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*infl.Gamma[gamma_row_idx][gamma_col_idx];
		},
		// odd-even
		(r, j, m, n, row_idx, col_idx) {
			//writeln("oe r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
			if(r == 0) {
				infl.L_s[row_idx][col_idx] = 0;
			} else {
				immutable double l = min(r, m);
				immutable double exp1 = abs(r - m);
				immutable double exp2 = abs(r + m);

				//M_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[row_idx + odd_states[Mo][0].length][col_idx + total_odd_states + odd_states[Mo][0].length];

				immutable gamma_col_idx = col_idx + infl.total_odd_states + even_states[infl.Me][0].length - infl.total_odd_sin_states;
				immutable gamma_row_idx = row_idx + odd_states[infl.Mo][0].length;

				
				//writeln("gamma_row: ", gamma_row_idx, " gamma_col: ", gamma_col_idx);
				infl.L_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*infl.Gamma[gamma_row_idx][gamma_col_idx];
			}
		},
		// even-odd
		(r, j, m, n, row_idx, col_idx) {
			//writeln("eo r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
			if(r == 0) {
				//Meo_s[row_idx][col_idx] = 0;
			} else {
				immutable l = min(r, m);
				immutable exp1 = abs(r - m);
				immutable exp2 = abs(r + m);

				immutable gamma_col_idx = col_idx + odd_states[infl.Mo][0].length;
				immutable gamma_row_idx = row_idx + infl.total_odd_states + even_states[infl.Me][0].length - infl.total_odd_sin_states;

				//writeln("gamma_row: ", gamma_row_idx, " gamma_col: ", gamma_col_idx);

				infl.L_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*infl.Gamma[gamma_row_idx][gamma_col_idx];
			}
		},
		// even-even
		(r, j, m, n, row_idx, col_idx) {
			//writeln("ee r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
			if(r == 0) {
				infl.L_s[row_idx][col_idx] = 0;
			} else {
				immutable l = min(r, m);
				immutable exp1 = abs(r - m);
				immutable exp2 = abs(r + m);

				immutable gamma_col_idx = col_idx + infl.total_odd_states + even_states[infl.Me][0].length - infl.total_odd_sin_states;
				immutable gamma_row_idx = row_idx + infl.total_odd_states + even_states[infl.Me][0].length - infl.total_odd_sin_states;

				//writeln("gamma_row: ", gamma_row_idx, " gamma_col: ", gamma_col_idx);
				infl.L_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*infl.Gamma[gamma_row_idx][gamma_col_idx];
			}
		}
	)(infl.Mo, infl.Me);

	//foreach(r_idx; 0..infl.total_states) {
	foreach(r_idx, ref L_c_row; infl.L_c) {
		infl.L_c_inv[r_idx][] = L_c_row[];
	}

	foreach(r_idx, ref L_s_row; infl.L_s) {
		infl.L_s_inv[r_idx][] = L_s_row[];
	}

	//writeln("infl.L_s: ", infl.L_s);
	//write("L_s: ");
	//infl.L_s.print_matlab;
	//writeln("infl.M_s: ", infl.M_s);
	// Invert L_c matrix
	int info = 0;
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, infl.total_states, infl.total_states, infl.L_c_inv[0].ptr, infl.total_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert L_c matrix");
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, infl.total_states, infl.L_c_inv[0].ptr, infl.total_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert L_c matrix");

	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, infl.total_sin_states, infl.total_sin_states, infl.L_s_inv[0].ptr, infl.total_sin_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert L_s matrix");
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, infl.total_sin_states, infl.L_s_inv[0].ptr, infl.total_sin_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert L_s matrix");

	//write("L_s^-1: ");
	//infl.L_s_inv.print_matlab;

	// Multiply L_c_inv*M_c
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_states, infl.total_states, infl.total_states, 1.0, infl.L_c_inv[0].ptr, infl.total_states, infl.M_c[0].ptr, infl.total_states, 0.0, infl.VLM_c[0].ptr, infl.total_states);

	// Multiply L_s_inv*M_s
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, infl.total_sin_states, 1.0, infl.L_s_inv[0].ptr, infl.total_sin_states, infl.M_s[0].ptr, infl.total_sin_states, 0.0, infl.VLM_s[0].ptr, infl.total_sin_states);

	

	immutable v_inf_sin_squared = (advance_ratio*infl.sin_chi)^^2.0;
	immutable v_inf_cos = advance_ratio*infl.cos_chi;
	immutable double Vt = sqrt(v_inf_sin_squared + (v_inf_cos + infl.average_inflow)^^2.0);
	immutable double V = (v_inf_sin_squared + (v_inf_cos + infl.average_inflow)*(v_inf_cos + 2.0*infl.average_inflow))/Vt;

	//immutable double Vt = sqrt(advance_ratio^^2.0 + (axial_advance_ratio + infl.average_inflow)^^2.0);
	//immutable double V = (advance_ratio^^2.0 + (axial_advance_ratio + infl.average_inflow)*(axial_advance_ratio + 2.0*infl.average_inflow))/Vt;

	//debug writeln("Vt: ", Vt);
	//debug writeln("V: ", V);

	infl.LM_c_scratch[] = infl.VLM_c[0][];// * alpha[];
	//immutable alpha_0_1 = infl.alpha_scratch.sum;
	//infl.average_inflow = sqrt(3.0)*alpha_0_1;

	// Build the non linearity matrix
	infl.VLM_c[0][] *= Vt;
	//infl.VLM_s[0][] *= Vt;
	foreach(idx, ref vlm; infl.VLM_c[1..$]) {
		vlm[] *= V;
		//infl.VLM_s[idx + 1][] *= V;
	}

	foreach(idx, ref vlm; infl.VLM_s[0..$]) {
		vlm[] *= V;
		//infl.VLM_s[idx + 1][] *= V;
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_states, infl.total_states, infl.total_states, 1.0, infl.M_c_inv[0].ptr, infl.total_states, infl.L_c[0].ptr, infl.total_states, 0.0, infl.VLM_c_inv[0].ptr, infl.total_states);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, infl.total_sin_states, 1.0, infl.M_s_inv[0].ptr, infl.total_sin_states, infl.L_s[0].ptr, infl.total_sin_states, 0.0, infl.VLM_s_inv[0].ptr, infl.total_sin_states);

	foreach(idx, ref vlM_c_inv; infl.VLM_c_inv) {
		vlM_c_inv[0] *= 1.0/Vt;
		vlM_c_inv[1..$] *= 1.0/V;

		//infl.VLM_s_inv[idx][0] *= 1.0/Vt;
		//infl.VLM_s_inv[idx][1..$] *= 1.0/V;
	}

	foreach(idx, ref vlM_s_inv; infl.VLM_s_inv) {
		//vlM_s_inv[0] *= 1.0/Vt;
		vlM_s_inv[0..$] *= 1.0/V;

		//infl.VLM_s_inv[idx][0] *= 1.0/Vt;
		//infl.VLM_s_inv[idx][1..$] *= 1.0/V;
	}

	//debug infl.VLM_c_inv.print_matlab;
	//debug writeln(infl.adjoint_mat);
	foreach(idx, ref vlM_c_inv; infl.VLM_c_inv) {
		vlM_c_inv[] *= infl.adjoint_mat[];
		//infl.VLM_s_inv[idx][] *= infl.adjoint_mat[];
	}

	foreach(idx, ref vlM_s_inv; infl.VLM_s_inv) {
		vlM_s_inv[] *= infl.adjoint_mat_sin[];
		//infl.VLM_s_inv[idx][] *= infl.adjoint_mat[];
	}

	//debug infl.VLM_c_inv.print_matlab;

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_states, infl.total_states, infl.total_states, 1.0, infl.VLM_c_inv[0].ptr, infl.total_states, infl.VLM_c[0].ptr, infl.total_states, 0.0, infl.QS_c_mat[0].ptr, infl.total_states);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, infl.total_sin_states, 1.0, infl.VLM_s_inv[0].ptr, infl.total_sin_states, infl.VLM_s[0].ptr, infl.total_sin_states, 0.0, infl.QS_s_mat[0].ptr, infl.total_sin_states);

}

void simple_harmonic_solution(ArrayContainer AC)(HuangPetersInflowT!AC infl, double advance_ratio, double axial_advance_ratio) {

	infl.build_vlm_matrix(advance_ratio, axial_advance_ratio);

	int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, infl.total_states, infl.total_states, infl.VLM_c[0].ptr, infl.total_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert VLM_c matrix");
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, infl.total_states, infl.VLM_c[0].ptr, infl.total_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert VLM_c matrix");


	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, infl.total_sin_states, infl.total_sin_states, infl.VLM_s[0].ptr, infl.total_sin_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert VLM_s matrix");
	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, infl.total_sin_states, infl.VLM_s[0].ptr, infl.total_sin_states, infl.ipiv.ptr);
	assert(info == 0, "Failed to invert VLM_s matrix");


	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.VLM_c[0].ptr, infl.total_states, infl.tau_c.ptr, 1, 0.0, infl.state_history[infl.get_circular_index(infl.curr_state)][0..infl.total_states].ptr, 1);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, 1.0, infl.VLM_c[0].ptr, infl.total_sin_states, infl.tau_s.ptr, 1, 0.0, infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states..2*infl.total_states + infl.total_sin_states].ptr, 1);
	
	infl.alpha = infl.state_history[infl.get_circular_index(infl.curr_state)][0..infl.total_states];
	infl.beta = infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states..2*infl.total_states + infl.total_sin_states];

	foreach(i; 0..infl.total_states) {
		infl.state_history[infl.get_circular_index(infl.curr_state)][infl.total_states + i] = 0;
		foreach(k; 0..infl.total_states) {
			infl.state_history[infl.get_circular_index(infl.curr_state)][infl.total_states + i] += infl.QS_c_mat[i][k]*infl.alpha[k];
		}
	}

	foreach(i; 0..infl.total_sin_states) {
		infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states + infl.total_sin_states + i] = 0;
		foreach(k; 0..infl.total_sin_states) {
			infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states + infl.total_sin_states + i] += infl.QS_s_mat[i][k]*infl.beta[k];
		}
	}

	infl.lambda = infl.state_history[infl.get_circular_index(infl.curr_state)][infl.total_states..2*infl.total_states];

	infl.lambda_s = infl.state_history[infl.get_circular_index(infl.curr_state)][2*infl.total_states + infl.total_sin_states..$];

	infl.a[] = 0;
	infl.b[] = 0;
	// Get MD state variables
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_odd_states, infl.total_odd_states, 1.0, infl.A_inv[0].ptr, infl.total_odd_states, infl.alpha.ptr, 1, 0.0, infl.a.ptr, 1);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_odd_sin_states, infl.total_odd_sin_states, 1.0, infl.A_s_inv[0].ptr, infl.total_odd_sin_states, infl.beta.ptr, 1, 0.0, infl.b.ptr, 1);

	// Get MD adjoint variables
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_odd_states, infl.total_odd_states, 1.0, infl.A_inv[0].ptr, infl.total_odd_states, infl.lambda.ptr, 1, 0.0, infl.delta.ptr, 1);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_odd_sin_states, infl.total_odd_sin_states, 1.0, infl.A_s_inv[0].ptr, infl.total_odd_sin_states, infl.lambda_s.ptr, 1, 0.0, infl.delta_s.ptr, 1);

	infl.a[infl.total_odd_states..$] = infl.alpha[infl.total_odd_states..$];
	infl.delta[infl.total_odd_states..$] = infl.lambda[infl.total_odd_states..$];

	infl.b[infl.total_odd_sin_states..$] = infl.beta[infl.total_odd_sin_states..$];
	infl.delta_s[infl.total_odd_sin_states..$] = infl.lambda_s[infl.total_odd_sin_states..$];

}

void system_derivative(ArrayContainer AC, RIS, RS)(double[] state_dot, double[] state, double t, double dt, HuangPetersInflowT!AC infl, auto ref RIS rotor, auto ref RS rotor_state, double advance_ratio, double axial_advance_ratio) {
	
	double[] alpha = state[0..infl.total_states];
	double[] lambda = state[infl.total_states..2*infl.total_states];

	double[] beta = state[2*infl.total_states..2*infl.total_states + infl.total_sin_states];
	double[] lambda_s = state[2*infl.total_states + infl.total_sin_states..$];

	build_vlm_matrix(infl, advance_ratio, axial_advance_ratio);

	// scratch * alpha = V*L_c_inv*M_c*alpha
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.VLM_c[0].ptr, infl.total_states, alpha.ptr, 1, 0.0, infl.alpha_scratch.ptr, 1);

	// scratch * beta = V*L_s_inv*M_s*beta
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, 1.0, infl.VLM_s[0].ptr, infl.total_sin_states, beta.ptr, 1, 0.0, infl.beta_scratch.ptr, 1);

	// scratch * delta = V*L_c_inv*M_c*delta
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.VLM_c[0].ptr, infl.total_states, lambda.ptr, 1, 0.0, infl.delta_scratch.ptr, 1);

	// scratch * delta_s = V*L_s_inv*M_s*delta_s
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, 1.0, infl.VLM_s[0].ptr, infl.total_sin_states, lambda_s.ptr, 1, 0.0, infl.delta_s_scratch.ptr, 1);
	
	// ADJ*tau_c
	//infl.tau_c_scratch[] = infl.tau_c[]*infl.adjoint_mat[];
	//infl.tau_c_scratch[] = infl.tau_c[]*infl.adjoint_mat[] - infl.delta_scratch[];
	//infl.tau_c_scratch[] = infl.delta_scratch[] - infl.tau_c[]*infl.adjoint_mat[];

	infl.alpha_scratch[] = infl.tau_c[] - infl.alpha_scratch[];
	infl.beta_scratch[] = infl.tau_s[] - infl.beta_scratch[];
	//infl.alpha_scratch[] = infl.alpha_scratch[] - infl.tau_c[];

	// M_c_inv*D*alpha_scratch = M_c_inv*D*(tau_c - V*L_c_inv*M_c*alpha)
	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.M_c_inv_D[0].ptr, infl.total_states, infl.alpha_scratch.ptr, 1, 0.0, state_dot[0..infl.total_states].ptr, 1);

	cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_sin_states, infl.total_sin_states, 1.0, infl.M_s_inv_D[0].ptr, infl.total_sin_states, infl.beta_scratch.ptr, 1, 0.0, state_dot[2*infl.total_states..2*infl.total_states + infl.total_sin_states].ptr, 1);

	// M_c_inv*D*tau_c_scratch = M_c_inv*D*(V*L_c_inv*M_c*delta - ADJ*tau_c)
	//cblas_dgemv(CblasRowMajor, CblasNoTrans, infl.total_states, infl.total_states, 1.0, infl.M_c_inv_D[0].ptr, infl.total_states, infl.tau_c_scratch.ptr, 1, 0.0, state_dot[infl.total_states..$].ptr, 1);
}

alias HuangPetersInflow = HuangPetersInflowT!(ArrayContainer.none);

class HuangPetersInflowT(ArrayContainer AC = ArrayContainer.none) {// : Inflow { 

	alias RG = RotorGeometryT!AC;

	private immutable long Mo;
	private immutable long Me;
	private long N;

	private Chunk[] mpsi_buff;
	private Chunk[] sin_mpsi_buff;

	private double[][] M_c; // Constant.
	private double[][] L_c_inv; // Modification on Gamma by skew angle then inverted each timestep.
	private double[][] L_c; // Modification on Gamma by skew angle then inverted each timestep.
	private double[][] M_s; // Constant.
	private double[][] L_s_inv; // Modification on Gamma by skew angle then inverted each timestep.
	private double[][] L_s; // Modification on Gamma by skew angle then inverted each timestep.
	private double[][] Gamma; // Constant throughout.
	private double[][] M_c_inv_D; // Damping. Constant
	private double[][] M_s_inv_D; // Damping. Constant
	private double[][] VLM_c; // This holds the result of the matrix mults each timestep.
	private double[][] VLM_c_inv; // This holds the result of the matrix mults each timestep.
	private double[][] VLM_s; // This holds the result of the matrix mults each timestep.
	private double[][] VLM_s_inv; // This holds the result of the matrix mults each timestep.
	double[][] QS_c_mat; // This holds the result of the matrix mults each timestep.
	double[][] QS_s_mat; // This holds the result of the matrix mults each timestep.
	private double[][] A_inv;
	private double[][] A;
	private double[][] A_s;
	private double[][] A_s_inv;
	private double[][] M_c_inv;
	private double[][] M_s_inv;
	private double[][] D;
	private double[][] D_s;

	//Chunk[] V_buffer;

	private int[] ipiv; // Used by lapack for matrix inversion. Keep it around instead of re-allocing each update

	private ptrdiff_t z_length;
	private ptrdiff_t time_history;
	private ptrdiff_t curr_state;
	private Chunk[] time_delay_alpha;
	private Chunk[] time_delay_alpha_2;
	private Chunk[] time_delay_a;
	private Chunk[] time_delay_a_2;
	private Chunk[] time_delay_beta;
	private Chunk[] time_delay_beta_2;
	private Chunk[] time_delay_b;
	private Chunk[] time_delay_b_2;
	private double[][] state_history;
	private double[] times;
	double[] tau_c; // Pressure coefficients
	double[] tau_s; // Pressure coefficients
	private double[] tau_c_scratch;
	private double[] tau_s_scratch;
	double[] alpha; // NH state variables
	double[] beta; // NH state variables
	private double[] a; // MD state variables
	private double[] b; // MD state variables
	private double[] alpha_scratch;
	private double[] beta_scratch;
	double[] delta; // MD adjoint variables
	double[] lambda; // NH adjoint variables
	double[] delta_s; // MD adjoint variables
	double[] lambda_s; // NH adjoint variables
	private double[] delta_scratch;
	private double[] delta_s_scratch;
	private double[] adjoint_mat;
	private double[] adjoint_mat_sin;
	private double[] LM_c_scratch;
	//private double[] LM_s_scratch;
	private double[][] K_table;
	private double[] average_inflow_array;
	//private size_t ai_idx;

	int total_states;
	int total_odd_states;
	int total_even_states;

	int total_sin_states;
	int total_odd_sin_states;
	int total_even_sin_states;

	private double average_inflow;

	private double[][] P_coefficients;
	private double[][] P_coefficients_nh;
	private Chunk[][] Qmn_bar;
	private Chunk[][] P_md_mn_bar;
	private Chunk[][] P_nh_mn_bar;
	
	alias Integrator = ForwardEuler!double;
	//alias Integrator = RK4!double;

	Integrator integrator;

	private Chunk[] blade_scratch;
	private Chunk[] blade_scratch_s;

	private RG* rotor;

	private double advance_ratio;
	private double axial_advance_ratio;
	private double omega;

	private double chi;
	private double sin_chi;
	private double cos_chi;
	private double tan_chi;
	private double cot_chi;

	size_t n_r = 16;
	size_t n_psi = 6;
	double v_0;

	Chunk[] contraction_array;
	double[] contraction_array_alias;
	Chunk[] contraction_z_array;
	double[] contraction_z_array_alias;

	bool contraction_mapping = false;
	bool debug_coords = false;

	double[] zero;
	Chunk[] c_zero;

	@nogc private ptrdiff_t get_circular_index(ptrdiff_t idx) {
		return ((idx % time_history) + time_history) % time_history;
	}

	this(long _Mo, long _Me, RG* _rotor, double dt) {
		size_t len = round(2.0*PI/(dt*235.325)).to!size_t;
		//ai_idx = 0;
		average_inflow_array = new double[len];
		average_inflow_array[] = 0;
		rotor = _rotor;
		size_t num_chunks = rotor.blades[0].chunks.length;

		openblas_set_num_threads(1);

		z_length = 48;
		contraction_z_array_alias = linspace(-5.0, 5.0, z_length);
		contraction_z_array = cast(Chunk[])contraction_z_array_alias.ptr[0..z_length];

		contraction_array_alias = new double[z_length];
		contraction_array = cast(Chunk[])contraction_array_alias.ptr[0..z_length];

		blade_scratch = new Chunk[num_chunks];
		blade_scratch_s = new Chunk[num_chunks];

		Mo = _Mo;
		Me = _Me;

		mpsi_buff = new Chunk[max(Mo, Me) + 1];
		sin_mpsi_buff = new Chunk[max(Mo, Me) + 1];

		N = 0;

		total_states = 0;
		total_even_states = 0;
		total_odd_states = 0;

		iterate_odds!((m, n, idx) {

			size_t num_coefficients = n - m + 2;
			P_coefficients ~= new double[num_coefficients];
			P_coefficients[$-1][0] = (-1.0)^^m.to!double*((-1.0)^^m.to!double)*(2.0^^n.to!double);
			foreach(k; m..n + 1) {
				P_coefficients[$-1][k - m + 1] = (k.factorial/(k - m).factorial)*(n.falling_factorial(k)/k.factorial)*((0.5*(n.to!double + k.to!double - 1.0)).falling_factorial(n)/n.factorial );
			}

			N = max(n + 1, N);
			immutable rho = sqrt(1.0/(2.0*n.to!double + 1.0)*(n + m).factorial/(n - m).factorial);
			P_coefficients[$-1][0] /= rho;

			//writeln("idx: ", idx, "; m: ", m, ", n: ", n);
			total_states++;
			total_odd_states++;
		})(Mo, 0);

		iterate_evens!((m, n, idx) {
			size_t num_coefficients = n - m + 2;
			P_coefficients ~= new double[num_coefficients];
			P_coefficients[$-1][0] = (-1.0)^^m.to!double*((-1.0)^^m.to!double)*(2.0^^n.to!double);
			foreach(k; m..n + 1) {
				P_coefficients[$-1][k - m + 1] = (k.factorial/(k - m).factorial)*(n.falling_factorial(k)/k.factorial)*((0.5*(n.to!double + k.to!double - 1.0)).falling_factorial(n)/n.factorial );
			}

			N = max(n + 1, N);

			immutable rho = sqrt(1.0/(2.0*n.to!double + 1.0)*(n + m).factorial/(n - m).factorial);
			P_coefficients[$-1][0] /= rho;

			//writeln("idx: ", idx, "; m: ", m, ", n: ", n);
			total_states++;
			total_even_states++;
		})(Me, total_odd_states);

		iterate_odds!((r, j, idx) {

			size_t num_coefficients = j - r + 2;
			P_coefficients_nh ~= new double[num_coefficients];
			P_coefficients_nh[$-1][0] = sqrt((2.0*j.to!double + 1.0)*H(r, j));

			foreach(qidx, q; iota(r, j, 2).enumerate) {
				P_coefficients_nh[$-1][q - r + 1] = (-1.0)^^(0.5*(q.to!double - r.to!double))*double_factorial(j + q);
				P_coefficients_nh[$-1][q - r + 1] /= double_factorial(q - r)*double_factorial(q + r)*double_factorial(j - q - 1);
			}
		})(Mo, 0);

		iterate_evens!((r, j, idx) {

			size_t num_coefficients = j - r + 2;
			P_coefficients_nh ~= new double[num_coefficients];
			P_coefficients_nh[$-1][0] = sqrt((2.0*j.to!double + 1.0)*H(r, j));

			foreach(qidx, q; iota(r, j, 2).enumerate) {
				P_coefficients_nh[$-1][q - r + 1] = (-1.0)^^(0.5*(q.to!double - r.to!double))*double_factorial(j + q);
				P_coefficients_nh[$-1][q - r + 1] /= double_factorial(q - r)*double_factorial(q + r)*double_factorial(j - q - 1);
			}
		})(Me, 0);

		total_even_sin_states = total_even_states - even_states[Me][0].length.to!int;
		total_odd_sin_states = total_odd_states - odd_states[Mo][0].length.to!int;

		total_sin_states = total_odd_sin_states + total_even_sin_states;

		debug writeln("total_states: ", total_states);
		debug writeln("total_odd_states: ", total_odd_states);
		debug writeln("total_even_states: ", total_even_states);
		debug writeln("total_sin_states: ", total_sin_states);
		debug writeln("total_odd_sin_states: ", total_odd_sin_states);
		debug writeln("total_even_sin_states: ", total_even_sin_states);
		
		integrator = new Integrator(2*total_states + 2*total_sin_states);

		//V_buffer = new Chunk[total_states];

		double[][] Meo = allocate_dense(total_even_states, total_odd_states);
		double[][] Meo_new = allocate_dense(total_even_states, total_odd_states);
		double[][] Meo_s = allocate_dense(total_even_sin_states, total_odd_sin_states);
		double[][] Meo_s_new = allocate_dense(total_even_sin_states, total_odd_sin_states);
		K_table = allocate_dense(max(Me, Mo) + 1, N);
		D = allocate_dense(total_states, total_states);
		D_s = allocate_dense(total_sin_states, total_sin_states);
		QS_c_mat = allocate_dense(total_states, total_states);
		M_c_inv = allocate_dense(total_states, total_states);
		QS_s_mat = allocate_dense(total_sin_states, total_sin_states);
		M_s_inv = allocate_dense(total_sin_states, total_sin_states);
		A_inv = allocate_dense(total_odd_states, total_odd_states);
		A_s_inv = allocate_dense(total_odd_sin_states, total_odd_sin_states);
		A_s = allocate_dense(total_odd_sin_states, total_odd_sin_states);
		A = allocate_dense(total_odd_states, total_odd_states);
		VLM_c = allocate_dense(total_states, total_states);
		VLM_c_inv = allocate_dense(total_states, total_states);
		M_c = allocate_dense(total_states, total_states);
		M_c_inv_D = allocate_dense(total_states, total_states);
		L_c_inv = allocate_dense(total_states, total_states);
		L_c = allocate_dense(total_states, total_states);
		VLM_s = allocate_dense(total_sin_states, total_sin_states);
		VLM_s_inv = allocate_dense(total_sin_states, total_sin_states);
		M_s = allocate_dense(total_sin_states, total_sin_states);
		M_s_inv_D = allocate_dense(total_sin_states, total_sin_states);
		L_s_inv = allocate_dense(total_sin_states, total_sin_states);
		L_s = allocate_dense(total_sin_states, total_sin_states);
		Gamma = allocate_dense(total_states, total_states);

		time_history = 1_000_000;
		state_history = allocate_dense(time_history, 2*total_states + 2*total_sin_states);
		
		//Qmn_bar = new Chunk[][](max(Me, Mo) + 1, N);
		Qmn_bar = allocate_dense_chunk(max(Me, Mo) + 1, N);
		P_nh_mn_bar = allocate_dense_chunk(max(Me, Mo) + 1, N);
		P_md_mn_bar = allocate_dense_chunk(max(Me, Mo) + 1, N);

		times = new double[time_history];
		
		zero = new double[total_sin_states];
		c_zero = new Chunk[total_sin_states];
		time_delay_alpha = new Chunk[2*total_states];
		time_delay_alpha_2 = new Chunk[2*total_states];
		time_delay_a = new Chunk[2*total_states];
		time_delay_a_2 = new Chunk[2*total_states];

		time_delay_beta = new Chunk[2*total_sin_states];
		time_delay_beta_2 = new Chunk[2*total_sin_states];
		time_delay_b = new Chunk[2*total_sin_states];
		time_delay_b_2 = new Chunk[2*total_sin_states];

		a = new double[total_states];
		b = new double[total_sin_states];
		delta = new double[total_states];
		delta_s = new double[total_sin_states];
		alpha_scratch = new double[total_states];
		beta_scratch = new double[total_sin_states];
		delta_scratch = new double[total_states];
		delta_s_scratch = new double[total_sin_states];
		tau_c = new double[total_states];
		tau_c_scratch = new double[total_states];
		tau_s = new double[total_sin_states];
		tau_s_scratch = new double[total_sin_states];
		adjoint_mat = new double[total_states];
		adjoint_mat_sin = new double[total_sin_states];
		LM_c_scratch = new double[total_states];
		//LM_s_scratch = new double[total_sin_states];

		zero[] = 0;
		foreach(ref ch; c_zero) {
			ch[] = 0;
		}
		alpha_scratch[] = 0;
		beta_scratch[] = 0;
		delta_scratch[] = 0;
		delta_s_scratch[] = 0;
		tau_c[] = 0;
		tau_c_scratch[] = 0;
		tau_s[] = 0;
		tau_s_scratch[] = 0;
		adjoint_mat[] = 0;
		adjoint_mat_sin[] = 0;

		state_history.zero_matrix;

		alpha = state_history[0][0..total_states];
		lambda = state_history[0][total_states..2*total_states];

		beta = state_history[0][2*total_states..2*total_states + total_sin_states];
		lambda_s = state_history[0][2*total_states + total_sin_states..$];

		A_inv.zero_matrix;
		A_s_inv.zero_matrix;
		Meo.zero_matrix;
		Meo_new.zero_matrix;
		Meo_s.zero_matrix;
		Meo_s_new.zero_matrix;
		M_c_inv.zero_matrix;
		M_c.zero_matrix;
		M_s_inv.zero_matrix;
		M_s.zero_matrix;
		D.zero_matrix;
		D_s.zero_matrix;
		Gamma.zero_matrix;
		L_c_inv.zero_matrix;
		VLM_c.zero_matrix;
		VLM_c_inv.zero_matrix;
		QS_c_mat.zero_matrix;
		L_s_inv.zero_matrix;
		VLM_s.zero_matrix;
		VLM_s_inv.zero_matrix;
		QS_s_mat.zero_matrix;

		iterate_even_odd!(
			(r, j, m, n, row_idx, col_idx) {
				if(r == m) {
					if((j == n - 1) || (j == n + 1)) {
						
						Meo[row_idx][col_idx] = 1;
						Meo[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
					}
				}
			}
		)(Mo, Me, 0, 0);

		iterate_odd_odd!((r, j, m, n, row_idx, col_idx) {
			if(r == m) {
				A_inv[row_idx][col_idx] = 2.0*(-1.0)^^(0.5*(n + j - 2*r).to!double);
				A_inv[row_idx][col_idx] *= sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
				A_inv[row_idx][col_idx] /= sqrt(H(r, n)*H(r, j))*(n + j).to!double*(n.to!double + j.to!double + 2.0)*((n.to!double - j.to!double)^^2.0 - 1.0);
			}
		})(Mo, 0, 0);

		foreach(r_idx; 0..total_odd_states) {
			A[r_idx][] = A_inv[r_idx][];
		}

		iterate_odd_odd_sin!((r, j, m, n, row_idx, col_idx) {
			if(r == m) {
				A_s_inv[row_idx][col_idx] = 2.0*(-1.0)^^(0.5*(n + j - 2*r).to!double);
				A_s_inv[row_idx][col_idx] *= sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
				A_s_inv[row_idx][col_idx] /= sqrt(H(r, n)*H(r, j))*(n + j).to!double*(n.to!double + j.to!double + 2.0)*((n.to!double - j.to!double)^^2.0 - 1.0);
			}
		})(Mo, 0, 0);

		foreach(r_idx; 0..total_odd_sin_states) {
			A_s[r_idx][] = A_s_inv[r_idx][];
		}

		int info = 0;
		ipiv = new int[total_odd_states];

		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_odd_states, total_odd_states, A_inv[0].ptr, total_odd_states, ipiv.ptr);
		assert(info == 0, "Failed to invert A matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_odd_states, A_inv[0].ptr, total_odd_states, ipiv.ptr);
		assert(info == 0, "Failed to invert A matrix");

		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_odd_sin_states, total_odd_sin_states, A_s_inv[0].ptr, total_odd_sin_states, ipiv.ptr);
		assert(info == 0, "Failed to invert A_s matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_odd_sin_states, A_s_inv[0].ptr, total_odd_sin_states, ipiv.ptr);
		assert(info == 0, "Failed to invert A_s matrix");

		//A_s_inv.print_matlab;

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, total_even_states, total_odd_states, total_odd_states, 1.0, Meo[0].ptr, total_odd_states, A_inv[0].ptr, total_odd_states, 0.0, Meo_new[0].ptr, total_odd_states);

		iterate_whole_matrix!(
			// odd-odd
			(r, j, m, n, row_idx, col_idx) {
				//writeln("oo r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				if(row_idx == col_idx) {
					M_c[row_idx][col_idx] = 1;
				}

				if(r == m) {
					D[row_idx][col_idx] = kronecker_delta(j, n)/K(m, n);
				}

				if(((r + m) % 2 != 0) && ((j == (n - 1)) || (j == (n + 1)))) {
					Gamma[row_idx][col_idx] = sgn(r - m);
					Gamma[row_idx][col_idx] /= sqrt(K(m, n)*K(r, j))*sqrt((2*n + 1).to!double*(2*j + 1).to!double);
				}

				if((r + m) % 2 == 0) {
					Gamma[row_idx][col_idx] = 2.0*(-1.0)^^(0.5*(n + j - 2.0*r));
					Gamma[row_idx][col_idx] *= sqrt((2.0*n + 1.0)*(2.0*j + 1.0));
					Gamma[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*(n + j).to!double*(n + j + 2).to!double*((n - j).to!double^^2.0 - 1.0);
				}
			},
			// odd-even
			(r, j, m, n, row_idx, col_idx) {
				//writeln("oe r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				if(r == m) {
					if((j == n - 1) || (j == n + 1)) {
						M_c[row_idx][col_idx] = 1;
						M_c[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*sqrt((2*n + 1).to!double*(2*j + 1).to!double);
					}

					D[row_idx][col_idx] = 2.0*sqrt((2*j + 1).to!double*(2*n + 1).to!double);
					D[row_idx][col_idx] *= (-1.0)^^(0.5*(j + 3*n - 1).to!double);
					D[row_idx][col_idx] /= PI*sqrt(H(m, n)*H(m, j))*(j + n + 1).to!double*(j - n).to!double;
				}

				if((r + m) % 2 != 0) {
					Gamma[row_idx][col_idx] = sgn(r - m)*4.0*(-1.0)^^(0.5*(3.0*n + j + 2.0*m - 2.0*r));
					Gamma[row_idx][col_idx] *= sqrt((2.0*n + 1.0)*(2.0*j + 1.0));
					Gamma[row_idx][col_idx] /= PI*sqrt(H(m, n)*H(r, j))*(n + j)*(n + j + 2.0)*((n - j)^^2.0 - 1.0);
				}

				if(((r + m) % 2 == 0) && ((j == (n - 1)) || (j == (n + 1)))) {
					Gamma[row_idx][col_idx] = 1.0;
					Gamma[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*sqrt((2*n + 1).to!double*(2*j + 1).to!double);
				}
			},
			// even-odd
			(r, j, m, n, row_idx, col_idx) {
				//writeln("eo r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				if(r == m) {
					D[row_idx][col_idx] = 2.0*sqrt((2.0*j.to!double + 1)*(2.0*n.to!double + 1));
					D[row_idx][col_idx] *= (-1.0)^^(0.5*(j + 3*n - 1).to!double);
					D[row_idx][col_idx] /= PI*sqrt(H(m, n)*H(m, j))*(j + n + 1).to!double*(j - n).to!double;
				}

				M_c[row_idx][col_idx] = Meo_new[row_idx - total_odd_states][col_idx];

				if((r + m) % 2 != 0) {
					Gamma[row_idx][col_idx] = sgn(r - m)*4.0*(-1.0)^^(0.5*(3*n + j + 2*m - 2*r).to!double);
					Gamma[row_idx][col_idx] *= sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
					Gamma[row_idx][col_idx] /= PI*sqrt(H(m, n)*H(r, j))*(n + j).to!double*(n + j + 2).to!double*((n - j).to!double^^2.0 - 1.0);
				}

				if(((r + m) % 2 == 0) && ((j == (n - 1)) || (j == (n + 1)))) {
					Gamma[row_idx][col_idx] = 1.0;
					Gamma[row_idx][col_idx] /= sqrt(H(m, n)*H(r, j))*sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
				}
			},
			// even-even
			(r, j, m, n, row_idx, col_idx) {
				//writeln("ee r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				if(r == m) {
					M_c[row_idx][col_idx] = 8.0*(-1.0)^^(0.5*(n + j - 2*m + 2).to!double);
					M_c[row_idx][col_idx] *= sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
					M_c[row_idx][col_idx] /= PI*PI*sqrt(H(m, n)*H(m, j))*(n + j).to!double*(n + j + 2).to!double*((n - j).to!double^^2.0 - 1.0);

					D[row_idx][col_idx] = kronecker_delta(j, n)/K(m, n);
				}

				if(((r + m) % 2 != 0) && ((j == (n - 1)) || (j == (n + 1)))) {
					Gamma[row_idx][col_idx] = sgn(r - m);
					Gamma[row_idx][col_idx] /= sqrt(K(m, n)*K(r, j))*sqrt((2.0*n.to!double + 1.0)*(2.0*j.to!double + 1.0));
				}

				if((r + m) % 2 == 0) {
					Gamma[row_idx][col_idx] = 8.0*(-1.0)^^(0.5*(n + j - 2.0*r + 2.0));
					Gamma[row_idx][col_idx] *= sqrt((2*n + 1).to!double*(2*j + 1).to!double);
					Gamma[row_idx][col_idx] /= PI*PI*sqrt(H(m, n)*H(r, j))*(n + j).to!double*(n + j + 2).to!double*((n - j).to!double^^2.0 - 1.0);
				}
			}
		)(Mo, Me);

		//writeln;

		double X = 0;
		iterate_even_odd_sin!(
			(r, j, m, n, row_idx, col_idx) {
				//writeln("eo r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				if(r == 0) {
					Meo_s[row_idx][col_idx] = 0;
				} else {
					immutable l = min(r, m);
					immutable exp1 = abs(r - m);
					immutable exp2 = abs(r + m);

					immutable gamma_col_idx = col_idx + odd_states[Mo][0].length;
					immutable gamma_row_idx = row_idx + total_odd_states + even_states[Me][0].length;

					//writeln("gamma_row: ", gamma_row_idx, " gamma_col: ", gamma_col_idx);

					Meo_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[gamma_row_idx][gamma_col_idx];
				}
			}
		)(Mo, Me, 0, 0);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, total_even_sin_states, total_odd_sin_states, total_odd_sin_states, 1.0, Meo_s[0].ptr, total_odd_sin_states, A_s_inv[0].ptr, total_odd_sin_states, 0.0, Meo_s_new[0].ptr, total_odd_sin_states);

		//writeln;
		//writeln;

		iterate_whole_matrix_sin!(
			// odd-odd
			(r, j, m, n, row_idx, col_idx) {

				//writeln("oo r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);

				immutable gamma_col_idx = col_idx + odd_states[Mo][0].length;
				immutable gamma_row_idx = row_idx + odd_states[Mo][0].length;

				//writeln("gamma_row: ", gamma_row_idx, " gamma_col: ", gamma_col_idx);

				if(row_idx == col_idx) {
					M_s[row_idx][col_idx] = 1;
				} else {
					M_s[row_idx][col_idx] = 0;
				}
				
				// if(r == 0) {
				// 	M_s[row_idx][col_idx] = 0;
				// } else {
				// 	immutable double l = min(r, m);
				// 	immutable double exp1 = abs(r - m);
				// 	immutable double exp2 = abs(r + m);

				// 	M_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[row_idx][col_idx];
				// }
			},
			// odd-even
			(r, j, m, n, row_idx, col_idx) {
				//writeln("oe r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				if(r == 0) {
					M_s[row_idx][col_idx] = 0;
				} else {
					immutable double l = min(r, m);
					immutable double exp1 = abs(r - m);
					immutable double exp2 = abs(r + m);

					//M_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[row_idx + odd_states[Mo][0].length][col_idx + total_odd_states + odd_states[Mo][0].length];

					immutable gamma_col_idx = col_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;
					immutable gamma_row_idx = row_idx + odd_states[Mo][0].length;

					
					//writeln("gamma_row: ", gamma_row_idx, " gamma_col: ", gamma_col_idx);
					M_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[gamma_row_idx][gamma_col_idx];
				}
			},
			// even-odd
			(r, j, m, n, row_idx, col_idx) {
				//writeln("eo r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				M_s[row_idx][col_idx] = Meo_s_new[row_idx - total_odd_sin_states][col_idx];
				//M_s[row_idx][col_idx] = 0;
				// if(r == 0) {
				// 	Meo_s[row_idx][col_idx] = 0;
				// } else {
				// 	immutable l = min(r, m);
				// 	immutable exp1 = abs(r - m);
				// 	immutable exp2 = abs(r + m);

				// 	Meo_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[row_idx][col_idx];
				// }
			},
			// even-even
			(r, j, m, n, row_idx, col_idx) {
				//writeln("ee r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				if(r == 0) {
					M_s[row_idx][col_idx] = 0;
				} else {
					immutable l = min(r, m);
					immutable exp1 = abs(r - m);
					immutable exp2 = abs(r + m);

					immutable gamma_col_idx = col_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;
					immutable gamma_row_idx = row_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;

					//writeln("gamma_row: ", gamma_row_idx, " gamma_col: ", gamma_col_idx);
					M_s[row_idx][col_idx] = (X^^exp1 - ((-1.0)^^l)*X^^exp2)*Gamma[gamma_row_idx][gamma_col_idx];
				}
			}
		)(Mo, Me);

		//Gamma[0][0] 
		// We'll leave it at this size, as this is the same size as the L_c matrix will want.
		ipiv = new int[total_states];

		// We need both M_c and M_c^-1 for timestep updates so copy M_c to M_c_inv and invert.
		foreach(r_idx; 0..total_states) {
			M_c_inv[r_idx][] = M_c[r_idx][];
		}

		foreach(r_idx; 0..total_sin_states) {
			M_s_inv[r_idx][] = M_s[r_idx][];
		}

		//write("M_s: ");
		//M_s.print_matlab;
		
		iterate_whole_matrix_sin!(
			// odd-odd
			(r, j, m, n, row_idx, col_idx) {
				immutable gamma_col_idx = col_idx + odd_states[Mo][0].length;
				immutable gamma_row_idx = row_idx + odd_states[Mo][0].length;

				D_s[row_idx][col_idx] = D[gamma_row_idx][gamma_col_idx];
			},
			// odd-even
			(r, j, m, n, row_idx, col_idx) {
				//writeln("oe r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				immutable gamma_col_idx = col_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;
				immutable gamma_row_idx = row_idx + odd_states[Mo][0].length;

				D_s[row_idx][col_idx] = D[gamma_row_idx][gamma_col_idx];
			},
			// even-odd
			(r, j, m, n, row_idx, col_idx) {
				//writeln("eo r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				immutable gamma_col_idx = col_idx + odd_states[Mo][0].length;
				immutable gamma_row_idx = row_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;

				D_s[row_idx][col_idx] = D[gamma_row_idx][gamma_col_idx];
			},
			// even-even
			(r, j, m, n, row_idx, col_idx) {
				//writeln("ee r: ", r, " m: ", m, " j: ", j, " n: ", n, " row: ", row_idx, " col: ", col_idx);
				immutable gamma_col_idx = col_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;
				immutable gamma_row_idx = row_idx + total_odd_states + even_states[Me][0].length - total_odd_sin_states;

				//writeln("gamma_row: ", gamma_row_idx, " gamma_col: ", gamma_col_idx);
				D_s[row_idx][col_idx] = D[gamma_row_idx][gamma_col_idx];
			}
		)(Mo, Me);

		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_states, total_states, M_c_inv[0].ptr, total_states, ipiv.ptr);
		assert(info == 0, "Failed to invert M_c matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_states, M_c_inv[0].ptr, total_states, ipiv.ptr);
		assert(info == 0, "Failed to invert M_c matrix");


		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, total_sin_states, total_sin_states, M_s_inv[0].ptr, total_sin_states, ipiv.ptr);
		assert(info == 0, "Failed to invert M_s matrix");
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, total_sin_states, M_s_inv[0].ptr, total_sin_states, ipiv.ptr);
		assert(info == 0, "Failed to invert M_s matrix");

		//write("\nM_c_inv: ");
		//M_c_inv.print_matlab;

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, total_states, total_states, total_states, 1.0, M_c_inv[0].ptr, total_states, D[0].ptr, total_states, 0.0, M_c_inv_D[0].ptr, total_states);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, total_sin_states, total_sin_states, total_sin_states, 1.0, M_s_inv[0].ptr, total_sin_states, D_s[0].ptr, total_sin_states, 0.0, M_s_inv_D[0].ptr, total_sin_states);

		iterate_odds!((m, n, idx) {
			K_table[m][n] = K(m, n);
			adjoint_mat[idx] = (-1.0)^^(n.to!double + 1.0);
		})(Mo, 0);

		iterate_evens!((m, n, idx) {
			K_table[m][n] = K(m, n);
			adjoint_mat[idx] = (-1.0)^^(n.to!double + 1.0);
		})(Me, total_odd_states);


		iterate_odds_sin!((m, n, idx) {
			adjoint_mat_sin[idx] = (-1.0)^^(n.to!double + 1.0);
		})(Mo, 0);

		iterate_evens_sin!((m, n, idx) {
			adjoint_mat_sin[idx] = (-1.0)^^(n.to!double + 1.0);
		})(Me, total_odd_sin_states);

		foreach(n; 1..Qmn_bar[0].length - 1) {
			immutable long m = 0;
			K_table[m][n] = K(m, n);
		}

		foreach(m; 0..Qmn_bar.length - 1) {
			foreach(n; 1..Qmn_bar[m].length) {
				K_table[m][n] = K(m, n);
			}
		}		

		curr_state = 0;
		times[] = 0;
		average_inflow = 0.01;
		chi = 0;
		tau_c[0] = 0.001;
		tau_s[0] = 0.001;
		tau_s[] = [-0.00155803, 0.000844799, -0.00282596, 0.000751855, -0.00304, -0.00204455, -0.000332659, -0.000887072];
		//tau_s[0] = 0.001;

		simple_harmonic_solution(this, 0.0, 0.0);

		//A_s_inv.print_matlab;

	}

	@nogc auto find_bracket(double t) {
		
		auto ordered_times = times[get_circular_index(curr_state) + 1..$].chain(times[0..get_circular_index(curr_state) + 1]);
		auto delta_b = time_history - get_circular_index(curr_state);
		ptrdiff_t l = 0;
		ptrdiff_t R = time_history - 1;

		while(l <= R) {
			auto m = (l + R)/2;

			if(m == (time_history - 1)) {
				if(ordered_times[m] != t) {
					return -1;
				} else {
					return m >= delta_b ? m - delta_b + 1 : m - get_circular_index(curr_state) + 1;
				}
			} else if ((ordered_times[m] <= t) && (ordered_times[m + 1] > t)) {
				return m >= delta_b ? m - delta_b + 1 : m - get_circular_index(curr_state) + 1;
			} else if(ordered_times[m] <= t) {
				l = m + 1;
			} else if(ordered_times[m] > t) {
				R = m - 1;
			}
		}

		return -1;
	}

	@nogc auto find_z_bracket(double z) {
		ptrdiff_t l = 0;
		ptrdiff_t R = z_length - 1;

		while(l <= R) {
			auto m = (l + R)/2;

			if(m == (z_length - 1)) {
				return m;
			} else if ((contraction_z_array_alias[m] <= z) && (contraction_z_array_alias[m + 1] > z)) {
				return m;
			} else if(contraction_z_array_alias[m] <= z) {
				l = m + 1;
			} else if(contraction_z_array_alias[m] > z) {
				R = m - 1;
			}
		}

		return 0;
	}

	@nogc void interpolate_to_time(immutable Chunk t, ref Chunk[] time_delay_alpha_buffer, ref Chunk[] time_delay_a_buffer, ref Chunk[] time_delay_beta_buffer, ref Chunk[] time_delay_b_buffer) {

		ptrdiff_t[chunk_size] upper_bounds;
		ptrdiff_t[chunk_size] lower_bounds;

		bool all_nan = t[].map!(a => a.isNaN).fold!((res, a) => res && a)(true);
		if(all_nan) {
			return;
		}

		foreach(idx, ref _t; t[].map!(a => a.isNaN ? 0 : a).enumerate) {
			lower_bounds[idx] = find_bracket(_t);
			upper_bounds[idx] = get_circular_index(lower_bounds[idx] + 1);
		}

		foreach(idx; 0..2*total_states) {
			Chunk t_lower;
			Chunk t_upper;
			Chunk state_u;
			Chunk state_l;
			foreach(c_idx; 0..chunk_size) {
				if(lower_bounds[c_idx] < 0) {
					state_u[c_idx] = 0;
					state_l[c_idx] = 0;
					t_lower[c_idx] = 0;
					t_upper[c_idx] = 1;
				} else {
					t_lower[c_idx] = times[lower_bounds[c_idx]];
					t_upper[c_idx] = times[upper_bounds[c_idx]];

					state_u[c_idx] = state_history[upper_bounds[c_idx]][idx];
					state_l[c_idx] = state_history[lower_bounds[c_idx]][idx];
				}
			}

			immutable Chunk dt = t_upper[] - t_lower[];
			immutable Chunk curr_t = times[get_circular_index(curr_state)];
			immutable Chunk t_delta = curr_t[] - t[];
			immutable Chunk t_delta_abs = abs(t_delta);
			time_delay_alpha_buffer[idx][] = 0;
			
			immutable Chunk time_delta = t[] - t_lower[];
			immutable Chunk state_delta = state_u[] - state_l[];
			time_delay_alpha_buffer[idx][] = state_l[] + time_delta[]*state_delta[]/dt[];

			foreach(c_idx; 0..chunk_size) {
				if(lower_bounds[c_idx] >= 0) {
					if(t_delta_abs[c_idx] <= dt[c_idx]) {
						time_delay_alpha_buffer[idx][c_idx] = state_history[get_circular_index(curr_state)][idx];
					}
				}
			}
		}

		foreach(idx; 0..2*total_sin_states) {
			Chunk t_lower;
			Chunk t_upper;
			Chunk state_u;
			Chunk state_l;
			foreach(c_idx; 0..chunk_size) {
				if(lower_bounds[c_idx] < 0) {
					state_u[c_idx] = 0;
					state_l[c_idx] = 0;
					t_lower[c_idx] = 0;
					t_upper[c_idx] = 1;
				} else {
					t_lower[c_idx] = times[lower_bounds[c_idx]];
					t_upper[c_idx] = times[upper_bounds[c_idx]];

					state_u[c_idx] = state_history[upper_bounds[c_idx]][2*total_states + idx];
					state_l[c_idx] = state_history[lower_bounds[c_idx]][2*total_states + idx];
				}
			}

			immutable Chunk dt = t_upper[] - t_lower[];
			immutable Chunk curr_t = times[get_circular_index(curr_state)];
			immutable Chunk t_delta = curr_t[] - t[];
			immutable Chunk t_delta_abs = abs(t_delta);
			time_delay_beta_buffer[idx][] = 0;
			
			immutable Chunk time_delta = t[] - t_lower[];
			immutable Chunk state_delta = state_u[] - state_l[];
			time_delay_beta_buffer[idx][] = state_l[] + time_delta[]*state_delta[]/dt[];

			foreach(c_idx; 0..chunk_size) {
				if(lower_bounds[c_idx] >= 0) {
					if(t_delta_abs[c_idx] <= dt[c_idx]) {
						time_delay_beta_buffer[idx][c_idx] = state_history[get_circular_index(curr_state)][2*total_states + idx];
					}
				}
			}
		}
		
		foreach(i; 0..total_odd_states) {
			time_delay_a_buffer[i][] = 0;
			time_delay_a_buffer[total_states + i][] = 0;
			foreach(k; 0..total_odd_states) {
				Chunk tmp1 = (time_delay_a_buffer[i][]);
				Chunk tmp2 = A_inv[i][k]*time_delay_alpha_buffer[k][];

				time_delay_a_buffer[i][] = (tmp1[] + tmp2[]);

				tmp1 = (time_delay_a_buffer[total_states + i][]);
				tmp2 = A_inv[i][k]*time_delay_alpha_buffer[total_states + k][];
				time_delay_a_buffer[total_states + i][] = (tmp1[] + tmp2[]);
			}
		}

		foreach(i; 0..total_odd_sin_states) {
			time_delay_b_buffer[i][] = 0;
			time_delay_b_buffer[total_sin_states + i][] = 0;
			foreach(k; 0..total_odd_sin_states) {
				Chunk tmp1 = (time_delay_b_buffer[i][]);
				Chunk tmp2 = A_s_inv[i][k]*time_delay_beta_buffer[k][];

				time_delay_b_buffer[i][] = (tmp1[] + tmp2[]);

				tmp1 = (time_delay_b_buffer[total_sin_states + i][]);
				tmp2 = A_s_inv[i][k]*time_delay_beta_buffer[total_sin_states + k][];
				time_delay_b_buffer[total_sin_states + i][] = (tmp1[] + tmp2[]);
			}
		}

		foreach(i; 0..total_even_states) {
			time_delay_a_buffer[i + total_odd_states][] = time_delay_alpha_buffer[i + total_odd_states][];
			time_delay_a_buffer[i + total_odd_states + total_states][] = time_delay_alpha_buffer[i + total_odd_states + total_states][];
		}

		foreach(i; 0..total_even_sin_states) {
			time_delay_b_buffer[i + total_odd_sin_states][] = time_delay_beta_buffer[i + total_odd_sin_states][];
			time_delay_b_buffer[i + total_odd_sin_states + total_sin_states][] = time_delay_beta_buffer[i + total_odd_sin_states + total_sin_states][];
		}
	}

	@nogc Chunk interpolate_contraction_ratio(immutable Chunk z) {
		Chunk K;

		size_t[chunk_size] lower_bounds;
		size_t[chunk_size] upper_bounds;

		Chunk z_upper;
		Chunk z_lower;
		Chunk K_upper;
		Chunk K_lower;
		foreach(idx, ref _z; z[].map!(a => a.isNaN ? 0 : a).enumerate) {
			lower_bounds[idx] = find_z_bracket(_z);
			if(lower_bounds[idx] == (z_length - 1)) {
				lower_bounds[idx]--;
				upper_bounds[idx] = lower_bounds[idx] + 1;
			} else {
				upper_bounds[idx] = lower_bounds[idx] + 1;
			}

			z_lower[idx] = contraction_z_array_alias[lower_bounds[idx]];
			z_upper[idx] = contraction_z_array_alias[upper_bounds[idx]];

			K_upper[idx] = contraction_array_alias[upper_bounds[idx]];
			K_lower[idx] = contraction_array_alias[lower_bounds[idx]];

			immutable Chunk dz = z_upper[] - z_lower[];
			immutable Chunk z_delta = z[] - z_lower[];
			immutable Chunk K_delta = K_upper[] - K_lower[];
			K[] = (K_lower[] + z_delta[]*K_delta[]/dz[]);
		}

		return K;
	}

	@nogc double wake_skew() {
		return chi;
	}

	@nogc double get_average_inflow() {
		return average_inflow;
	}

	void update(double C_T, ref RotorInputStateT!AC rotor, ref RotorStateT!AC rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		update_impl(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
	}

	void update(double C_T, RotorInputStateT!AC* rotor, RotorStateT!AC* rotor_state, double advance_ratio, double axial_advance_ratio, double dt) {
		update_impl(C_T, rotor, rotor_state, advance_ratio, axial_advance_ratio, dt);
	}

	private void compute_loading(RS)(auto ref RS rotor_state, double radius, double V_inf) {

		tau_c[] = 0;
		tau_s[] = 0;

		immutable omega_sgn = sgn(omega);

		foreach(b_idx, ref blade_state; rotor_state.blade_states) {
			size_t sin_idx = 0;
			auto _idx = iterate_odds!(
				(m, n, idx) {
					foreach(c_idx, ref chunk; rotor.blades[b_idx].chunks) {

						//Chunk atan_num = -omega_sgn*chunk.xi[];
						//immutable Chunk psi_r = atan2(atan_num, chunk.r);
						immutable Chunk mpsi = m.to!double*(blade_state.azimuth /++ psi_r[]+/ /++ PI/12.0+/);

						Chunk cos_mpsi;
						Chunk sin_mpsi;
						if(m == 0) {
							cos_mpsi[] = 1.0/(2.0*PI);
							sin_mpsi[] = 0;
						} else {
							immutable Chunk[2] sin_cos = sincos(mpsi)[];
							cos_mpsi[] = 1.0/(PI)*sin_cos[1][];
							sin_mpsi[] = 1.0/(PI)*sin_cos[0][];
						}

						Chunk nu = 1.0 - chunk.r[]*chunk.r[];
						nu = sqrt(nu);
						
						immutable Chunk Pmn = associated_legendre_polynomial_nh(m, n, nu, P_coefficients_nh[idx]);
						blade_scratch[c_idx][] = blade_state.chunks[c_idx].dC_T[]*Pmn[]*cos_mpsi[];
						if(m != 0) {
							blade_scratch_s[c_idx][] = blade_state.chunks[c_idx].dC_T[]*Pmn[]*sin_mpsi[];
						}
					}

					tau_c[idx] += integrate_trapaziodal(blade_scratch, rotor.blades[b_idx]);
					if(m != 0) {
						tau_s[sin_idx] += integrate_trapaziodal(blade_scratch_s, rotor.blades[b_idx]);
						sin_idx++;
					}

				}
			)(Mo, 0);

			iterate_evens!(
				(m, n, idx) {
					foreach(c_idx, ref chunk; rotor.blades[b_idx].chunks) {
						//Chunk atan_num = -omega_sgn*chunk.xi[];
						//immutable Chunk psi_r = atan2(atan_num, chunk.r);
						immutable Chunk mpsi = m.to!double*(blade_state.azimuth /++ psi_r[]+/ /++ PI/12.0+/);

						Chunk cos_mpsi;
						Chunk sin_mpsi;
						if(m == 0) {
							cos_mpsi[] = 1.0/(2.0*PI);
							sin_mpsi[] = 0;
						} else {
							immutable Chunk[2] sin_cos = sincos(mpsi)[];
							cos_mpsi[] = 1.0/(PI)*sin_cos[1][];
							sin_mpsi[] = 1.0/(PI)*sin_cos[0][];
						}
						
						Chunk nu = 1.0 - chunk.r[]*chunk.r[];
						nu = sqrt(nu);
						
						immutable Chunk Pmn = associated_legendre_polynomial(m, n, nu, P_coefficients[idx]);

						blade_scratch[c_idx][] = blade_state.chunks[c_idx].dC_T[]*Pmn[]*cos_mpsi[];
						if(m != 0) {
							blade_scratch_s[c_idx][] = blade_state.chunks[c_idx].dC_T[]*Pmn[]*sin_mpsi[];
						}
					}
					tau_c[idx] += integrate_trapaziodal(blade_scratch, rotor.blades[b_idx]);

					if(m != 0) {
						tau_s[sin_idx] += integrate_trapaziodal(blade_scratch_s, rotor.blades[b_idx]);
						sin_idx++;
					}

				}
			)(Me, _idx);
		}
		/+
		tau_c[] = 0;
		tau_s[] = 0;

		foreach(b_idx, ref blade_state; rotor_state.blade_states) {

			/+foreach(m; 0..max(Mo, Me) + 1) {
				//immutable Chunk mpsi = m.to!double*coords.psi[];
				immutable mpsi = m.to!double*blade_state.azimuth;
				mpsi_buff[m] = cos(mpsi);
			}+/
			size_t sin_idx = 0;
			auto _idx = iterate_odds!(
				(m, n, idx) {
					//immutable mpsi = m.to!double*(blade_state.azimuth + PI/4.0);
					
					//immutable mpsi = m.to!double*(blade_state.azimuth + PI);
					//immutable Chunk cos_mpsi = m == 0 ? 1.0/(4.0*PI) : 1.0/(2.0*PI)*mpsi_buff[m][0];
					//immutable Chunk cos_mpsi = m == 0 ? 1.0/(4.0*PI) : 1.0/(2.0*PI)*cos(mpsi);
					//immutable Chunk cos_mpsi = m == 0 ? 0.5 : cos(mpsi);
					//immutable Chunk cos_mpsi = m == 0 ? 1.0/(4.0) : 1.0/(2.0)*cos(mpsi);
					Chunk cos_mpsi = m == 0 ? 1.0/(1.0*PI) : 2.0/(PI)*cos(mpsi);
					Chunk sin_mpsi = m == 0 ? 0.0 : 2.0/(PI)*sin(mpsi);
					//Chunk cos_mpsi = m == 0 ? 1.0/(2.0*PI) : 1.0/(PI)*cos(mpsi);

					/+if((m != 0) && (m != 4)) {
						cos_mpsi[] = -cos_mpsi[];
					}+/
					
					foreach(c_idx, ref chunk; rotor.blades[b_idx].chunks) {
					//foreach(c_idx; 0..rotor.blades[b_idx].chunks.length) {
						//BladeGeometryChunk* chunk = rotor.blades[b_idx].chunks[c_idx];
						//immutable Chunk atan_arg = -omega_sgn*chunk.xi[]/chunk.r[];
						Chunk atan_num = -omega_sgn*chunk.xi[];
						immutable Chunk psi_r = atan2(atan_num, chunk.r);
						immutable Chunk mpsi = m.to!double*(blade_state.azimuth + psi_r[]);

						Chunk cos_mpsi;// = m == 0 ? 1.0/(1.0*PI) : 2.0/(PI)*cos(mpsi)[];
						if(m == 0) {
							cos_mpsi[] = 1.0/(2.0*PI);
						} else {
							cos_mpsi[] = 1.0/(PI)*cos(mpsi)[];
						}

						Chunk nu = 1.0 - chunk.r[]*chunk.r[];
						nu = sqrt(nu);
						
						immutable Chunk Pmn = associated_legendre_polynomial_nh(m, n, nu, P_coefficients_nh[idx]);
						blade_scratch[c_idx][] = blade_state.chunks[c_idx].dC_T[]*Pmn[]*cos_mpsi[];
						if(m != 0) {
							blade_scratch_s[c_idx][] = blade_state.chunks[c_idx].dC_T[]*Pmn[]*sin_mpsi[];
						}
					}

					tau_c[idx] += integrate_trapaziodal(blade_scratch, rotor.blades[b_idx]);
					if(m != 0) {
						tau_s[sin_idx] += integrate_trapaziodal(blade_scratch_s, rotor.blades[b_idx]);
						sin_idx++;
					}

				}
			)(Mo, 0);

			iterate_evens!(
				(m, n, idx) {
					//immutable mpsi = m.to!double*(blade_state.azimuth);
					//immutable mpsi = m.to!double*(blade_state.azimuth + PI/4.0);
					//immutable mpsi = m.to!double*(blade_state.azimuth + PI);
					//immutable Chunk cos_mpsi = m == 0 ? 1.0/(4.0*PI) : 1.0/(2.0*PI)*mpsi_buff[m][0];
					//immutable Chunk cos_mpsi = m == 0 ? 0.5 : cos(mpsi);
					//immutable Chunk cos_mpsi = m == 0 ? 1.0/(4.0) : 1.0/(2.0)*cos(mpsi);
					//immutable Chunk cos_mpsi = m == 0 ? 1.0/(4.0*PI) : 1.0/(2.0*PI)*cos(mpsi);
					//Chunk cos_mpsi = m == 0 ? 1.0/(2.0*PI) : 1.0/(PI)*cos(mpsi);
					Chunk cos_mpsi = m == 0 ? 1.0/(1.0*PI) : 2.0/(PI)*cos(mpsi);
					Chunk sin_mpsi = m == 0 ? 0.0 : 2.0/(PI)*sin(mpsi);
					//Chunk cos_mpsi = m == 0 ? 1.0/(1.0*PI) : 2.0/(PI)*cos(mpsi);

					/+if((m != 0) && (m != 4)) {
						cos_mpsi[] = -cos_mpsi[];
					}+/
					foreach(c_idx, ref chunk; rotor.blades[b_idx].chunks) {
					//foreach(c_idx; 0..rotor.blades[b_idx].chunks.length) {
						//BladeGeometryChunk* chunk = rotor.blades[b_idx].chunks[c_idx];
						//immutable Chunk atan_arg = -omega_sgn*chunk.xi[]/chunk.r[];
						//immutable Chunk psi_r = atan(atan_arg);
						Chunk atan_num = -omega_sgn*chunk.xi[];
						immutable Chunk psi_r = atan2(atan_num, chunk.r);
						immutable Chunk mpsi = m.to!double*(blade_state.azimuth + psi_r[]);

						Chunk cos_mpsi;// = m == 0 ? 1.0/(1.0*PI) : 2.0/(PI)*cos(mpsi)[];
						if(m == 0) {
							cos_mpsi[] = 1.0/(2.0*PI);
						} else {
							cos_mpsi[] = 1.0/(PI)*cos(mpsi)[];
						}
						
						Chunk nu = 1.0 - chunk.r[]*chunk.r[];
						nu = sqrt(nu);
						
						immutable Chunk Pmn = associated_legendre_polynomial(m, n, nu, P_coefficients[idx]);

						blade_scratch[c_idx][] = blade_state.chunks[c_idx].dC_T[]*Pmn[]*cos_mpsi[];
						if(m != 0) {
							blade_scratch_s[c_idx][] = blade_state.chunks[c_idx].dC_T[]*Pmn[]*sin_mpsi[];
						}
					}
					tau_c[idx] += integrate_trapaziodal(blade_scratch, rotor.blades[b_idx]);
					if(m != 0) {
						tau_s[sin_idx] += integrate_trapaziodal(blade_scratch_s, rotor.blades[b_idx]);
						sin_idx++;
					}
				}
			)(Me, _idx);
		}
		+/
	}

	// size_t n_r = 16;
	// size_t n_psi = 6;
	// double v_0;

	// Chunk[] contraction_array;
	// double[] contraction_array_alias;
	// Chunk[] contraction_z_array;
	// double[] contraction_z_array_alias;

	private void update_impl(RIS, RS)(double C_T, auto ref RIS rotor_input, auto ref RS rotor_state, double _advance_ratio, double _axial_advance_ratio, double dt) {

		omega = rotor_input.angular_velocity;
		immutable V_inf = rotor_input.freestream_velocity;
		immutable t_scale = abs(omega);//V_inf/rotor.radius;

		advance_ratio = _advance_ratio;//*abs(omega)*rotor.radius/V_inf;
		axial_advance_ratio = _axial_advance_ratio;//*abs(omega)*rotor.radius/V_inf;
		
		auto time = dt*curr_state.to!double*t_scale;//abs(omega);

		times[get_circular_index(curr_state + 1)] = time;

		compute_loading(rotor_state, rotor.radius, V_inf);
		
		integrator.step!(system_derivative)(state_history[get_circular_index(curr_state + 1)], state_history[get_circular_index(curr_state)], time, dt*t_scale, this, rotor_input, rotor_state, advance_ratio, axial_advance_ratio);

		alpha = state_history[get_circular_index(curr_state + 1)][0..total_states];
		beta = state_history[get_circular_index(curr_state + 1)][2*total_states..2*total_states + total_sin_states];

		//b[] = 0;
		//beta[] = 0;
		// Get MD state variables
		cblas_dgemv(CblasRowMajor, CblasNoTrans, total_odd_states, total_odd_states, 1.0, A_inv[0].ptr, total_odd_states, alpha.ptr, 1, 0.0, a.ptr, 1);

		cblas_dgemv(CblasRowMajor, CblasNoTrans, total_odd_sin_states, total_odd_sin_states, 1.0, A_s_inv[0].ptr, total_odd_sin_states, beta.ptr, 1, 0.0, b.ptr, 1);

		a[total_odd_states..$] = alpha[total_odd_states..$];
		b[total_odd_sin_states..$] = beta[total_odd_sin_states..$];

		LM_c_scratch[] *= alpha[];

		immutable alpha_0_1 = LM_c_scratch.sum;

		average_inflow = sqrt(3.0)*alpha_0_1;

		foreach(i; 0..total_states) {
			state_history[get_circular_index(curr_state + 1)][total_states + i] = 0;
			foreach(k; 0..total_states) {
				state_history[get_circular_index(curr_state + 1)][total_states + i] += QS_c_mat[i][k]*alpha[k];
			}
		}
/+
		size_t sin_idx = 0;
		auto _idx = iterate_odds!(
			(m, n, idx) {
				if(m != 0) {
					state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states + sin_idx] = (1.0)^^(n.to!double + 0.0)*state_history[get_circular_index(curr_state + 1)][total_states + idx];
					//delta_s[sin_idx] = (1.0)^^(n.to!double + 0.0)*delta[idx];
					sin_idx++;
				}
			}
		)(Mo, 0);

		iterate_evens!(
			(m, n, idx) {
				if(m != 0) {
					state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states + sin_idx] = (1.0)^^(n.to!double + 0.0)*state_history[get_circular_index(curr_state + 1)][total_states + idx];
					//lambda_s[sin_idx] = (1.0)^^(n.to!double + 0.0)*lambda[idx];
					//delta_s[sin_idx] = (1.0)^^(n.to!double + 0.0)*delta[idx];
					sin_idx++;
				}
			}
		)(Me, _idx);
+/
		/+size_t sin_idx = 0;
		auto _idx = iterate_odds!(
			(m, n, idx) {
				if(m != 0) {
					state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states + sin_idx] = state_history[get_circular_index(curr_state + 1)][total_states + idx];
					sin_idx++;
				}
			}
		)(Mo, 0);

		iterate_evens!(
			(m, n, idx) {
				if(m != 0) {
					state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states + sin_idx] = state_history[get_circular_index(curr_state + 1)][total_states + idx];
					sin_idx++;
				}
			}
		)(Me, _idx);+/

		foreach(i; 0..total_sin_states) {
			state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states + i] = 0;
			foreach(k; 0..total_sin_states) {
				state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states + i] += QS_s_mat[i][k]*beta[k];
			}
		}

		lambda = state_history[get_circular_index(curr_state + 1)][total_states..2*total_states];
		lambda_s = state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states..$];

		// Get MD adjoint variables
		cblas_dgemv(CblasRowMajor, CblasNoTrans, total_odd_states, total_odd_states, 1.0, A_inv[0].ptr, total_odd_states, lambda.ptr, 1, 0.0, delta.ptr, 1);

		cblas_dgemv(CblasRowMajor, CblasNoTrans, total_odd_sin_states, total_odd_sin_states, 1.0, A_s_inv[0].ptr, total_odd_sin_states, lambda_s.ptr, 1, 0.0, delta_s.ptr, 1);
		
		delta[total_odd_states..$] = lambda[total_odd_states..$];
		delta_s[total_odd_sin_states..$] = lambda_s[total_odd_sin_states..$];
		//delta_s[] = 0;
		//lambda_s[] = 0;
		/+size_t sin_idx = 0;
		auto _idx = iterate_odds!(
			(m, n, idx) {
				if(m != 0) {
					lambda_s[sin_idx] = (1.0)^^(n.to!double + 1.0)*lambda[idx];
					delta_s[sin_idx] = (1.0)^^(n.to!double + 1.0)*delta[idx];
					sin_idx++;
				}
			}
		)(Mo, 0);

		iterate_evens!(
			(m, n, idx) {
				if(m != 0) {
					lambda_s[sin_idx] = (1.0)^^(n.to!double + 1.0)*lambda[idx];
					delta_s[sin_idx] = (1.0)^^(n.to!double + 1.0)*delta[idx];
					sin_idx++;
				}
			}
		)(Me, _idx);+/

		curr_state++;

		//delta_s[] = 0;
		//lambda_s[] = 0;
		//b[] = 0;
		//beta[] = 0;
		//state_history[get_circular_index(curr_state + 1)][2*total_states + total_sin_states..$] = 0;
		v_0 = compute_inflow_average_at_disk;

		auto v_inf = sqrt(advance_ratio*advance_ratio + axial_advance_ratio*axial_advance_ratio);

		foreach(z_idx, ref z_chunk; contraction_z_array) {
			auto v_z = compute_inflow_average_at_z(z_chunk);
			immutable Chunk v_0_v_z = (v_inf + v_0)/(v_inf + v_z[]);
			contraction_array[z_idx] = sqrt(v_0_v_z);
		}
	}

	double compute_inflow_average_at_disk() {
		Chunk v_z = 0;

		double d_psi = 2.0*PI/n_psi.to!double;
		double d_r = 1.0/n_r.to!double;
		immutable Chunk z = 0;

		foreach(psi_i; 0..n_psi) {
			immutable double psi = 2.0*PI*psi_i.to!double/n_psi.to!double;
			immutable double cos_psi = cos(psi);
			immutable double sin_psi = sin(psi);
			foreach(ref r_idx_chunk; iota(0, n_r, 1).chunks(chunk_size)) {
				immutable Chunk r_idx_d = r_idx_chunk.map!(a => a.to!double).staticArray!Chunk;
				immutable Chunk r_chunk = r_idx_d[]/n_r.to!double;
				immutable Chunk x_r = r_chunk[]*cos_psi;
				immutable Chunk y_r = r_chunk[]*sin_psi;
				immutable v = inflow_at_impl(this, x_r, y_r, z);
				v_z[] += v[];
			}
		}

		return v_z.sum*d_r*d_psi*1.0/PI;
	}

	Chunk compute_inflow_average_at_z(immutable Chunk z) {
		Chunk v_z = 0;

		immutable Chunk x_c = -z[]*tan_chi;

		double d_psi = 2.0*PI/n_psi.to!double;
		double d_r = 1.0/n_r.to!double;

		foreach(psi_i; 0..n_psi) {
			immutable double psi = 2.0*PI*psi_i.to!double/n_psi.to!double;
			immutable double cos_psi = cos(psi);
			immutable double sin_psi = sin(psi);
			foreach(r_j; 0..n_r) {
				immutable Chunk x_r = (r_j.to!double/n_r.to!double)*cos_psi;
				immutable Chunk y_r = (r_j.to!double/n_r.to!double)*sin_psi;
				immutable Chunk x = x_c[] + x_r[];
				immutable v = inflow_at_impl(this, x, y_r, z);
				v_z[] += v[];
			}
		}

		v_z[] *= d_r*d_psi*1.0/PI;
		return v_z;
	}

	@nogc Chunk compute_contraction_multiplier(immutable Chunk x, immutable Chunk y, immutable Chunk z) {
		immutable K = interpolate_contraction_ratio(z);

		immutable Chunk x_0 = z[]*tan_chi;
		immutable Chunk x_r = x_0[] + x[];

		immutable Chunk r_0_squared = x_r[]*x_r[] + y[]*y[];
		immutable Chunk r_0 = sqrt(r_0_squared);
		
		Chunk k_bar = 0;
		foreach(r_idx, ref _r_0; r_0[]) {
			if(_r_0 <= 1.0) {
				k_bar[r_idx] = K[r_idx];
			} else {
				k_bar[r_idx] = (K[r_idx] + _r_0 - 1.0)/_r_0;
			}
		}

		return k_bar;
	}

	Chunk inflow_at(immutable Chunk x, immutable Chunk y, immutable Chunk z, immutable Chunk x_e, double angle_of_attack) {
		// Rotate into the coord frame the huang model expects.
		if(!contraction_mapping) {
			//immutable Chunk neg_y = -y[];
			//immutable V = omega > 0 ? inflow_at_impl(this, x, neg_y, z) : inflow_at_impl(this, x, y, z);
			immutable V = inflow_at_impl(this, x, y, z);
			return V;
		} else {

			immutable k_bar = compute_contraction_multiplier(x, y, z);			

			immutable Chunk x_c = x[]/k_bar[];
			immutable Chunk y_c = y[]/k_bar[];
			immutable Chunk neg_y = -y_c[];
			//immutable V = omega > 0 ? inflow_at_impl(this, x_c, neg_y, z) : inflow_at_impl(this, x_c, y_c, z);
			immutable V = inflow_at_impl(this, x_c, y_c, z);
			return V;
		}
	}

	@nogc Chunk inflow_at(immutable Chunk r, immutable double cos_azimuth, immutable double sin_azimuth) {
		Chunk i;
		//immutable V = inflow_at_impl(this, x, y, z);
		return i;
	}
}
