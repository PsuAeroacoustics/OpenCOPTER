module validation.hart_ii.pressure;

import plt = matplotlibd.pyplot;

import std.algorithm;
import std.array;
import std.conv;
import std.getopt;
import std.file;
import std.path;
import std.range;
import std.stdio;
import std.string;

struct Dataset {
	Channel[] channels;
}

struct Function {
	float[] data;
}

struct Channel {
	
	float[] time;
	Function[] functions;
}

auto read_dataset(string dir) {
	int[] header_buffer = new int[4];

	auto pressure_file = File(dir.buildPath("pressure.fn"), "rb");
	scope(exit) pressure_file.close;

	Dataset dataset;

	header_buffer = pressure_file.rawRead(header_buffer);

	int i_max = header_buffer[0];
	int j_max = header_buffer[1];
	int samples = header_buffer[2];
	int functions = header_buffer[3];

	size_t total_obs = i_max*j_max;

	dataset.channels = new Channel[total_obs];

	foreach(ref channel; dataset.channels) {
		channel.time = new float[samples];
		channel.functions = new Function[functions - 1];
		foreach(ref func; channel.functions) {
			func.data = new float[samples];
		}
	}

	float[] sample_buff = new float[1];

	foreach(s_idx; 0..samples) {
		foreach(ref channel; dataset.channels) {
			sample_buff = pressure_file.rawRead(sample_buff);
			channel.time[s_idx] = sample_buff[0];
		}
	}

	foreach(f_idx; 0..(functions - 1)) {
		foreach(s_idx; 0..samples) {
			foreach(ref channel; dataset.channels) {
				sample_buff = pressure_file.rawRead(sample_buff);
				channel.functions[f_idx].data[s_idx] = sample_buff[0];
			}
		}
	}

	return dataset;
}

void plot_mic(T, M)(auto ref T time, ref Channel channel, auto ref M measured, bool plot_thickness, bool plot_loading, bool plot_total, string title) {

	string[] legend;

	plt.figure;
	if(plot_thickness) {
		plt.plot(time, channel.functions[0].data);
		legend ~= "Thickness";
	}
	if(plot_loading) {
		plt.plot(time, channel.functions[1].data);
		legend ~= "Loading";
	}
	if(plot_total) {
		plt.plot(time, channel.functions[2].data);
		legend ~= "Total";
	}
	//plt.plot(measured[0], measured[1]);
	legend ~= "Measured";
	plt.title(title);

	plt.legend(legend);
}

int main(string[] args) {
	
	string dir = "./";

	bool plot_loading = false;
	bool plot_thickness = false;
	bool plot_total = true;

	arraySep = ",";
	auto help_information = getopt(
		args,
		std.getopt.config.bundling,
		"directory|d", "Directory where pressure data is. Default = "~dir, &dir,
		"loading|l", "Plot loading data", &plot_loading,
		"thickness|k", "Plot thickness data", &plot_thickness,
		"total|t", "Plot total noise data", &plot_total
	);

	if(help_information.helpWanted) {
		defaultGetoptPrinter("Program options:", help_information.options);
		import core.stdc.stdlib : exit;
		exit(-1);
	}

	auto dataset = read_dataset(dir);

	auto m2_txt = readText("3014_M2.txt").split("\n")[1..$];
	auto m3_txt = readText("3014_M3.txt").split("\n")[1..$];
	auto m4_txt = readText("3014_M4.txt").split("\n")[1..$];

	//writeln(m2_txt);

	auto m2_measured = m2_txt.map!(line => line.split("\t").map!(a => a.to!double).array).filter!(a => a.length > 0).array.transposed;
	auto m3_measured = m3_txt.map!(line => line.split("\t").map!(a => a.to!double).array).filter!(a => a.length > 0).array.transposed;
	auto m4_measured = m4_txt.map!(line => line.split("\t").map!(a => a.to!double).array).filter!(a => a.length > 0).array.transposed;

	//writeln(m2_measured[0]);

	
	//m2_measured[0][] += dataset.channels[0].time[0]/0.0267;
	//plt.figure;
	//plt.plot(m2_measured[0], m2_measured[1]);
	//plt.title("M2");

	auto time = dataset.channels[0].time.map!(x => (x - dataset.channels[0].time[0])/0.0267);

	//writeln(dataset.channels[0].functions[1].data);
	plot_mic(time, dataset.channels[0], m2_measured, plot_thickness, plot_loading, plot_total, "M2");
	plot_mic(time, dataset.channels[1], m3_measured, plot_thickness, plot_loading, plot_total, "M3");
	plot_mic(time, dataset.channels[2], m4_measured, plot_thickness, plot_loading, plot_total, "M4");
	/+plt.figure;
	//plt.plot(time, dataset.channels[0].functions[0].data, time, dataset.channels[0].functions[1].data, time, dataset.channels[0].functions[2].data, m2_measured[0], m2_measured[1]);
	if(plot_thickness) {
		plt.plot(time, dataset.channels[0].functions[0].data);
	}
	if(plot_loading) {
		plt.plot(time, dataset.channels[0].functions[1].data);
	}
	if(plot_total) {
		plt.plot(time, dataset.channels[0].functions[2].data);
	}
	plt.plot(m2_measured[0], m2_measured[1]);

	plt.title("M2");
	plt.legend(["Thickness", "Loading", "Total", "Measured"]);
	
	plt.figure;
	plt.plot(time, dataset.channels[1].functions[0].data, time, dataset.channels[1].functions[1].data, time, dataset.channels[1].functions[2].data, m3_measured[0], m3_measured[1]);
	plt.title("M3");
	plt.legend(["Thickness", "Loading", "Total", "Measured"]);

	plt.figure;
	plt.plot(time, dataset.channels[2].functions[0].data, time, dataset.channels[2].functions[1].data, time, dataset.channels[2].functions[2].data, m4_measured[0], m4_measured[1]);
	plt.title("M4");
	plt.legend(["Thickness", "Loading", "Total", "Measured"]);

	plt.show;+/
	plt.show;
	return 0;
}
