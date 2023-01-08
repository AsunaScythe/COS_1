#pragma once


#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>

std::vector<std::pair<double,double>> Fourier(const std::vector<double>& signal);

std::vector<double> reFourier(const std::vector<std::pair<double, double>>& spectrum);

void print_spectrum(const std::vector<std::pair<double, double>>& spectrum);

void signal_to_file(std::vector<double> signal, std::vector<double> reSignal, std::string name);

void spectrum_to_file(const std::vector<std::pair<double, double>> spectrum, std::string name);

std::vector<double> impuls(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& input);

std::vector<double> bandpass_filter(const std::vector<double>& signal, int frequency, double min_freq, double max_freq);

std::vector<double> convolution(const std::vector<double>& x, const std::vector<double>& h);

std::vector<double> aboba(std::vector<std::pair<double, double>> spec);