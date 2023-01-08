#include "Fourier.h"
#include <fstream>

std::vector<std::pair<double,double>> Fourier(const std::vector<double>& signal) {
    std::vector<std::pair<double, double>> spectrum(signal.size() / 2,{0,0});
    for (int k = 0; k < spectrum.size(); k++) 
        for (int i = 0; i < signal.size(); i++) {
            spectrum[k].first += signal[i] * cos(2 * M_PI / signal.size() * i * k);
            spectrum[k].second += signal[i] * sin(2 * M_PI / signal.size() * i * k);
        }   
    return spectrum;
}

 std::vector<double> reFourier(const std::vector<std::pair<double, double>> & spectrum) {
    std::vector<double> signal(spectrum.size()*2,0);
    for (int i = 1; i < spectrum.size(); i++) 
        for (int j = 0; j < signal.size(); j++) 
            signal[j] += spectrum[i].first * cos(M_PI * 2 / signal.size() * i * j) + spectrum[i].second * sin(M_PI * 2 / signal.size() * i * j);

    for (int j = 0; j < signal.size(); j++)
        signal[j] += spectrum[0].first/2;
    for (int j = 0; j < signal.size(); j++)
        signal[j] /= spectrum.size();
    return signal;
}

 void print_spectrum(const std::vector<std::pair<double, double>>& spectrum) {
     for (int i = 0; i < spectrum.size(); i++)
         if ((abs(spectrum[i].second) + abs(spectrum[i].first)) > 0.00000000001) {
             std::cout << i << " sin = " << spectrum[i].second << " cos = " << spectrum[i].first << std::endl;
         }
 }

 void signal_to_file(std::vector<double> signal, std::vector<double> reSignal, std::string name) {
     std::ofstream file(name);
     file << "x,y,z" << std::endl;
     for (int i = 0; i < signal.size(); i++)
         file << i << ',' << signal[i] << ',' << reSignal[i] << std::endl;
     
 }

 void spectrum_to_file(const std::vector<std::pair<double, double>> spectrum, std::string name) {
     std::ofstream file(name);
     file << "x,y" << std::endl;
     for (int i = 0; i < spectrum.size(); i++)
         file << i << ','<< spectrum[i].first << ',' << spectrum[i].second << std::endl;
     
 }

std::vector<double> impuls(const std::vector<double> & a,const std::vector<double> & b,const std::vector<double> & input) {
    std::vector<double> answer(input.size(),0);
     for (int i = 0; i < answer.size(); i++)
         for (int j = 0; j < a.size(); j++) 
             if (i - j >= 0)
                 answer[i] += -answer[i - j] * a[j] + input[i - j] * b[j];
     return answer;
 }

std::vector<double> convolution(const std::vector<double>& x, const std::vector<double>& h) {
    std::vector<double> answer(x.size(),0);
    for (int i = 0; i < answer.size(); i++)
        for (int j = 0; j < i + 1; j++)
            answer[i] += x[j] * h[i-j];
    return answer;
}

/*std::vector<double> bandpass_filter(const std::vector<double>& signal, int frequency, double min_freq, double max_freq) {
    unsigned N = signal.size();
    double time = double(N) / frequency;
    std::vector<double> filtered_signal(N, 0);
    for (int i = int(min_freq * time); i < max_freq * time; i++) {
        std::vector<double> temp = Fourier(signal, i);
        for (int j = 0; j < N; j++) {
            filtered_signal[j] += temp[0] * cos(M_PI * 2 / N * i * j) + temp[1] * sin(M_PI * 2 / N * i * j);
        }
    }
    return filtered_signal;
}
*/

std::vector<double> aboba(std::vector<std::pair<double, double>> spec) {
    std::vector<double> answer(spec.size(), 0);
    for (int i = 0; i < spec.size(); i++)
        answer[i] = atan(spec[i].first / spec[i].second);
    return answer;



}

