

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include "Fourier.h"
#include <vector>
#include <fstream>


int main()
{

    int N = 100;
    std::vector<double> F1(N), F2(N), F3(N);

    for (int i = 0; i < N; i++) {
        F1[i] = sin(2 * M_PI * i / N * 20);
        F2[i] = sin(2 * M_PI * i / N * 60);
        F3[i] = sin(2 * M_PI * i / N * 120);
    }

    std::vector<std::pair<double, double>> spectrum_F1 = Fourier(F1),
        spectrum_F2 = Fourier(F2),
        spectrum_F3 = Fourier(F3);
    std::cout << "spectrum F1" << std::endl;
    print_spectrum(spectrum_F1);
    spectrum_to_file(spectrum_F1, "spectrum_F1.csv");
    signal_to_file(F1, reFourier(spectrum_F1), "F1.csv");

    std::cout << "spectrum F2" << std::endl;
    print_spectrum(spectrum_F2);
    spectrum_to_file(spectrum_F2, "spectrum_F2.csv");
    signal_to_file(F2, reFourier(spectrum_F2), "F2.csv");

    std::cout << "spectrum F3" << std::endl;
    print_spectrum(spectrum_F3);
    spectrum_to_file(spectrum_F3, "spectrum_F3.csv");
    signal_to_file(F3, reFourier(spectrum_F3), "F3.csv");

    std::vector<double> F_SIN(N), F_COS(N);

    for (int i = 0; i < N; i++) {
        F_SIN[i] = sin(2 * M_PI * i / N);
        F_COS[i] = cos(2 * M_PI * i / N);
    }

    std::vector<std::pair<double, double>> spectrum_F_SIN = Fourier(F_SIN),
        spectrum_F_COS = Fourier(F_COS);

    std::cout << "spectrum F_SIN" << std::endl;
    print_spectrum(spectrum_F_SIN);
    spectrum_to_file(spectrum_F_SIN, "spectrum_F_SIN.csv");
    signal_to_file(F_SIN, reFourier(spectrum_F_SIN), "F_SIN.csv");

    std::cout << "spectrum F_COS" << std::endl;
    print_spectrum(spectrum_F_COS);
    spectrum_to_file(spectrum_F_COS, "spectrum_F_COS.csv");
    signal_to_file(F_COS, reFourier(spectrum_F_COS), "F_COS.csv");

    std::vector<double> F_SIN_HALF(N), F_COS_HALF(N);

    for (int i = 0; i < N; i++) {
        F_SIN_HALF[i] = sin(2 * M_PI * i / N / 2);
        F_COS_HALF[i] = cos(2 * M_PI * i / N / 2);

    }

    std::vector<std::pair<double, double>> spectrum_F_SIN_HALF = Fourier(F_SIN_HALF),
        spectrum_F_COS_HALF = Fourier(F_COS_HALF);

    std::cout << "spectrum F_SIN_HALF" << std::endl;
    print_spectrum(spectrum_F_SIN_HALF);
    spectrum_to_file(spectrum_F_SIN_HALF, "spectrum_F_SIN_HALF.csv");
    signal_to_file(F_SIN_HALF, reFourier(spectrum_F_SIN_HALF), "F_SIN_HALF.csv");


    std::cout << "spectrum F_COS_HALF" << std::endl;
    print_spectrum(spectrum_F_COS_HALF);
    spectrum_to_file(spectrum_F_COS_HALF, "spectrum_F_COS_HALF.csv");
    signal_to_file(F_COS_HALF, reFourier(spectrum_F_COS_HALF), "F_COS_HALF.csv");

    std::vector<double> random_signal(N);
    for (int i = 0; i < random_signal.size(); i++)
        random_signal[i] = rand() % 200 - 100;

    std::vector<std::pair<double, double>> spectrum_random = Fourier(random_signal);

    std::cout << "spectrum RANDOM" << std::endl;
    print_spectrum(spectrum_random);
    spectrum_to_file(spectrum_random, "spectrum_random.csv");
    signal_to_file(random_signal, reFourier(spectrum_random), "RANDOM.csv");


    std::vector<double> system_a = { 1, -1.1, 0.6 },
        system_b = { 1.5,0 ,0 };

    std::vector<double> input_50(50, 0),
        input_100(100, 0),
        input_200(1000, 0);

    input_50[0] = 1; input_100[0] = 1; input_200[0] = 1;

    std::vector<double> impulse_50 = impuls(system_a, system_b, input_50),
        impulse_100 = impuls(system_a, system_b, input_100),
        impulse_200 = impuls(system_a, system_b, input_200);

    std::vector<std::pair<double, double>> spectrum_50 = Fourier(impulse_50),
        spectrum_100 = Fourier(impulse_100),
        spectrum_200 = Fourier(impulse_200);

    spectrum_to_file(spectrum_50, "spectrum_50.csv");
    spectrum_to_file(spectrum_100, "spectrum_100.csv");
    spectrum_to_file(spectrum_200, "spectrum_200.csv");
    signal_to_file(impulse_200, reFourier(spectrum_200), "impulse_200.csv");


    std::vector<std::pair<double, double>> temp(50);
    for (int i = 0; i < spectrum_random.size(); i++) {
        temp[i].first = spectrum_random[i].first * spectrum_100[i].first - spectrum_random[i].second * spectrum_100[i].second;
        temp[i].second = spectrum_random[i].first * spectrum_100[i].second + spectrum_random[i].second * spectrum_100[i].first;
    }

    std::vector<double> answer_1 = impuls(system_a, system_b, random_signal),
        answer_2 = reFourier(temp),
        answer_3 = convolution(random_signal, impulse_100);



    signal_to_file(answer_1, answer_2, "diff.csv");

    std::vector<double> test(N);
    for (int i = 0; i < N; i++)
        test[i] = sin(M_PI * 2 * i / N);
    auto answer = aboba(Fourier(test));

    signal_to_file(answer, answer, "aboba.csv");

    long signal_size = 32000;
    std::vector<double> signal(signal_size),
        filtred(signal_size);

    std::ofstream a("in.csv", std::ios_base::out);
    std::ofstream b("out.csv", std::ios_base::out);

    a << "x,y" << std::endl;
    b << "x,y" << std::endl;
    for (int i = 0; i < signal_size; i++)
        signal[i] = sin(M_PI * 2 * 14000 * i / 32000) + sin(M_PI * 2 * 4000 * i / 32000);

    //sin(M_PI * 2 * 4000 * i / 32000) +
    //signal[i] += sin(M_PI * 2 * 14000 * i / 32000);
   // for (int i = signal_size / 2; i < signal_size; i++)
     //   signal[i] = 0;



    __int64 z1[3] = { 0,0,0 },
        z2[3] = { 0,0,0 },
        z3[3] = { 0,0,0 },
        z4[3] = { 0,0,0 },
        b0[3] = { 14626,19524,27008 },
        b1[3] = { 27560,27217,30767 },
        b2[3] = { 14626,19524,27008 },
        a0[3] = { -21464,-13285,-24554 },
        a1[3] = { -12741,-21937,-30209 },
        out[4] = { 0,0,0,0 },
        temp1 = 0,
        shift[3] = {0,1,1};

    for (int i = 0; i < signal_size; i++) {
        out[0] = (int)((signal[i]) * (1 << 16));
        for (int j = 0; j < 3; j++) {
            temp1 = (out[j] + z1[j]);
            out[j + 1] = ((temp1 * b0[j]) >> (16 - shift[j])) + z4[j];
            z1[j] = ((temp1 * a0[j]) >> (16 - shift[j])) + z2[j];
            z4[j] = ((temp1 * b1[j]) >> (16 - shift[j])) + z3[j];

            z2[j] = ((temp1 * a1[j]) >> (16 - shift[j]));
            z3[j] = ((temp1 * b2[j]) >> (16 - shift[j]));
        }
        filtred[i] = ((double)(out[3])) / (1 << 16);
    }

    for (int i = 0; i < signal_size; i++) {
        a << i << "," << signal[i] << std::endl;
        b << i << "," << filtred[i] << std::endl;
    }


    spectrum_to_file(Fourier(signal), "spec.csv");
    spectrum_to_file(Fourier(filtred), "spec_f.csv");


}


