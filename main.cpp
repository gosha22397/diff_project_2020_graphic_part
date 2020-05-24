#include <iostream>
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>

using namespace std;

double return_N_plus_1 (double N_n, double h, double P_n,
                        double alpha, double beta, double delta, int number) {
    if (number > 0) {
        double N_n_1 = N_n + h * N_n * (1 - N_n - P_n / (N_n + alpha * P_n));
        double P_n_1 = P_n + h * beta * P_n * (delta - P_n / N_n);
        return return_N_plus_1(N_n_1, h, P_n_1, alpha, beta, delta, number - 1);
    } else {
        return N_n;
    }
}

double return_P_plus_1 (double N_n, double h, double P_n,
                        double alpha, double beta, double delta, int number) {
    if (number > 0) {
        double N_n_1 = N_n + h * N_n * (1 - N_n - P_n / (N_n + alpha * P_n));
        double P_n_1 = P_n + h * beta * P_n * (delta - P_n / N_n);
        return return_P_plus_1(N_n_1, h, P_n_1, alpha, beta, delta, number - 1);
    } else {
        return P_n;
    }
}

double return_n (double t) {
    return 0.2 +
           t * (0.042352941176470586) +
           t * t * (0.000161204966415631996743334011805414207205373498880521066);// +
           //t * t * t * (-0.01512000109907781793049112222803657927054156392744715371);
}

double return_p (double t) {
    return 0.2 +
           t * (-0.892) +
           t * t * (2.095042352941176470588235294117647058823529411764705882352);// +
           //t * t * t * (-3.16814847140240179116629350702218603704457561571341339303);
}

int main() {
    cairo_surface_t *surface; // Определяем поверхность для рисования
    cairo_t *cr;              // Определяем источник

    double alpha = 0.7;
    double beta = 0.9;
    double delta = 0.6;
    double h = 0.1;

    double P_min = 0.2;
    double P_max = 0.2;
    double N_min = 0.2;
    double N_max = 0.2;

    for (double number = 0; number < 5; number += 0.0001) {
        double N_n = return_n(number);
        double P_n = return_p(number);
        //double N_n = return_N_plus_1(0.2, h, 0.2, alpha, beta, delta, number);
        //double P_n = return_P_plus_1(0.2, h, 0.2, alpha, beta, delta, number);
        if (P_n < P_min) {
            P_min = P_n;
        }
        if (N_n < N_min) {
            N_min = N_n;
        }
        if (P_n > P_max) {
            P_max = P_n;
        }
        if (N_n > N_max) {
            N_max = N_n;
        }
    }

    double N_size_x = N_max - N_min;
    double P_size_y = P_max - P_min;
    double coef = 100 / min(N_size_x, P_size_y);
    surface = cairo_pdf_surface_create("result.pdf", N_size_x, P_size_y); // Создаем поверхность для рисования в виде файла PDF
    /*
    double N_coef = N_size_x * coef;
    double P_coef = P_size_y * coef;
    if (N_coef == 100) {
        P_coef = 100 / P_coef;
    } else {
        N_coef = 100 / N_coef;
    }
    cout << N_size_x * coef << ' ' << P_size_y * coef << '\n';
    */
    for (double number = 0; number < 5; number += 0.0001) {
        double N_n = return_n(number);
        double P_n = return_p(number);
        //double N_n = return_N_plus_1(0.2, h, 0.2, alpha, beta, delta, number);
        //double P_n = return_P_plus_1(0.2, h, 0.2, alpha, beta, delta, number);
        /*
        cout << "Now step: " << number << " - result N_step = " << N_n
             << "\t - result P_step = " << P_n << '\n';
        */
        cr = cairo_create(surface);
        cairo_set_source_rgb(cr, 112.0/255, 128.0/255, 144.0/255);
        cairo_set_line_width(cr, 0.001);
        cairo_arc(cr,
                  ((N_n - N_min) / (N_max - N_min)) * N_size_x,
                  ((P_n - P_min) / (P_max - P_min)) * P_size_y, 0.2, 0, 2*3.14159265359);
        cairo_fill(cr);
        cairo_destroy(cr);
    }

    cairo_surface_destroy(surface);
    return 0;
}
