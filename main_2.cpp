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

    for (int number = 0; number < 250; ++number) {
        double N_n = return_N_plus_1(0.2, h, 0.2, alpha, beta, delta, number);
        double P_n = return_P_plus_1(0.2, h, 0.2, alpha, beta, delta, number);
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
    surface = cairo_pdf_surface_create("result.pdf", N_size_x * coef, P_size_y * coef); // Создаем поверхность для рисования в виде файла PDF
    cout << N_size_x * coef << ' ' << P_size_y * coef << '\n';
    for (int number = 0; number < 250; ++number) {
        double N_n = return_N_plus_1(0.2, h, 0.2, alpha, beta, delta, number);
        double P_n = return_P_plus_1(0.2, h, 0.2, alpha, beta, delta, number);
        /*
        cout << "Now step: " << number << " - result N_step = " << N_n
             << "\t - result P_step = " << P_n << '\n';
        */
        cr = cairo_create(surface);
        cairo_set_source_rgb(cr, 112.0/255, 128.0/255, 144.0/255);
        cairo_set_line_width(cr, 0.001);
        cairo_arc(cr,
                  ((N_n - N_min) / (N_max - N_min)) * N_size_x * coef,
                  ((P_n - P_min) / (P_max - P_min)) * P_size_y * coef, 0.25, 0, 2*3.14159265359);
        cairo_fill(cr);
        cairo_destroy(cr);
    }

    cairo_surface_destroy(surface);
    return 0;
}
