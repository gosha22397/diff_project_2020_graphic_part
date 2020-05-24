#include <math.h>
#include <iostream>
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>

using namespace std;

double N_k_NSFD_hyt (double N_k, double h, double P_k,
                     double a, double b, double y, int number) {
    if (number > 0) {
        double N_k_1 = N_k * (1 + h + h * (N_k + a * P_k)) * (N_k + a * P_k) / ((1 + 2 * h * N_k + a * h * P_k) * (N_k + a * P_k) + h * P_k);
        double P_k_1 = P_k * N_k * (1 + b * y * h) / (N_k + b * h * P_k);
        return N_k_NSFD_hyt(N_k_1, h, P_k_1, a, b, y, number - 1);
    } else {
        return N_k;
    }
}

double P_k_NSFD_hyt (double N_k, double h, double P_k,
                     double a, double b, double y, int number) {
    if (number > 0) {
        double N_k_1 = N_k * (1 + h + h * (N_k + a * P_k)) * (N_k + a * P_k) / ((1 + 2 * h * N_k + a * h * P_k) * (N_k + a * P_k) + h * P_k);
        double P_k_1 = P_k * N_k * (1 + b * y * h) / (N_k + b * h * P_k);
        return P_k_NSFD_hyt(N_k_1, h, P_k_1, a, b, y, number - 1);
    } else {
        return P_k;
    }
}

double N_k_NSFD_only (double N_k, double h, double P_k,
                      double a, double b, double y, int number) {
    if (number > 0) {
        double N_k_1 = (h * N_k + N_k) * (N_k + a * P_k) / ((1 + h * N_k) * (N_k + a * P_k) + h * P_k);
        double P_k_1 = P_k * N_k * (1 + b * y * h) / (N_k + b * h * P_k);
        return N_k_NSFD_only(N_k_1, h, P_k_1, a, b, y, number - 1);
    } else {
        return N_k;
    }
}

double P_k_NSFD_only (double N_k, double h, double P_k,
                      double a, double b, double y, int number) {
    if (number > 0) {
        double N_k_1 = (h * N_k + N_k) * (N_k + a * P_k) / ((1 + h * N_k) * (N_k + a * P_k) + h * P_k);
        double P_k_1 = P_k * N_k * (1 + b * y * h) / (N_k + b * h * P_k);
        return P_k_NSFD_only(N_k_1, h, P_k_1, a, b, y, number - 1);
    } else {
        return P_k;
    }
}

double N_k_eyler (double N_k, double h, double P_k,
                  double a, double b, double y, int number) {
    if (number > 0) {
        double N_k_1 = N_k + h * N_k * (1 - N_k - P_k / (N_k + a * P_k));
        double P_k_1 = P_k + h * b * P_k * (y - P_k / N_k);
        return N_k_eyler(N_k_1, h, P_k_1, a, b, y, number - 1);
    } else {
        return N_k;
    }
}

double P_k_eyler (double N_k, double h, double P_k,
                  double a, double b, double y, int number) {
    if (number > 0) {
        double N_k_1 = N_k + h * N_k * (1 - N_k - P_k / (N_k + a * P_k));
        double P_k_1 = P_k + h * b * P_k * (y - P_k / N_k);
        return P_k_eyler(N_k_1, h, P_k_1, a, b, y, number - 1);
    } else {
        return P_k;
    }
}


int main() {
    cairo_surface_t *surface; // Определяем поверхность для рисования
    cairo_t *cr;              // Определяем источник




    bool NSFD_hyt = 1;
    bool NSFD_only = 1;
    bool eyler = 0;

    double N_size_x_min = 0.1;
    double N_size_x_pl = 0.8;

    double P_size_y_min = 0.1;
    double P_size_y_pl = 0.5;

    double a = 0.7;
    double b = 0.9;
    double y = 0.6;
    double h = 30000;

    double dot_size = 0.1;
    double gs = 100;





    double N_size = N_size_x_pl - N_size_x_min;
    double P_size = P_size_y_pl - P_size_y_min;

    double first_N = 0.2;
    double first_P = 0.2;
    /*
    double U_min = first_u;
    double U_max = first_u;

    for (double number = 0; (number * h) < h_max; number += 0.01) {
        double U_k = return_u(number, c);
        if (U_k < U_min) {
            U_min = U_k;
        }
        if (U_k > U_max) {
            U_max = U_k;
        }
    }
    for (int number = 0; (number * h) < h_max; number += 1) {
        double U_k = return_u_1_mic(first_u, h, a, b, number);
        if (U_k < U_min) {
            U_min = U_k;
        }
        if (U_k > U_max) {
            U_max = U_k;
        }
        U_k = return_u_1_sch(first_u, h, a, b, number);
        if (U_k < U_min) {
            U_min = U_k;
        }
        if (U_k > U_max) {
            U_max = U_k;
        }
        U_k = return_u_1_stan(first_u, h, a, b, number);
        if (U_k < U_min) {
            U_min = U_k;
        }
        if (U_k > U_max) {
            U_max = U_k;
        }
    }
    double U_size = U_max - U_min;
    double graph_size = 2;
    */

    surface = cairo_pdf_surface_create("result.pdf", N_size * gs, P_size * gs); // Создаем поверхность для рисования в виде файла PDF

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

    cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 166.0/255, 166.0/255, 166.0/255);
    cairo_set_line_width(cr, 0.0001);
    cairo_move_to(cr, abs(N_size_x_min) * gs, 0 * gs);
    cairo_line_to(cr, abs(N_size_x_min) * gs, P_size * gs);
    cairo_stroke(cr);
    cairo_destroy(cr);

    /*
    cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 166.0/255, 166.0/255, 166.0/255);
    cairo_set_line_width(cr, 0.0001);
    cairo_move_to(cr, (abs(N_size_x_min) + 0.2) * gs, 0 * gs);
    cairo_line_to(cr, (abs(N_size_x_min) + 0.2) * gs, P_size * gs);
    cairo_stroke(cr);
    cairo_destroy(cr);

    cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 166.0/255, 166.0/255, 166.0/255);
    cairo_set_line_width(cr, 0.0001);
    cairo_move_to(cr, 0 * gs, P_size_y_pl * gs);
    cairo_line_to(cr, N_size * gs, P_size_y_pl * gs);
    cairo_stroke(cr);
    cairo_destroy(cr);
    */

    cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 166.0/255, 166.0/255, 166.0/255);
    cairo_set_line_width(cr, 0.0001);
    cairo_move_to(cr, 0 * gs, (P_size_y_pl - 0.2) * gs);
    cairo_line_to(cr, N_size * gs, (P_size_y_pl - 0.2) * gs);
    cairo_stroke(cr);
    cairo_destroy(cr);

    /*
    cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 166.0/255, 166.0/255, 166.0/255);
    cairo_set_line_width(cr, 0.0001);
    cairo_move_to(cr, 0 * gs, (U_size_y_pl - 1) * gs);
    cairo_line_to(cr, h_size * gs, (U_size_y_pl - 1) * gs);
    cairo_stroke(cr);
    cairo_destroy(cr);
    */

    if (NSFD_hyt) {
        for (int number = 0; number < 250; number += 1) {
            double N_k_1 = N_k_NSFD_hyt(first_N, h, first_P, a, b, y, number);
            double P_k_1 = P_k_NSFD_hyt(first_N, h, first_P, a, b, y, number);
            double N_k_2 = N_k_NSFD_hyt(first_N, h, first_P, a, b, y, number + 1);
            double P_k_2 = P_k_NSFD_hyt(first_N, h, first_P, a, b, y, number + 1);
            cr = cairo_create(surface);
            cairo_set_source_rgb(cr, 0.0/255, 0.0/255, 255.0/255);
            cairo_set_line_width(cr, 0.001);
            cairo_arc(cr,
                      (N_k_1 - N_size_x_min) * gs,
                      (P_size - (P_k_1 - P_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
            cairo_fill(cr);
            cairo_destroy(cr);

            cr = cairo_create(surface); // рисование линии
            cairo_set_source_rgb(cr, 0.0/255, 0.0/255, 255.0/255);
            cairo_set_line_width(cr, dot_size);
            cairo_move_to(cr, (N_k_1 - N_size_x_min) * gs, (P_size - (P_k_1 - P_size_y_min)) * gs);
            cairo_line_to(cr, (N_k_2 - N_size_x_min) * gs, (P_size - (P_k_2 - P_size_y_min)) * gs);
            cairo_stroke(cr);
            cairo_destroy(cr);
        }
    }

    if (NSFD_only) {
        for (int number = 0; number < 250; number += 1) {
            double N_k_1 = N_k_NSFD_only(first_N, h, first_P, a, b, y, number);
            double P_k_1 = P_k_NSFD_only(first_N, h, first_P, a, b, y, number);
            double N_k_2 = N_k_NSFD_only(first_N, h, first_P, a, b, y, number + 1);
            double P_k_2 = P_k_NSFD_only(first_N, h, first_P, a, b, y, number + 1);
            cr = cairo_create(surface);
            cairo_set_source_rgb(cr, 0.0/255, 255.0/255, 0.0/255);
            cairo_set_line_width(cr, 0.001);
            cairo_arc(cr,
                      (N_k_1 - N_size_x_min) * gs,
                      (P_size - (P_k_1 - P_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
            cairo_fill(cr);
            cairo_destroy(cr);

            cr = cairo_create(surface); // рисование линии
            cairo_set_source_rgb(cr, 0.0/255, 255.0/255, 0.0/255);
            cairo_set_line_width(cr, dot_size);
            cairo_move_to(cr, (N_k_1 - N_size_x_min) * gs, (P_size - (P_k_1 - P_size_y_min)) * gs);
            cairo_line_to(cr, (N_k_2 - N_size_x_min) * gs, (P_size - (P_k_2 - P_size_y_min)) * gs);
            cairo_stroke(cr);
            cairo_destroy(cr);
        }
    }

    if (eyler) {
        for (int number = 0; number < 250; number += 1) {
            double N_k_1 = N_k_eyler(first_N, h, first_P, a, b, y, number);
            double P_k_1 = P_k_eyler(first_N, h, first_P, a, b, y, number);
            double N_k_2 = N_k_eyler(first_N, h, first_P, a, b, y, number + 1);
            double P_k_2 = P_k_eyler(first_N, h, first_P, a, b, y, number + 1);
            cr = cairo_create(surface);
            cairo_set_source_rgb(cr, 255.0/255, 0.0/255, 0.0/255);
            cairo_set_line_width(cr, 0.001);
            cairo_arc(cr,
                      (N_k_1 - N_size_x_min) * gs,
                      (P_size - (P_k_1 - P_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
            cairo_fill(cr);
            cairo_destroy(cr);

            cr = cairo_create(surface); // рисование линии
            cairo_set_source_rgb(cr, 255.0/255, 0.0/255, 0.0/255);
            cairo_set_line_width(cr, dot_size);
            cairo_move_to(cr, (N_k_1 - N_size_x_min) * gs, (P_size - (P_k_1 - P_size_y_min)) * gs);
            cairo_line_to(cr, (N_k_2 - N_size_x_min) * gs, (P_size - (P_k_2 - P_size_y_min)) * gs);
            cairo_stroke(cr);
            cairo_destroy(cr);
        }
    }


    cairo_surface_destroy(surface);
    return 0;
}
