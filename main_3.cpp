#include <math.h>
#include <iostream>
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>

using namespace std;

double return_u_1_standart(double U_k, double h, double a, double b, int number) {
    if (number > 0) {
        double U_k_1 = U_k + h * (a * U_k - b * U_k * U_k);
        return return_u_1_standart(U_k_1, h, a, b, number - 1);
    } else {
        return U_k;
    }
}

double return_u_1_non_stan_sch(double U_k, double h, double a, double b, int number) {
    if (number > 0) {
        double U_k_1 = (U_k + a * h * U_k) / (1 + h * b * U_k);
        return return_u_1_non_stan_sch(U_k_1, h, a, b, number - 1);
    } else {
        return U_k;
    }
}

double return_u_1_mickens(double U_k, double h, double a, double b, int number) {
    if (number > 0) {
        double U_k_1 = (U_k + a * ((exp(a * h) - 1) / a) * U_k) / (1 + b * U_k * ((exp(a * h) - 1) / a));
        return return_u_1_mickens(U_k_1, h, a, b, number - 1);
    } else {
        return U_k;
    }
}

double return_u(double t, double c) {
    if ((exp(t) + c) != 0) {
        return exp(t) / (exp(t) + c);
    }
    return 0;
}

int main() {
    cairo_surface_t *surface; // Определяем поверхность для рисования
    cairo_t *cr;              // Определяем источник




    double a = 1;
    double b = 1;

    bool graph = 1;
    bool stand_disk = 1;
    bool mickens = 0;
    bool non_stan_scheme = 0;

    double h_size_x_min = -0.25; // -0.25
    double h_size_x_pl = 3; // 3

    double U_size_y_min = -30; // -0.25
    double U_size_y_pl = 15; // 5

    double h = 0.31;
    double c = 1;
    double dot_size = 0.1;

    double gs = 10;




    double h_size = h_size_x_pl + abs(h_size_x_min);
    double U_size = U_size_y_pl + abs(U_size_y_min);

    double first_u = return_u(0, c);

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

    surface = cairo_pdf_surface_create("result.pdf", h_size * gs, U_size * gs); // Создаем поверхность для рисования в виде файла PDF

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
    cairo_move_to(cr, abs(h_size_x_min) * gs, 0 * gs);
    cairo_line_to(cr, abs(h_size_x_min) * gs, U_size * gs);
    cairo_stroke(cr);
    cairo_destroy(cr);

    cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 166.0/255, 166.0/255, 166.0/255);
    cairo_set_line_width(cr, 0.0001);
    cairo_move_to(cr, 0 * gs, U_size_y_pl * gs);
    cairo_line_to(cr, h_size * gs, U_size_y_pl * gs);
    cairo_stroke(cr);
    cairo_destroy(cr);

    cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 166.0/255, 166.0/255, 166.0/255);
    cairo_set_line_width(cr, 0.0001);
    cairo_move_to(cr, 0 * gs, (U_size_y_pl - 1) * gs);
    cairo_line_to(cr, h_size * gs, (U_size_y_pl - 1) * gs);
    cairo_stroke(cr);
    cairo_destroy(cr);

    if (graph) {
        for (double x = 0.1; x < 0; x += 0.1) {
            for (double number = 0; number < h_size_x_pl; number += 0.1) {
                double U_k_1 = return_u(number, (1.0 / x) - 1);
                double U_k_2 = return_u(number + 0.1, (1.0 / x) - 1);
                cr = cairo_create(surface);
                cairo_set_source_rgb(cr, 112.0/255, 128.0/255, 144.0/255);
                cairo_set_line_width(cr, 0.001);
                cairo_arc(cr,
                          (number - h_size_x_min) * gs,
                          (U_size - (U_k_1 - U_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
                cairo_fill(cr);
                cairo_destroy(cr);

                cr = cairo_create(surface); // рисование линии
                cairo_set_source_rgb(cr, 112.0/255, 128.0/255, 144.0/255);
                cairo_set_line_width(cr, dot_size);
                cairo_move_to(cr, (number - h_size_x_min) * gs, (U_size - (U_k_1 - U_size_y_min)) * gs);
                cairo_line_to(cr, ((number + 0.1) - h_size_x_min) * gs, (U_size - (U_k_2 - U_size_y_min)) * gs);
                cairo_stroke(cr);
                cairo_destroy(cr);
            }
        }
        for (double c = 0; c > -0.9; c -= 0.1) {
            for (double number = 0; number < h_size_x_pl; number += 0.1) {
                double U_k_1 = return_u(number, c);
                double U_k_2 = return_u(number + 0.1, c);
                cr = cairo_create(surface);
                cairo_set_source_rgb(cr, 112.0/255, 128.0/255, 144.0/255);
                cairo_set_line_width(cr, 0.001);
                cairo_arc(cr,
                          (number - h_size_x_min) * gs,
                          (U_size - (U_k_1 - U_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
                cairo_fill(cr);
                cairo_destroy(cr);

                cr = cairo_create(surface); // рисование линии
                cairo_set_source_rgb(cr, 112.0/255, 128.0/255, 144.0/255);
                cairo_set_line_width(cr, dot_size);
                cairo_move_to(cr, (number - h_size_x_min) * gs, (U_size - (U_k_1 - U_size_y_min)) * gs);
                cairo_line_to(cr, ((number + 0.1) - h_size_x_min) * gs, (U_size - (U_k_2 - U_size_y_min)) * gs);
                cairo_stroke(cr);
                cairo_destroy(cr);
            }
        }
    }

    if (stand_disk) {
        for (double x = 0.1; x < 0; x += 0.1) {
            first_u = return_u(0, (1.0 / x) - 1);
            for (int number = 0; (number * h) < h_size_x_pl; number += 1) {
                double U_k_1 = return_u_1_standart(first_u, h, a, b, number);
                double U_k_2 = return_u_1_standart(first_u, h, a, b, number + 1);
                cr = cairo_create(surface);
                cairo_set_source_rgb(cr, 0, 1, 0);
                cairo_set_line_width(cr, 0.001);
                cairo_arc(cr,
                          (number * h - h_size_x_min) * gs,
                          (U_size - (U_k_1 - U_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
                cairo_fill(cr);
                cairo_destroy(cr);

                if ((U_size - (U_k_2 - U_size_y_min)) < U_size) {
                    cr = cairo_create(surface); // рисование линии
                    cairo_set_source_rgb(cr, 0, 1, 0);
                    cairo_set_line_width(cr, dot_size);
                    cairo_move_to(cr, (number * h - h_size_x_min) * gs, (U_size - (U_k_1 - U_size_y_min)) * gs);
                    cairo_line_to(cr, ((number + 1) * h - h_size_x_min) * gs, (U_size - (U_k_2 - U_size_y_min)) * gs);
                    cairo_stroke(cr);
                    cairo_destroy(cr);
                }
            }
        }
        for (double c = 0; c > -0.9; c -= 0.1) {
            first_u = return_u(0, c);
            for (int number = 0; (number * h) < h_size_x_pl; number += 1) {
                double U_k_1 = return_u_1_standart(first_u, h, a, b, number);
                double U_k_2 = return_u_1_standart(first_u, h, a, b, number + 1);
                cr = cairo_create(surface);
                cairo_set_source_rgb(cr, 0, 1, 0);
                cairo_set_line_width(cr, 0.001);
                cairo_arc(cr,
                          (number * h - h_size_x_min) * gs,
                          (U_size - (U_k_1 - U_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
                cairo_fill(cr);
                cairo_destroy(cr);

                if ((U_size - (U_k_2 - U_size_y_min)) < U_size) {
                    cr = cairo_create(surface); // рисование линии
                    cairo_set_source_rgb(cr, 0, 1, 0);
                    cairo_set_line_width(cr, dot_size);
                    cairo_move_to(cr, (number * h - h_size_x_min) * gs, (U_size - (U_k_1 - U_size_y_min)) * gs);
                    cairo_line_to(cr, ((number + 1) * h - h_size_x_min) * gs, (U_size - (U_k_2 - U_size_y_min)) * gs);
                    cairo_stroke(cr);
                    cairo_destroy(cr);
                }
            }
        }
    }

    if (mickens) {
        for (double x = 0.1; x < 1; x += 0.1) {
            first_u = return_u(0, (1.0 / x) - 1);
            for (int number = 0; (number * h) < h_size_x_pl; number += 1) {
                double U_k_1 = return_u_1_mickens(first_u, h, a, b, number);
                double U_k_2 = return_u_1_mickens(first_u, h, a, b, number + 1);
                cr = cairo_create(surface);
                cairo_set_source_rgb(cr, 255.0/255, 0.0/255, 0.0/255);
                cairo_set_line_width(cr, 0.001);
                cairo_arc(cr,
                          (number * h - h_size_x_min) * gs,
                          (U_size - (U_k_1 - U_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
                cairo_fill(cr);
                cairo_destroy(cr);

                if ((U_size - (U_k_2 - U_size_y_min)) < U_size) {
                    cr = cairo_create(surface); // рисование линии
                    cairo_set_source_rgb(cr, 1, 0, 0);
                    cairo_set_line_width(cr, dot_size);
                    cairo_move_to(cr, (number * h - h_size_x_min) * gs, (U_size - (U_k_1 - U_size_y_min)) * gs);
                    cairo_line_to(cr, ((number + 1) * h - h_size_x_min) * gs, (U_size - (U_k_2 - U_size_y_min)) * gs);
                    cairo_stroke(cr);
                    cairo_destroy(cr);
                }
            }
        }
        for (double c = 0; c > -0.9; c -= 0.1) {
            first_u = return_u(0, c);
            for (int number = 0; (number * h) < h_size_x_pl; number += 1) {
                double U_k_1 = return_u_1_mickens(first_u, h, a, b, number);
                double U_k_2 = return_u_1_mickens(first_u, h, a, b, number + 1);
                cr = cairo_create(surface);
                cairo_set_source_rgb(cr, 255.0/255, 0.0/255, 0.0/255);
                cairo_set_line_width(cr, 0.001);
                cairo_arc(cr,
                          (number * h - h_size_x_min) * gs,
                          (U_size - (U_k_1 - U_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
                cairo_fill(cr);
                cairo_destroy(cr);

                if ((U_size - (U_k_2 - U_size_y_min)) < U_size) {
                    cr = cairo_create(surface); // рисование линии
                    cairo_set_source_rgb(cr, 1, 0, 0);
                    cairo_set_line_width(cr, dot_size);
                    cairo_move_to(cr, (number * h - h_size_x_min) * gs, (U_size - (U_k_1 - U_size_y_min)) * gs);
                    cairo_line_to(cr, ((number + 1) * h - h_size_x_min) * gs, (U_size - (U_k_2 - U_size_y_min)) * gs);
                    cairo_stroke(cr);
                    cairo_destroy(cr);
                }
            }
        }
    }

    if (non_stan_scheme) {
        for (double x = 0.1; x < 1; x += 0.1) {
            first_u = return_u(0, (1.0 / x) - 1);
            for (int number = 0; (number * h) < h_size_x_pl; number += 1) {
                double U_k_1 = return_u_1_non_stan_sch(first_u, h, a, b, number);
                double U_k_2 = return_u_1_non_stan_sch(first_u, h, a, b, number + 1);
                cr = cairo_create(surface);
                cairo_set_source_rgb(cr, 0.0/255, 0.0/255, 255.0/255);
                cairo_set_line_width(cr, 0.001);
                cairo_arc(cr,
                          (number * h - h_size_x_min) * gs,
                          (U_size - (U_k_1 - U_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
                cairo_fill(cr);
                cairo_destroy(cr);

                if ((U_size - (U_k_2 - U_size_y_min)) < U_size) {
                    cr = cairo_create(surface); // рисование линии
                    cairo_set_source_rgb(cr, 0.0/255, 0.0/255, 255.0/255);
                    cairo_set_line_width(cr, dot_size);
                    cairo_move_to(cr, (number * h - h_size_x_min) * gs, (U_size - (U_k_1 - U_size_y_min)) * gs);
                    cairo_line_to(cr, ((number + 1) * h - h_size_x_min) * gs, (U_size - (U_k_2 - U_size_y_min)) * gs);
                    cairo_stroke(cr);
                    cairo_destroy(cr);
                }
            }
        }
        for (double c = 0; c > -0.9; c -= 0.1) {
            first_u = return_u(0, c);
            for (int number = 0; (number * h) < h_size_x_pl; number += 1) {
                double U_k_1 = return_u_1_non_stan_sch(first_u, h, a, b, number);
                double U_k_2 = return_u_1_non_stan_sch(first_u, h, a, b, number + 1);
                cr = cairo_create(surface);
                cairo_set_source_rgb(cr, 0.0/255, 0.0/255, 255.0/255);
                cairo_set_line_width(cr, 0.001);
                cairo_arc(cr,
                          (number * h - h_size_x_min) * gs,
                          (U_size - (U_k_1 - U_size_y_min)) * gs, dot_size, 0, 2*3.14159265359);
                cairo_fill(cr);
                cairo_destroy(cr);

                if ((U_size - (U_k_2 - U_size_y_min)) < U_size) {
                    cr = cairo_create(surface); // рисование линии
                    cairo_set_source_rgb(cr, 0.0/255, 0.0/255, 255.0/255);
                    cairo_set_line_width(cr, dot_size);
                    cairo_move_to(cr, (number * h - h_size_x_min) * gs, (U_size - (U_k_1 - U_size_y_min)) * gs);
                    cairo_line_to(cr, ((number + 1) * h - h_size_x_min) * gs, (U_size - (U_k_2 - U_size_y_min)) * gs);
                    cairo_stroke(cr);
                    cairo_destroy(cr);
                }
            }
        }
    }

    cairo_surface_destroy(surface);
    return 0;
}
