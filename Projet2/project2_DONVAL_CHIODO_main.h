#ifndef MY_HEADER_MAIN  /* Evite la définition multiple */
#define MY_HEADER_MAIN
struct Param {
    double g;
    double gamma;
    double delta_x;
    double delta_y;
    double delta_t;
    double t_max;
    double A;
    double f;
    unsigned int S;
    unsigned int s;
    double r_threshold;
    int v_condition;
};

struct Param extract_parametres(char* file);
#endif