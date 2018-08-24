#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>

#include "input_data.h"

enum selection_mth_t {
    energy,
    norm
};

enum gauge_t {
    velocity,
    dipole
};

class Settings {
   public:
    static void create(Input_data& in);

    static Settings* get() { return instance; }

    bool _compute_bound;
    bool _compute_cont;

    double _energy_i;
    double _energy_f;

    std::string _path_in;
    std::string _file_basis;
    std::string _file_1E;
    std::string _file_2E;
    std::string _file_norm;
    std::string _file_HFv_i;
    std::string _file_CI_i;
    std::string _file_HFv_f;

    std::string _continuum_id;

    selection_mth_t _selection_m;
    gauge_t _gauge;
    bool _force_orth;

    bool _write;
    std::string _file_out;

   private:
    Settings() {}
    Settings(const Settings& old);
    const Settings& operator=(const Settings& old);
    ~Settings() {}

    static Settings* instance;

   private:
    bool _compute_bound;
    bool _compute_cont;

    double _energy_i;
    double _energy_f;

    std::string _path_in;
    std::string _file_basis;
    std::string _file_1E;
    std::string _file_2E;
    std::string _file_norm;
    std::string _file_HFv_i;
    std::string _file_CI_i;
    std::string _file_HFv_f;

    std::string _continuum_id;

    selection_mth_t _selection_m;
    gauge_t _gauge;
    bool _force_orth;

    bool _write;
    std::string _path_out;
    std::string _file_out;
};

#endif