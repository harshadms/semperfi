/*!
 * \file dll_pll_conf_fpga.h
 * \brief Class that contains all the configuration parameters for generic
 * tracking block based on a DLL and a PLL for the FPGA.
 * \author Marc Majoral, 2019. mmajoral(at)cttc.cat
 * \author Javier Arribas, 2018. jarribas(at)cttc.es
 *
 * Class that contains all the configuration parameters for generic tracking block based on a DLL and a PLL.
 *
 * -----------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2020  (see AUTHORS file for a list of contributors)
 *
 * GNSS-SDR is a software defined Global Navigation
 *          Satellite Systems receiver
 *
 * This file is part of GNSS-SDR.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * -----------------------------------------------------------------------------
 */

#ifndef GNSS_SDR_DLL_PLL_CONF_FPGA_H
#define GNSS_SDR_DLL_PLL_CONF_FPGA_H

#include "configuration_interface.h"
#include <cstdint>
#include <string>

class Dll_Pll_Conf_Fpga
{
public:
    Dll_Pll_Conf_Fpga();
    void SetFromConfiguration(const ConfigurationInterface* configuration, const std::string& role);

    /* DLL/PLL tracking configuration */
    std::string device_name;
    std::string dump_filename;

    double fs_in;
    double carrier_lock_th;

    float pll_pull_in_bw_hz;
    float dll_pull_in_bw_hz;
    float fll_bw_hz;
    float pll_bw_hz;
    float dll_bw_hz;
    float pll_bw_narrow_hz;
    float dll_bw_narrow_hz;
    float early_late_space_chips;
    float very_early_late_space_chips;
    float early_late_space_narrow_chips;
    float very_early_late_space_narrow_chips;
    float slope;
    float spc;
    float y_intercept;
    float cn0_smoother_alpha;
    float carrier_lock_test_smoother_alpha;

    uint32_t pull_in_time_s;  // signed integer, when pull in time is not yet reached it has to be compared against a negative number
    uint32_t bit_synchronization_time_limit_s;
    uint32_t vector_length;
    uint32_t smoother_length;
    uint32_t dev_file_num;
    uint32_t num_prev_assigned_ch;
    uint32_t code_length_chips;
    uint32_t code_samples_per_chip;
    uint32_t extend_fpga_integration_periods;
    uint32_t fpga_integration_period;

    int32_t fll_filter_order;
    int32_t pll_filter_order;
    int32_t dll_filter_order;
    int32_t extend_correlation_symbols;
    int32_t cn0_samples;
    int32_t cn0_min;
    int32_t max_code_lock_fail;
    int32_t max_carrier_lock_fail;
    int32_t cn0_smoother_samples;
    int32_t carrier_lock_test_smoother_samples;
    // int32_t max_lock_fail;

    int32_t* ca_codes;
    int32_t* data_codes;

    char signal[3];
    char system;

    bool extended_correlation_in_fpga;
    bool track_pilot;
    bool enable_doppler_correction;
    bool enable_fll_pull_in;
    bool enable_fll_steady_state;
    bool carrier_aiding;
    bool high_dyn;
    bool dump;
    bool dump_mat;
};

#endif
