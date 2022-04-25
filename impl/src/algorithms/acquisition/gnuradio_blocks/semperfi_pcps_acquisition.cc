/*!
 * \file semperfi_pcps_acquisition.cc
 * \brief This class implements a Parallel Code Phase Search Acquisition
 * \authors <ul>
 *          <li> Javier Arribas, 2011. jarribas(at)cttc.es
 *          <li> Luis Esteve, 2012. luis(at)epsilon-formacion.com
 *          <li> Marc Molina, 2013. marc.molina.pena@gmail.com
 *          <li> Cillian O'Driscoll, 2017. cillian(at)ieee.org
 *          </ul>
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

// includes for signal generator
#include "GPS_L1_CA.h"
#include "gps_sdr_signal_processing.h"
#include "Galileo_E1.h"
#include "Galileo_E5a.h"
#include "Galileo_E5b.h"
#include "galileo_e1_signal_processing.h"
#include "galileo_e5_signal_processing.h"

#include "semperfi_pcps_acquisition.h"
#include "GLONASS_L1_L2_CA.h"  // for GLONASS_PRN
#include "MATH_CONSTANTS.h"    // for TWO_PI
#include "gnss_frequencies.h"
#include "gnss_sdr_create_directory.h"
#include "gnss_sdr_make_unique.h"
#include "gnss_synchro.h"
#if HAS_STD_FILESYSTEM
#if HAS_STD_FILESYSTEM_EXPERIMENTAL
#include <experimental/filesystem>
#else
#include <filesystem>
#endif
#else
#include <boost/filesystem/path.hpp>
#endif
#include <boost/math/special_functions/gamma.hpp>
#include <gnuradio/io_signature.h>
#include <matio.h>
#include <pmt/pmt.h>        // for from_long
#include <pmt/pmt_sugar.h>  // for mp
#include <volk/volk.h>
#include <volk_gnsssdr/volk_gnsssdr.h>
#include <algorithm>  // for fill_n, min
#include <array>
#include <cmath>    // for floor, fmod, rint, ceil
#include <cstring>  // for memcpy
#include <iostream>
#include <map> 
#include "persistence1d.hpp"

#if HAS_STD_FILESYSTEM
#if HAS_STD_FILESYSTEM_EXPERIMENTAL
namespace fs = std::experimental::filesystem;
#else
namespace fs = std::filesystem;
#endif
#else
namespace fs = boost::filesystem;
#endif


semperfi_pcps_acquisition_sptr semperfi_pcps_make_acquisition(const Acq_Conf& conf_)
{
    return semperfi_pcps_acquisition_sptr(new semperfi_pcps_acquisition(conf_));
}


semperfi_pcps_acquisition::semperfi_pcps_acquisition(const Acq_Conf& conf_) : gr::block("semperfi_pcps_acquisition",
                                                                gr::io_signature::make(1, 1, conf_.it_size),
                                                                gr::io_signature::make(0, 0, conf_.it_size))
{
    this->message_port_register_out(pmt::mp("events"));

    d_acq_parameters = conf_;
    d_sample_counter = 0ULL;  // SAMPLE COUNTER
    d_active = false;
    d_positive_acq = 0;
    d_acq_state = 0;
    d_doppler_bias = 0;
    d_num_noncoherent_integrations_counter = 0U;
    d_consumed_samples = d_acq_parameters.sampled_ms * d_acq_parameters.samples_per_ms * (d_acq_parameters.bit_transition_flag ? 2.0 : 1.0);
    if (d_acq_parameters.sampled_ms == d_acq_parameters.ms_per_code)
        {
            d_fft_size = d_consumed_samples;
        }
    else
        {
            d_fft_size = d_consumed_samples * 2;
        }
    // d_fft_size = next power of two?  ////
    d_mag = 0;
    d_input_power = 0.0;
    d_num_doppler_bins = 0U;
    d_threshold = 0.0;
    d_doppler_step = d_acq_parameters.doppler_step;
    d_doppler_center = 0U;
    d_doppler_center_step_two = 0.0;
    d_test_statistics = 0.0;
    d_channel = 0U;
    if (conf_.it_size == sizeof(gr_complex))
    {
        d_cshort = false;
    }
    else
    {
        d_cshort = true;
    }

    // COD:
    // Experimenting with the overlap/save technique for handling bit trannsitions
    // The problem: Circular correlation is asynchronous with the received code.
    // In effect the first code phase used in the correlation is the current
    // estimate of the code phase at the start of the input buffer. If this is 1/2
    // of the code period a bit transition would move all the signal energy into
    // adjacent frequency bands at +/- 1/T where T is the integration time.
    //
    // We can avoid this by doing linear correlation, effectively doubling the
    // size of the input buffer and padding the code with zeros.
    // if (d_acq_parameters.bit_transition_flag)
    // {
    //  d_fft_size = d_consumed_samples * 2;
    //  d_acq_parameters.max_dwells = 1;  // Activation of d_acq_parameters.bit_transition_flag invalidates the value of d_acq_parameters.max_dwells
    // }

    d_tmp_buffer = volk_gnsssdr::vector<float>(d_fft_size);
    d_fft_codes = volk_gnsssdr::vector<std::complex<float>>(d_fft_size);
    d_input_signal = volk_gnsssdr::vector<std::complex<float>>(d_fft_size);

    // Direct FFT
    d_fft_if = std::make_unique<gr::fft::fft_complex>(d_fft_size, true);

    // Inverse FFT
    d_ifft = std::make_unique<gr::fft::fft_complex>(d_fft_size, false);

    d_gnss_synchro = nullptr;
    d_worker_active = false;
    d_data_buffer = volk_gnsssdr::vector<std::complex<float>>(d_consumed_samples);

    // Spoofing and recovery
    d_recovery_signal_buff = volk_gnsssdr::vector<std::complex<float>>(d_consumed_samples);
    d_codes_generated = false;
    d_code_delay_diff = 0;
    d_itr = 0;
    d_global_itr = 0;
    d_legit_code_delay = 0;
    d_adv_code_delay = 0;
    d_reset_time = true;

    // generate_codes()

    if (d_cshort)
        {
            d_data_buffer_sc = volk_gnsssdr::vector<lv_16sc_t>(d_consumed_samples);
        }

    d_grid = arma::fmat();
    d_narrow_grid = arma::fmat();
    d_step_two = false;
    d_num_doppler_bins_step2 = d_acq_parameters.num_doppler_bins_step2;

    d_samplesPerChip = d_acq_parameters.samples_per_chip;
    d_buffer_count = 0U;
    d_use_CFAR_algorithm_flag = d_acq_parameters.use_CFAR_algorithm_flag;
    d_dump_number = 0LL;
    d_dump_channel = d_acq_parameters.dump_channel;
    d_dump = d_acq_parameters.dump;
    d_dump_filename = d_acq_parameters.dump_filename;

    semperfi_start = std::chrono::high_resolution_clock::now();

    if (d_dump)
        {
            std::string dump_path;
            // Get path
            if (d_dump_filename.find_last_of('/') != std::string::npos)
                {
                    const std::string dump_filename_ = d_dump_filename.substr(d_dump_filename.find_last_of('/') + 1);
                    dump_path = d_dump_filename.substr(0, d_dump_filename.find_last_of('/'));
                    d_dump_filename = dump_filename_;
                }
            else
                {
                    dump_path = std::string(".");
                }
            if (d_dump_filename.empty())
                {
                    d_dump_filename = "acquisition";
                }
            // remove extension if any
            if (d_dump_filename.substr(1).find_last_of('.') != std::string::npos)
                {
                    d_dump_filename = d_dump_filename.substr(0, d_dump_filename.find_last_of('.'));
                }
            d_dump_filename = dump_path + fs::path::preferred_separator + d_dump_filename;
            // create directory
            if (!gnss_sdr_create_directory(dump_path))
                {
                    std::cerr << "GNSS-SDR cannot create dump file for the Acquisition block. Wrong permissions?\n";
                    d_dump = false;
                }
        }
}


void semperfi_pcps_acquisition::set_resampler_latency(uint32_t latency_samples)
{
    gr::thread::scoped_lock lock(d_setlock);  // require mutex with work function called by the scheduler
    d_acq_parameters.resampler_latency_samples = latency_samples;
}


void semperfi_pcps_acquisition::set_local_code(std::complex<float>* code)
{
    // This will check if it's fdma, if yes will update the intermediate frequency and the doppler grid
    if (is_fdma())
        {
            update_grid_doppler_wipeoffs();
        }
    // COD
    // Here we want to create a buffer that looks like this:
    // [ 0 0 0 ... 0 c_0 c_1 ... c_L]
    // where c_i is the local code and there are L zeros and L chips
    gr::thread::scoped_lock lock(d_setlock);  // require mutex with work function called by the scheduler
    if (d_acq_parameters.bit_transition_flag)
        {
            const int32_t offset = d_fft_size / 2;
            std::fill_n(d_fft_if->get_inbuf(), offset, gr_complex(0.0, 0.0));
            memcpy(d_fft_if->get_inbuf() + offset, code, sizeof(gr_complex) * offset);
        }
    else
        {
            if (d_acq_parameters.sampled_ms == d_acq_parameters.ms_per_code)
                {
                    memcpy(d_fft_if->get_inbuf(), code, sizeof(gr_complex) * d_consumed_samples);
                }
            else
                {
                    std::fill_n(d_fft_if->get_inbuf(), d_fft_size - d_consumed_samples, gr_complex(0.0, 0.0));
                    memcpy(d_fft_if->get_inbuf() + d_consumed_samples, code, sizeof(gr_complex) * d_consumed_samples);
                }
        }

    d_fft_if->execute();  // We need the FFT of local code
    volk_32fc_conjugate_32fc(d_fft_codes.data(), d_fft_if->get_outbuf(), d_fft_size);
}


bool semperfi_pcps_acquisition::is_fdma()
{
    // reset the intermediate frequency
    d_doppler_bias = 0;
    // Dealing with FDMA system
    if (strcmp(d_gnss_synchro->Signal, "1G") == 0)
        {
            d_doppler_bias = static_cast<int32_t>(DFRQ1_GLO * GLONASS_PRN.at(d_gnss_synchro->PRN));
            DLOG(INFO) << "Trying to acquire SV PRN " << d_gnss_synchro->PRN << " with freq " << d_doppler_bias << " in Glonass Channel " << GLONASS_PRN.at(d_gnss_synchro->PRN) << '\n';
            return true;
        }
    if (strcmp(d_gnss_synchro->Signal, "2G") == 0)
        {
            d_doppler_bias += static_cast<int32_t>(DFRQ2_GLO * GLONASS_PRN.at(d_gnss_synchro->PRN));
            DLOG(INFO) << "Trying to acquire SV PRN " << d_gnss_synchro->PRN << " with freq " << d_doppler_bias << " in Glonass Channel " << GLONASS_PRN.at(d_gnss_synchro->PRN) << '\n';
            return true;
        }
    return false;
}


void semperfi_pcps_acquisition::update_local_carrier(own::span<gr_complex> carrier_vector, float freq)
{
    float phase_step_rad;
    if (d_acq_parameters.use_automatic_resampler)
        {
            phase_step_rad = static_cast<float>(TWO_PI) * freq / static_cast<float>(d_acq_parameters.resampled_fs);
        }
    else
        {
            phase_step_rad = static_cast<float>(TWO_PI) * freq / static_cast<float>(d_acq_parameters.fs_in);
        }
    std::array<float, 1> _phase{};
    volk_gnsssdr_s32f_sincos_32fc(carrier_vector.data(), -phase_step_rad, _phase.data(), carrier_vector.size());
}


void semperfi_pcps_acquisition::init()
{
    d_gnss_synchro->Flag_valid_acquisition = false;
    d_gnss_synchro->Flag_valid_symbol_output = false;
    d_gnss_synchro->Flag_valid_pseudorange = false;
    d_gnss_synchro->Flag_valid_word = false;
    d_gnss_synchro->Acq_doppler_step = 0U;
    d_gnss_synchro->Acq_delay_samples = 0.0;
    d_gnss_synchro->Acq_doppler_hz = 0.0;
    d_gnss_synchro->Acq_samplestamp_samples = 0ULL;
    d_mag = 0.0;
    d_input_power = 0.0;

    d_num_doppler_bins = static_cast<uint32_t>(std::ceil(static_cast<double>(static_cast<int32_t>(d_acq_parameters.doppler_max) - static_cast<int32_t>(-d_acq_parameters.doppler_max)) / static_cast<double>(d_doppler_step)));

    // Spoofing stuff
    d_spoofer_present = false;
    d_spoofer_detected = false;
    d_repeat_acq = false;
    d_restart_sic = false;
    d_itr = 0;

    // Create the carrier Doppler wipeoff signals
    if (d_grid_doppler_wipeoffs.empty())
        {
            d_grid_doppler_wipeoffs = volk_gnsssdr::vector<volk_gnsssdr::vector<std::complex<float>>>(d_num_doppler_bins, volk_gnsssdr::vector<std::complex<float>>(d_fft_size));
        }
    if (d_acq_parameters.make_2_steps && (d_grid_doppler_wipeoffs_step_two.empty()))
        {
            d_grid_doppler_wipeoffs_step_two = volk_gnsssdr::vector<volk_gnsssdr::vector<std::complex<float>>>(d_num_doppler_bins_step2, volk_gnsssdr::vector<std::complex<float>>(d_fft_size));
        }

    if (d_magnitude_grid.empty())
        {
            d_magnitude_grid = volk_gnsssdr::vector<volk_gnsssdr::vector<float>>(d_num_doppler_bins, volk_gnsssdr::vector<float>(d_fft_size));
        }

    for (uint32_t doppler_index = 0; doppler_index < d_num_doppler_bins; doppler_index++)
        {
            std::fill(d_magnitude_grid[doppler_index].begin(), d_magnitude_grid[doppler_index].end(), 0.0);
        }

    update_grid_doppler_wipeoffs();
    d_worker_active = false;

    if (d_dump)
        {
            const uint32_t effective_fft_size = (d_acq_parameters.bit_transition_flag ? (d_fft_size / 2) : d_fft_size);
            d_grid = arma::fmat(effective_fft_size, d_num_doppler_bins, arma::fill::zeros);
            d_narrow_grid = arma::fmat(effective_fft_size, d_num_doppler_bins_step2, arma::fill::zeros);
        }
}


void semperfi_pcps_acquisition::update_grid_doppler_wipeoffs()
{
    for (uint32_t doppler_index = 0; doppler_index < d_num_doppler_bins; doppler_index++)
        {
            const int32_t doppler = -static_cast<int32_t>(d_acq_parameters.doppler_max) + d_doppler_center + d_doppler_step * doppler_index;
            update_local_carrier(d_grid_doppler_wipeoffs[doppler_index], static_cast<float>(d_doppler_bias + doppler));
        }
}


void semperfi_pcps_acquisition::update_grid_doppler_wipeoffs_step2()
{
    for (uint32_t doppler_index = 0; doppler_index < d_num_doppler_bins_step2; doppler_index++)
        {
            const float doppler = (static_cast<float>(doppler_index) - static_cast<float>(floor(d_num_doppler_bins_step2 / 2.0))) * d_acq_parameters.doppler_step2;
            update_local_carrier(d_grid_doppler_wipeoffs_step_two[doppler_index], d_doppler_center_step_two + doppler);
        }
}


void semperfi_pcps_acquisition::set_state(int32_t state)
{
    gr::thread::scoped_lock lock(d_setlock);  // require mutex with work function called by the scheduler
    d_acq_state = state;
    if (d_acq_state == 1)
        {
            d_gnss_synchro->Acq_delay_samples = 0.0;
            d_gnss_synchro->Acq_doppler_hz = 0.0;
            d_gnss_synchro->Acq_samplestamp_samples = 0ULL;
            d_gnss_synchro->Acq_doppler_step = 0U;
            d_mag = 0.0;
            d_test_statistics = 0.0;
            d_active = true;
        }
    else if (d_acq_state == 0)
        {
        }
    else
        {
            LOG(ERROR) << "State can only be set to 0 or 1";
        }
}


void semperfi_pcps_acquisition::send_positive_acquisition()
{
    //std::cout << "394 pos acq";
    // Declare positive acquisition using a message port
    // 0=STOP_CHANNEL 1=ACQ_SUCCEES 2=ACQ_FAIL
    DLOG(INFO) << "positive acquisition"
               << ", satellite " << d_gnss_synchro->System << " " << d_gnss_synchro->PRN
               << ", sample_stamp " << d_sample_counter
               << ", test statistics value " << d_test_statistics
               << ", test statistics threshold " << d_threshold
               << ", code phase " << d_gnss_synchro->Acq_delay_samples
               << ", doppler " << d_gnss_synchro->Acq_doppler_hz
               << ", magnitude " << d_mag
               << ", input signal power " << d_input_power
               << ", Assist doppler_center " << d_doppler_center;
    d_positive_acq = 1;

    auto elapsed = std::chrono::high_resolution_clock::now() - semperfi_start;
    int time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    d_reset_time = true;
    

    std::cout << "\n****** SV " << d_gnss_synchro->PRN << " Legit params found in " << d_global_itr << " iterations; Time: " << time/1e6 << "  secs ******\n";
    d_global_itr = 0;
    //std::cout << "\nLegitimate code delay: " << d_gnss_synchro->Acq_delay_samples << " Doppler: " << d_gnss_synchro->Acq_doppler_hz << "\n";

    if (!d_channel_fsm.expired())
        {
            // the channel FSM is set, so, notify it directly the positive acquisition to minimize delays
            d_channel_fsm.lock()->Event_valid_acquisition();
        }
    else
        {
            this->message_port_pub(pmt::mp("events"), pmt::from_long(1));
        }
    dump_time();
}


void semperfi_pcps_acquisition::send_negative_acquisition()
{
    //std::cout << "430 Neg acq";
    // Declare negative acquisition using a message port
    // 0=STOP_CHANNEL 1=ACQ_SUCCEES 2=ACQ_FAIL
    DLOG(INFO) << "negative acquisition"
               << ", satellite " << d_gnss_synchro->System << " " << d_gnss_synchro->PRN
               << ", sample_stamp " << d_sample_counter
               << ", test statistics value " << d_test_statistics
               << ", test statistics threshold " << d_threshold
               << ", code phase " << d_gnss_synchro->Acq_delay_samples
               << ", doppler " << d_gnss_synchro->Acq_doppler_hz
               << ", magnitude " << d_mag
               << ", input signal power " << d_input_power;
    d_positive_acq = 0;
    //std::cout << "\nNegative acq PRN " << d_gnss_synchro -> PRN << "\n";
    this->message_port_pub(pmt::mp("events"), pmt::from_long(2));
}


void semperfi_pcps_acquisition::dump_time()
{
    std::chrono::high_resolution_clock::time_point curr_time =std::chrono::high_resolution_clock::now();
    int time = std::chrono::duration_cast<std::chrono::microseconds>(curr_time - semperfi_start).count();
    std::string filename = "/tmp/semperfi_time";
    //filename.append(d_acq_parameters.time_stats_filename);

    std::cout << "\nfilename: " << filename;
    std::ofstream outfile;

    outfile.open(filename, std::ios_base::app); // append instead of overwrite

    std::string csv_string = "\n";
    csv_string.append(std::to_string(d_gnss_synchro->PRN));
    csv_string.append(",");
    csv_string.append(std::to_string(time / 1e6));
    csv_string.append(",");
    csv_string.append(std::to_string(d_global_itr+1));

    //std::cout << "\ncsv_string : " << csv_string;
    outfile << csv_string;

}


void semperfi_pcps_acquisition::dump_results(int32_t effective_fft_size)
{

    d_dump_number++;
    std::string filename = d_dump_filename;
    filename.append("_");
    filename.append(1, d_gnss_synchro->System);
    filename.append("_");
    filename.append(1, d_gnss_synchro->Signal[0]);
    filename.append(1, d_gnss_synchro->Signal[1]);
    filename.append("_ch_");
    filename.append(std::to_string(d_channel));
    filename.append("_");
    filename.append(std::to_string(d_dump_number));
    filename.append("_sat_");
    filename.append(std::to_string(d_gnss_synchro->PRN));
    filename.append(".mat");

    mat_t* matfp = Mat_CreateVer(filename.c_str(), nullptr, MAT_FT_MAT73);
    if (matfp == nullptr)
        {
            std::cout << "Unable to create or open Acquisition dump file\n";
            // d_acq_parameters.dump = false;
        }
    else
        {
            /* Dump Standard Acquisition Parameters */
            std::array<size_t, 2> dims{static_cast<size_t>(effective_fft_size), static_cast<size_t>(d_num_doppler_bins)};
            matvar_t* matvar = Mat_VarCreate("acq_grid", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims.data(), d_grid.memptr(), 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            dims[0] = static_cast<size_t>(1);
            dims[1] = static_cast<size_t>(1);
            matvar = Mat_VarCreate("doppler_max", MAT_C_INT32, MAT_T_INT32, 1, dims.data(), &d_acq_parameters.doppler_max, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            matvar = Mat_VarCreate("doppler_step", MAT_C_INT32, MAT_T_INT32, 1, dims.data(), &d_doppler_step, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            matvar = Mat_VarCreate("d_positive_acq", MAT_C_INT32, MAT_T_INT32, 1, dims.data(), &d_positive_acq, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            auto aux = static_cast<float>(d_gnss_synchro->Acq_doppler_hz);
            matvar = Mat_VarCreate("acq_doppler_hz", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &aux, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            aux = static_cast<float>(d_gnss_synchro->Acq_delay_samples);
            matvar = Mat_VarCreate("acq_delay_samples", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &aux, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            matvar = Mat_VarCreate("test_statistic", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &d_test_statistics, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            matvar = Mat_VarCreate("threshold", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &d_threshold, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            matvar = Mat_VarCreate("input_power", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &d_input_power, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            matvar = Mat_VarCreate("sample_counter", MAT_C_UINT64, MAT_T_UINT64, 1, dims.data(), &d_sample_counter, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            matvar = Mat_VarCreate("PRN", MAT_C_UINT32, MAT_T_UINT32, 1, dims.data(), &d_gnss_synchro->PRN, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            matvar = Mat_VarCreate("num_dwells", MAT_C_INT32, MAT_T_INT32, 1, dims.data(), &d_num_noncoherent_integrations_counter, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

            if (d_acq_parameters.make_2_steps)
                {
                    dims[0] = static_cast<size_t>(effective_fft_size);
                    dims[1] = static_cast<size_t>(d_num_doppler_bins_step2);
                    matvar = Mat_VarCreate("acq_grid_narrow", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims.data(), d_narrow_grid.memptr(), 0);
                    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                    Mat_VarFree(matvar);

                    dims[0] = static_cast<size_t>(1);
                    dims[1] = static_cast<size_t>(1);
                    matvar = Mat_VarCreate("doppler_step_narrow", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &d_acq_parameters.doppler_step2, 0);
                    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                    Mat_VarFree(matvar);

                    aux = d_doppler_center_step_two - static_cast<float>(floor(d_num_doppler_bins_step2 / 2.0)) * d_acq_parameters.doppler_step2;
                    matvar = Mat_VarCreate("doppler_grid_narrow_min", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &aux, 0);
                    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                    Mat_VarFree(matvar);
                }

            /* Dump Spoofing Detection Acquisition Parameters */
            matvar = Mat_VarCreate("d_spoofer_detected", MAT_C_INT32, MAT_T_INT32, 1, dims.data(), &d_spoofer_detected, 0);
            Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
            Mat_VarFree(matvar);

                // std::cout << "\n================== DUMP SPOOFING "<< "SV " << d_gnss_synchro->PRN << " ==================";
                // std::cout << "\nSpoofing Detected: " << d_spoofer_detected;
                // std::cout << "\nSpoofing Present: " << d_spoofer_present;
                // std::cout << "\n=====================================================================\n";

            if (d_spoofer_detected)
            {
                matvar = Mat_VarCreate("pre_correction_doppler_hz", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &d_p_doppler_hz, 0);
                Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                Mat_VarFree(matvar);

                matvar = Mat_VarCreate("pre_correction_delay_samples", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &d_p_delay_samples, 0);
                Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                Mat_VarFree(matvar);

                matvar = Mat_VarCreate("pre_correction_input_power", MAT_C_SINGLE, MAT_T_SINGLE, 1, dims.data(), &d_p_input_power, 0);
                Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                Mat_VarFree(matvar);

                matvar = Mat_VarCreate("pre_correction_sample_counter", MAT_C_UINT64, MAT_T_UINT64, 1, dims.data(), &d_p_sample_counter, 0);
                Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                Mat_VarFree(matvar);

                matvar = Mat_VarCreate("pre_correction_num_dwells", MAT_C_INT32, MAT_T_INT32, 1, dims.data(), &d_p_nci_counter, 0);
                Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                Mat_VarFree(matvar);

                std::array<size_t, 2> dims{static_cast<size_t>(effective_fft_size), static_cast<size_t>(d_num_doppler_bins)};
                matvar_t* matvar = Mat_VarCreate("pre_correction_acq_grid", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims.data(), d_p_grid.memptr(), 0);
                Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                Mat_VarFree(matvar);

                if (d_acq_parameters.make_2_steps)
                {
                    dims[0] = static_cast<size_t>(effective_fft_size);
                    dims[1] = static_cast<size_t>(d_num_doppler_bins_step2);
                    matvar = Mat_VarCreate("pre_correction_acq_grid_narrow", MAT_C_SINGLE, MAT_T_SINGLE, 2, dims.data(), d_p_narrow_grid.memptr(), 0);
                    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);  // or MAT_COMPRESSION_NONE
                    Mat_VarFree(matvar);
                }


            }

            Mat_Close(matfp);
        }
}


float semperfi_pcps_acquisition::max_to_input_power_statistic(uint32_t& indext, int32_t& doppler, uint32_t num_doppler_bins, int32_t doppler_max, int32_t doppler_step)
{    
    float grid_maximum = 0.0;
    uint32_t index_doppler = 0U;
    uint32_t tmp_intex_t = 0U;
    uint32_t index_time = 0U;
    double code_phase = 0;
    double temp_code_phase = 0;
    const int32_t effective_fft_size = (d_acq_parameters.bit_transition_flag ? d_fft_size / 2 : d_fft_size);
    int peak_index = 0;
    int32_t aux_doppler = 0;

    calculate_threshold();

    struct Peak{
        int code_phase;
        int doppler;
        uint32_t indext;
        float mag;
    };

    std::map<float, Peak> peaks;
    
    // Find the correlation peak and the carrier frequency
    for (uint32_t i = 0; i < num_doppler_bins; i++)
    {
        volk_gnsssdr_32f_index_max_32u(&tmp_intex_t, d_magnitude_grid[i].data(), effective_fft_size);

        if (d_magnitude_grid[i][tmp_intex_t] > grid_maximum)
            {
                grid_maximum = d_magnitude_grid[i][tmp_intex_t];
                index_doppler = i;
                index_time = tmp_intex_t;

                code_phase = static_cast<double>(std::fmod(static_cast<float>(tmp_intex_t), d_acq_parameters.samples_per_code));
                
                d_test_statistics = grid_maximum / d_input_power;

                // Search bins with atleast 1 peak above threshold
                

                if (d_test_statistics > d_threshold)
                    {
                        p1d::Persistence1D p;
                        std::vector<float> dp;
                        for(unsigned int k = 0; k < effective_fft_size; k++)
                        {
                            dp.push_back(d_magnitude_grid[i][k]);
                        }
                        p.RunPersistence(dp);
                        std::vector <p1d::TPairedExtrema> Extrema;
                        p.GetPairedExtrema(Extrema, 0)l;

                        for (std::vector< p1d::TPairedExtrema >::iterator it = Extrema.begin(); it != Extrema.end(); it++)
                            {
                                int k = (*it).MaxIndex;
                                temp_code_phase = static_cast<double>(std::fmod(static_cast<float>(k), d_acq_parameters.samples_per_code));
                                if (temp_code_phase > (code_phase - d_acq_parameters.sep_min) && temp_code_phase < (code_phase + d_acq_parameters.sep_min))
                                {
                                    continue;
                                }
                                Peak peak;
                                peak.mag = d_magnitude_grid[i][k];
                                peak.doppler = -static_cast<int32_t>(doppler_max) + d_doppler_center + doppler_step * static_cast<int32_t>(i);

                                peak.code_phase = static_cast<double>(std::fmod(static_cast<float>(k), d_acq_parameters.samples_per_code));
                                peak.indext = k;
                                peaks[peak.mag] = peak;
                                //DLOG(INFO) << "APT: " << peak.mag << ", " << peak.code_phase << ", " << peak.doppler;
                            }
                    }
            }
    }
    
    indext = index_time;

    if (!d_step_two)
        {
            const auto index_opp = (index_doppler + d_num_doppler_bins / 2) % d_num_doppler_bins;
            d_input_power = static_cast<float>(std::accumulate(d_magnitude_grid[index_opp].data(), d_magnitude_grid[index_opp].data() + effective_fft_size, static_cast<float>(0.0)) / effective_fft_size / 2.0 / d_num_noncoherent_integrations_counter);
            doppler = -static_cast<int32_t>(doppler_max) + d_doppler_center + doppler_step * static_cast<int32_t>(index_doppler);    
        }
    else
        {
            doppler = static_cast<int32_t>(d_doppler_center_step_two + (static_cast<float>(index_doppler) - static_cast<float>(floor(d_num_doppler_bins_step2 / 2.0))) * d_acq_parameters.doppler_step2);
        }
  
    if (d_acq_parameters.spoofing_detection)
    {   
       

        std::map<float, Peak>::reverse_iterator rit;
        std::map<float, Peak>::reverse_iterator rit2;

        std::map<float, Peak> step1_peaks;

        for (rit = peaks.rbegin(); rit!=peaks.rend(); ++rit)
        {   
            if (d_threshold / 2 <= rit->second.mag / d_input_power)
            {
                step1_peaks[rit->second.mag] = rit -> second;
                d_spoofer_present = true;
                break;
            }
        }

        if (d_spoofer_present)
        {   
            double aux_code_phase = (--peaks.end())->second.code_phase;
            //std::cout << "\n ===== " << aux_code_phase; 
            index_doppler = (--peaks.end())->second.doppler;
            
            if (!d_step_two)
            {
                const auto index_opp = (index_doppler + d_num_doppler_bins / 2) % d_num_doppler_bins;
                d_input_power = static_cast<float>(std::accumulate(d_magnitude_grid[index_opp].data(), d_magnitude_grid[index_opp].data() + effective_fft_size, static_cast<float>(0.0)) / effective_fft_size / 2.0 / d_num_noncoherent_integrations_counter);
                aux_doppler = -static_cast<int32_t>(doppler_max) + d_doppler_center + doppler_step * static_cast<int32_t>(index_doppler);
            }
            else
            {
                aux_doppler = static_cast<int32_t>(d_doppler_center_step_two + (static_cast<float>(index_doppler) - static_cast<float>(floor(d_num_doppler_bins_step2 / 2.0))) * d_acq_parameters.doppler_step2);
            }
            // Check correct recovery by comparing old and new legitimate and adversarial code delay
            if (d_itr == 0)
            {
                d_legit_code_delay = aux_code_phase;
                d_adv_code_delay = code_phase;
                LOG(INFO) << "ACQDEBUG: LCP " << d_legit_code_delay << " ADVCP " <<code_phase;
            }

            if (abs(aux_code_phase - code_phase) >= d_acq_parameters.sep_min && abs(aux_code_phase - code_phase) <= d_acq_parameters.sep_max )
            {
                d_amp_est = (sqrt(grid_maximum) / effective_fft_size) / effective_fft_size;
                if (d_acq_parameters.verbose)
                {
                    std::cout << "\n============= POTENTIAL SPOOFING DETECTED "<< "SV " << d_gnss_synchro->PRN << " ITR: " << d_global_itr << " =============";
                    std::cout << "\nAdversarial code delay: " << code_phase << " Doppler: " << doppler;
                    std::cout << "\nLegitimate code delay: " << aux_code_phase << " Doppler: " << aux_doppler;
                    std::cout << "\nDifference: " << abs(aux_code_phase - code_phase);
                    std::cout << "\nSemperFi Iteration: " << d_global_itr;
                    std::cout << "\nCoeff: " << grid_maximum;
                    std::cout << "\nAmp_est: " << d_amp_est;
                    std::cout << "\n=====================================================================\n";
                }
                // Copy Pre-Acq Parameters for Dump
                d_spoofer_detected = true;
                d_code_phase = code_phase;
                d_perform_sic = true;

                if (d_itr == 0)
                {
                    //d_acq_parameters.sep_min = abs(aux_code_phase - code_phase);
                    d_p_doppler_hz = static_cast<float>(d_gnss_synchro->Acq_doppler_hz);
                    d_p_delay_samples = static_cast<float>(code_phase);
                    d_p_input_power = d_input_power;
                    d_p_sample_counter = d_sample_counter;
                    d_p_nci_counter = d_num_noncoherent_integrations_counter;
                    d_p_grid = d_grid;
                    d_p_narrow_grid = d_narrow_grid;
                }
            }
            else
            {
                //std::cout << "\n754 CP - " << abs(aux_code_phase - code_phase);
                if (d_itr == 0)
                {
                    d_spoofer_present = false;
                    d_restart_sic = false;
                }
                else{
                    d_restart_sic = true;
                }
            }
            
            // Temporarily declare no spoofing when auxiliary code delay == adversarial code delay - After cancellation, the weaker peak will be adversarial.
            if (code_phase == d_legit_code_delay)
            {
                d_spoofer_present = false;
                d_repeat_acq = false;
                d_recovered = true;
                //dump_time();
            }
            else
            {
                d_recovered = false;
            }

            if (d_acq_parameters.recovery || d_global_itr >= 30)
            {
                // Pass recovered parameters to Tracking
                indext = (--peaks.end())->second.indext;
                doppler = aux_doppler;
                d_spoofer_present = false;
                d_repeat_acq = false;
                d_recovered = true;
                //std::cout << "\n[!] Unable to attenuate. Passing on tracking params without cancellation";
            }
            d_itr++;
            //std::cout << "\n Global itr " << d_global_itr;
        }
       
    }
    return grid_maximum / d_input_power;
}


float semperfi_pcps_acquisition::first_vs_second_peak_statistic(uint32_t& indext, int32_t& doppler, uint32_t num_doppler_bins, int32_t doppler_max, int32_t doppler_step)
{
    // Look for correlation peaks in the results
    // Find the highest peak and compare it to the second highest peak
    // The second peak is chosen not closer than 1 chip to the highest peak

    float firstPeak = 0.0;
    uint32_t index_doppler = 0U;
    uint32_t tmp_intex_t = 0U;
    uint32_t index_time = 0U;

    // Find the correlation peak and the carrier frequency
    for (uint32_t i = 0; i < num_doppler_bins; i++)
        {
            volk_gnsssdr_32f_index_max_32u(&tmp_intex_t, d_magnitude_grid[i].data(), d_fft_size);
            if (d_magnitude_grid[i][tmp_intex_t] > firstPeak)
                {
                    firstPeak = d_magnitude_grid[i][tmp_intex_t];
                    index_doppler = i;
                    index_time = tmp_intex_t;
                    
                }
        }
    indext = index_time;


    if (!d_step_two)
        {
            doppler = -static_cast<int32_t>(doppler_max) + d_doppler_center + doppler_step * static_cast<int32_t>(index_doppler);
        }
    else
        {
            doppler = static_cast<int32_t>(d_doppler_center_step_two + (static_cast<float>(index_doppler) - static_cast<float>(floor(d_num_doppler_bins_step2 / 2.0))) * d_acq_parameters.doppler_step2);
        }

    // Find 1 chip wide code phase exclude range around the peak
    int32_t excludeRangeIndex1 = index_time - d_samplesPerChip;
    int32_t excludeRangeIndex2 = index_time + d_samplesPerChip;

    // Correct code phase exclude range if the range includes array boundaries
    if (excludeRangeIndex1 < 0)
        {
            excludeRangeIndex1 = d_fft_size + excludeRangeIndex1;
        }
    else if (excludeRangeIndex2 >= static_cast<int32_t>(d_fft_size))
        {
            excludeRangeIndex2 = excludeRangeIndex2 - d_fft_size;
        }

    int32_t idx = excludeRangeIndex1;
    memcpy(d_tmp_buffer.data(), d_magnitude_grid[index_doppler].data(), d_fft_size);
    do
        {
            d_tmp_buffer[idx] = 0.0;
            idx++;
            if (idx == static_cast<int32_t>(d_fft_size))
                {
                    idx = 0;
                }
        }
    while (idx != excludeRangeIndex2);

    // Find the second highest correlation peak in the same freq. bin ---
    volk_gnsssdr_32f_index_max_32u(&tmp_intex_t, d_tmp_buffer.data(), d_fft_size);
    const float secondPeak = d_tmp_buffer[tmp_intex_t];

    // Compute the test statistics and compare to the threshold
    return firstPeak / secondPeak;
}


void semperfi_pcps_acquisition::acquisition_core(uint64_t samp_count)
{
    gr::thread::scoped_lock lk(d_setlock);

    // Initialize acquisition algorithm
    int32_t doppler = 0;
    uint32_t indext = 0U;
    const int32_t effective_fft_size = (d_acq_parameters.bit_transition_flag ? d_fft_size / 2 : d_fft_size);
    if (d_cshort)
        {
            volk_gnsssdr_16ic_convert_32fc(d_data_buffer.data(), d_data_buffer_sc.data(), d_consumed_samples);
        }
    memcpy(d_input_signal.data(), d_data_buffer.data(), d_consumed_samples * sizeof(gr_complex));
    if (d_fft_size > d_consumed_samples)
        {
            for (uint32_t i = d_consumed_samples; i < d_fft_size; i++)
                {
                    d_input_signal[i] = gr_complex(0.0, 0.0);
                }
        }
    const gr_complex* in = d_input_signal.data();  // Get the input samples pointer

    d_mag = 0.0;
    d_num_noncoherent_integrations_counter++;

    DLOG(INFO) << "Channel: " << d_channel
               << " , doing acquisition of satellite: " << d_gnss_synchro->System << " " << d_gnss_synchro->PRN
               << " ,sample stamp: " << samp_count << ", threshold: "
               << d_threshold << ", doppler_max: " << d_acq_parameters.doppler_max
               << ", doppler_step: " << d_doppler_step
               << ", use_CFAR_algorithm_flag: " << (d_use_CFAR_algorithm_flag ? "true" : "false");

    lk.unlock();

    // Doppler frequency grid loop
    if (!d_step_two)
        {
            for (uint32_t doppler_index = 0; doppler_index < d_num_doppler_bins; doppler_index++)
                {
                    // Remove Doppler
                    volk_32fc_x2_multiply_32fc(d_fft_if->get_inbuf(), in, d_grid_doppler_wipeoffs[doppler_index].data(), d_fft_size);

                    // Perform the FFT-based convolution  (parallel time search)
                    // Compute the FFT of the carrier wiped--off incoming signal
                    d_fft_if->execute();

                    // Multiply carrier wiped--off, Fourier transformed incoming signal with the local FFT'd code reference
                    volk_32fc_x2_multiply_32fc(d_ifft->get_inbuf(), d_fft_if->get_outbuf(), d_fft_codes.data(), d_fft_size);

                    // Compute the inverse FFT
                    d_ifft->execute();

                    // Compute squared magnitude (and accumulate in case of non-coherent integration)
                    const size_t offset = (d_acq_parameters.bit_transition_flag ? effective_fft_size : 0);
                    if (d_num_noncoherent_integrations_counter == 1)
                        {
                            volk_32fc_magnitude_squared_32f(d_magnitude_grid[doppler_index].data(), d_ifft->get_outbuf() + offset, effective_fft_size);
                        }
                    else
                        {
                            volk_32fc_magnitude_squared_32f(d_tmp_buffer.data(), d_ifft->get_outbuf() + offset, effective_fft_size);
                            volk_32f_x2_add_32f(d_magnitude_grid[doppler_index].data(), d_magnitude_grid[doppler_index].data(), d_tmp_buffer.data(), effective_fft_size);
                        }
                    // Record results to file if required
                    if (d_dump)
                        {
                            memcpy(d_grid.colptr(doppler_index), d_magnitude_grid[doppler_index].data(), sizeof(float) * effective_fft_size);
                        }
                }

            // Compute the test statistic
            if (d_use_CFAR_algorithm_flag)
                {
                    d_test_statistics = max_to_input_power_statistic(indext, doppler, d_num_doppler_bins, d_acq_parameters.doppler_max, d_doppler_step);
                }
            else
                {
                    d_test_statistics = first_vs_second_peak_statistic(indext, doppler, d_num_doppler_bins, d_acq_parameters.doppler_max, d_doppler_step);
                }
            if (d_acq_parameters.use_automatic_resampler)
                {
                    // take into account the acquisition resampler ratio
                    d_gnss_synchro->Acq_delay_samples = static_cast<double>(std::fmod(static_cast<float>(indext), d_acq_parameters.samples_per_code)) * d_acq_parameters.resampler_ratio;
                    d_gnss_synchro->Acq_delay_samples -= static_cast<double>(d_acq_parameters.resampler_latency_samples);  // account the resampler filter latency
                    d_gnss_synchro->Acq_doppler_hz = static_cast<double>(doppler);
                    d_gnss_synchro->Acq_samplestamp_samples = rint(static_cast<double>(samp_count) * d_acq_parameters.resampler_ratio);
                }
            else
                {
                    d_gnss_synchro->Acq_delay_samples = static_cast<double>(std::fmod(static_cast<float>(indext), d_acq_parameters.samples_per_code));
                    d_gnss_synchro->Acq_doppler_hz = static_cast<double>(doppler);
                    d_gnss_synchro->Acq_samplestamp_samples = samp_count;
                }
        }
    else
        {
            for (uint32_t doppler_index = 0; doppler_index < d_num_doppler_bins_step2; doppler_index++)
                {
                    volk_32fc_x2_multiply_32fc(d_fft_if->get_inbuf(), in, d_grid_doppler_wipeoffs_step_two[doppler_index].data(), d_fft_size);

                    // Perform the FFT-based convolution  (parallel time search)
                    // Compute the FFT of the carrier wiped--off incoming signal
                    d_fft_if->execute();

                    // Multiply carrier wiped--off, Fourier transformed incoming signal
                    // with the local FFT'd code reference using SIMD operations with VOLK library
                    volk_32fc_x2_multiply_32fc(d_ifft->get_inbuf(), d_fft_if->get_outbuf(), d_fft_codes.data(), d_fft_size);

                    // compute the inverse FFT
                    d_ifft->execute();

                    const size_t offset = (d_acq_parameters.bit_transition_flag ? effective_fft_size : 0);
                    if (d_num_noncoherent_integrations_counter == 1)
                        {
                            volk_32fc_magnitude_squared_32f(d_magnitude_grid[doppler_index].data(), d_ifft->get_outbuf() + offset, effective_fft_size);
                        }
                    else
                        {
                            volk_32fc_magnitude_squared_32f(d_tmp_buffer.data(), d_ifft->get_outbuf() + offset, effective_fft_size);
                            volk_32f_x2_add_32f(d_magnitude_grid[doppler_index].data(), d_magnitude_grid[doppler_index].data(), d_tmp_buffer.data(), effective_fft_size);
                        }
                    // Record results to file if required
                    if (d_dump and d_channel == d_dump_channel)
                        {
                            memcpy(d_narrow_grid.colptr(doppler_index), d_magnitude_grid[doppler_index].data(), sizeof(float) * effective_fft_size);
                        }
                }
            // Compute the test statistic
            if (d_use_CFAR_algorithm_flag)
                {
                    d_test_statistics = max_to_input_power_statistic(indext, doppler, d_num_doppler_bins_step2, static_cast<int32_t>(d_doppler_center_step_two - (static_cast<float>(d_num_doppler_bins_step2) / 2.0) * d_acq_parameters.doppler_step2), d_acq_parameters.doppler_step2);
                }
            else
                {
                    d_test_statistics = first_vs_second_peak_statistic(indext, doppler, d_num_doppler_bins_step2, static_cast<int32_t>(d_doppler_center_step_two - (static_cast<float>(d_num_doppler_bins_step2) / 2.0) * d_acq_parameters.doppler_step2), d_acq_parameters.doppler_step2);
                }

            if (d_acq_parameters.use_automatic_resampler)
                {
                    // take into account the acquisition resampler ratio
                    d_gnss_synchro->Acq_delay_samples = static_cast<double>(std::fmod(static_cast<float>(indext), d_acq_parameters.samples_per_code)) * d_acq_parameters.resampler_ratio;
                    d_gnss_synchro->Acq_delay_samples -= static_cast<double>(d_acq_parameters.resampler_latency_samples);  // account the resampler filter latency
                    d_gnss_synchro->Acq_doppler_hz = static_cast<double>(doppler);
                    d_gnss_synchro->Acq_samplestamp_samples = rint(static_cast<double>(samp_count) * d_acq_parameters.resampler_ratio);
                    d_gnss_synchro->Acq_doppler_step = d_acq_parameters.doppler_step2;
                }
            else
                {
                    d_gnss_synchro->Acq_delay_samples = static_cast<double>(std::fmod(static_cast<float>(indext), d_acq_parameters.samples_per_code));
                    d_gnss_synchro->Acq_doppler_hz = static_cast<double>(doppler);
                    d_gnss_synchro->Acq_samplestamp_samples = samp_count;
                    d_gnss_synchro->Acq_doppler_step = d_acq_parameters.doppler_step2;
                }
        }

    lk.lock();

    if (!d_recovered && d_spoofer_present)
    {
        // If recovery not successful, check for max iterations.
        if ((d_global_itr % d_acq_parameters.max_itr == 0 && d_global_itr > 0) || (d_restart_sic))//((d_itr %  d_acq_parameters.max_itr == 0) && d_itr != 0)
        {
            auto elapsed = std::chrono::high_resolution_clock::now() - semperfi_start;
            //int time = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
            //std::cout << "\nMax iterations. Restarting SIC for SV " << d_gnss_synchro->PRN << " Time: " << time / 1e6;
            d_test_statistics = 0;
            d_perform_sic = false;
            d_num_noncoherent_integrations_counter = d_acq_parameters.max_dwells;
            d_restart_sic = false;
        }
        else
        {
            d_perform_sic = true;
            return;
        }
    }
    

    if (!d_acq_parameters.bit_transition_flag)
        {
            if (d_test_statistics > d_threshold)
                {
                    d_active = false;
                    if (d_acq_parameters.make_2_steps)
                        {
                            if (d_step_two)
                                {
                                    send_positive_acquisition();
                                    d_step_two = false;
                                    d_acq_state = 0;  // Positive acquisition
                                }
                            else
                                {
                                    d_step_two = true;  // Clear input buffer and make small grid acquisition
                                    d_num_noncoherent_integrations_counter = 0;
                                    calculate_threshold();
                                    d_positive_acq = 0;
                                    d_acq_state = 0;

                                    if (d_spoofer_present)
                                    {
                                        d_perform_sic = true;
                                    }
                                   
                                }
                            calculate_threshold();
                        }
                    else
                        {
                            send_positive_acquisition();
                            d_acq_state = 0;  // Positive acquisition
                        }
                }
            else
                {
                    d_buffer_count = 0;
                    d_acq_state = 1;
                }

            if (d_num_noncoherent_integrations_counter == d_acq_parameters.max_dwells)
                {
                    if (d_acq_state != 0)
                        {
                            send_negative_acquisition();
                            //std::cout << "\n1073 Neg acq ";
                        }
                    d_acq_state = 0;
                    d_active = false;
                    const bool was_step_two = d_step_two;
                    d_step_two = false;
                    if (was_step_two)
                        {
                            calculate_threshold();
                        }
                }
        }
    else
        {
            d_active = false;
            if (d_test_statistics > d_threshold)
                {
                    if (d_acq_parameters.make_2_steps)
                        {
                            if (d_step_two)
                                {
                                    send_positive_acquisition();
                                    d_step_two = false;
                                    d_acq_state = 0;  // Positive acquisition
                                }
                            else
                                {
                                    d_step_two = true;  // Clear input buffer and make small grid acquisition
                                    d_num_noncoherent_integrations_counter = 0U;
                                    calculate_threshold();
                                    d_acq_state = 0;
                                }
                            calculate_threshold();
                        }
                    else
                        {
                            send_positive_acquisition();
                            d_acq_state = 0;  // Positive acquisition
                        }
                }
            else
                {
                    d_acq_state = 0;  // Negative acquisition
                    const bool was_step_two = d_step_two;
                    d_step_two = false;
                    if (was_step_two)
                        {
                            calculate_threshold();
                        }
                    send_negative_acquisition();
                    //std::cout << "\n1127 Neg acq ";
                }
        }
    d_worker_active = false;

    if ((d_num_noncoherent_integrations_counter == d_acq_parameters.max_dwells) or (d_positive_acq == 1))
        {
            // Record results to file if required
            if (d_dump)
                {
                    semperfi_pcps_acquisition::dump_results(effective_fft_size);
                }
            d_num_noncoherent_integrations_counter = 0U;
            d_positive_acq = 0;
        }
}


// Called by gnuradio to enable drivers, etc for i/o devices.
bool semperfi_pcps_acquisition::start()
{
    d_sample_counter = 0ULL;
    calculate_threshold();
    return true;
}


void semperfi_pcps_acquisition::calculate_threshold()
{
    const float pfa = (d_step_two ? d_acq_parameters.pfa2 : d_acq_parameters.pfa);

    if (pfa <= 0.0)
        {
            return;
        }

    const auto effective_fft_size = static_cast<int>(d_acq_parameters.bit_transition_flag ? (d_fft_size / 2) : d_fft_size);
    const int num_doppler_bins = (d_step_two ? d_num_doppler_bins_step2 : d_num_doppler_bins);

    const int num_bins = effective_fft_size * num_doppler_bins;

    d_threshold = static_cast<float>(2.0 * boost::math::gamma_p_inv(2.0 * d_acq_parameters.max_dwells, std::pow(1.0 - pfa, 1.0 / static_cast<float>(num_bins))));
}


int semperfi_pcps_acquisition::general_work(int noutput_items __attribute__((unused)),
    gr_vector_int& ninput_items,
    gr_vector_const_void_star& input_items,
    gr_vector_void_star& output_items __attribute__((unused)))
{
    /*
     * By J.Arribas, L.Esteve and M.Molina
     * Acquisition strategy (Kay Borre book + CFAR threshold):
     * 1. Compute the input signal power estimation
     * 2. Doppler serial search loop
     * 3. Perform the FFT-based circular convolution (parallel time search)
     * 4. Record the maximum peak and the associated synchronization parameters
     * 5. Compute the test statistics and compare to the threshold
     * 6. Declare positive or negative acquisition using a message port
     */
    gr::thread::scoped_lock lk(d_setlock);
    if (!d_active or d_worker_active)
        {
            if (!d_acq_parameters.blocking_on_standby)
                {
                    d_sample_counter += static_cast<uint64_t>(ninput_items[0]);
                    consume_each(ninput_items[0]);
                }
            if (d_step_two)
                {
                    d_doppler_center_step_two = static_cast<float>(d_gnss_synchro->Acq_doppler_hz);
                    update_grid_doppler_wipeoffs_step2();
                    d_acq_state = 0;
                    d_active = true;
                }
            return 0;
        }

    do
    {
        // Acqusition states
        switch (d_acq_state)
        {
            case 0:
            {
                // Restart acquisition variables
                // std::cout << "\n1203 State1";

                d_gnss_synchro->Acq_delay_samples = 0.0;
                d_gnss_synchro->Acq_doppler_hz = 0.0;
                d_gnss_synchro->Acq_samplestamp_samples = 0ULL;
                d_gnss_synchro->Acq_doppler_step = 0U;
                d_mag = 0.0;
                d_acq_state = 1;
                d_buffer_count = 0U;
                d_itr = 0;

                d_adv_code_delay = 0;
                d_legit_code_delay = 0;

                d_codes_generated = false;
                d_recovered = false;
                d_repeat_acq = false;
                d_spoofer_present = false;
                d_spoofer_detected = false;
                //d_global_itr = 0;
                d_phase_set = false;

                if (!d_acq_parameters.blocking_on_standby)
                {
                    d_sample_counter += static_cast<uint64_t>(ninput_items[0]);  // sample counter
                    consume_each(ninput_items[0]);
                }

                // if (d_reset_time)
                // {
                //     semperfi_start = std::chrono::high_resolution_clock::now();
                //     d_reset_time = false;
                // }
                break;
            }
            case 1:
            {
                // std::cout << "\n1242 State2";
                uint32_t buff_increment;
                if (d_cshort)
                    {
                        const auto* in = reinterpret_cast<const lv_16sc_t*>(input_items[0]);  // Get the input samples pointer
                        if ((ninput_items[0] + d_buffer_count) <= d_consumed_samples)
                            {
                                buff_increment = ninput_items[0];
                            }
                        else
                            {
                                buff_increment = d_consumed_samples - d_buffer_count;
                            }
                        memcpy(&d_data_buffer_sc[d_buffer_count], in, sizeof(lv_16sc_t) * buff_increment);
                    }
                else
                    {
                        const auto* in = reinterpret_cast<const gr_complex*>(input_items[0]);  // Get the input samples pointer
                        if ((ninput_items[0] + d_buffer_count) <= d_consumed_samples)
                            {
                                buff_increment = ninput_items[0];
                            }
                        else
                            {
                                buff_increment = d_consumed_samples - d_buffer_count;
                            }
                        memcpy(&d_data_buffer[d_buffer_count], in, sizeof(gr_complex) * buff_increment);
                    }

                // If buffer will be full in next iteration
                if (d_buffer_count >= d_consumed_samples)
                    {
                        //semperfi_start = std::chrono::high_resolution_clock::now();
                        d_acq_state = 2;
                    }

                d_buffer_count += buff_increment;
                d_sample_counter += static_cast<uint64_t>(buff_increment);
                consume_each(buff_increment);
                d_repeat_acq = false;
                break;
            }
            case 2:
            // Perform acquisition and spoofing detection
            {
                //std::cout << "\n1287 State3";
                // Copy the data to the core and let it know that new data is available
                if (d_acq_parameters.blocking)
                {
                    lk.unlock();
                    d_buffer_count = d_acq_samples_count;
                    acquisition_core(d_sample_counter);
                    d_acq_samples_count = d_buffer_count;
                }
                else
                {
                    gr::thread::thread d_worker(&semperfi_pcps_acquisition::acquisition_core, this, d_sample_counter);
                    d_worker_active = true;
                }

                if (d_perform_sic)
                {
                    // Set ACQ block state to cancellation and recovery
                    d_acq_state = 3;
                    d_spoofer_present = false;
                }
                
                d_buffer_count = 0U;
                d_repeat_acq = false;
                consume_each(0);
                d_global_itr++;
                break;
            }

            case 3:
            // Cancellation and Recovery state
            {
                d_acq_samples_count = d_data_buffer.size();
                
                // Signal generator stuff ------------------------------------------------------------------------
                if (!d_codes_generated)
                {
                    d_codes_generated = true;
                    signal_gen_init();
                    generate_codes();
                }
                
                generate_signal(gr_complex(1, 0));

                // Signal generator stuff ------------------------------------------------------------------------

                // Recovery stuff ------------------------------------------------------------------------
                
                // Add the generated signal and the incoming signal
                unsigned int alignment = volk_get_alignment();

                const auto* in = reinterpret_cast<const gr_complex*>(&d_data_buffer[0]);
                const auto* rec = reinterpret_cast<const gr_complex*>(&d_recovery_signal_buff[0]);

                std::map<float, double> amps;

                if (!d_phase_set)
                {
                    for (int degree = 0; degree < 181; degree++)
                    {
                        lv_32fc_t* shifted_rec = (lv_32fc_t*)volk_malloc(sizeof(gr_complex)*d_acq_samples_count, alignment);
                        lv_32fc_t* out = (lv_32fc_t*)volk_malloc(sizeof(gr_complex)*d_acq_samples_count, alignment);

                        double rad = (degree * M_PI) / 180;

                        float *mean = (float *)volk_malloc(sizeof(float), alignment);
                        float *stddev = (float *)volk_malloc(sizeof(float), alignment);
                        float* magnitude = (float*)volk_malloc(sizeof(float)*d_acq_samples_count, alignment);

                        gr_complex multiplier = gr_complex(cos(rad), sin(rad));

                        volk_32fc_s32fc_multiply_32fc(shifted_rec, rec, multiplier, d_acq_samples_count);

                        volk_32fc_x2_add_32fc(out, in, shifted_rec, d_acq_samples_count);

                        volk_32fc_magnitude_32f(magnitude, out, d_acq_samples_count);

                        volk_32f_stddev_and_mean_32f_x2(stddev, mean, magnitude, d_acq_samples_count);

                        amps[*mean] = rad;
                        volk_free(out);
                        volk_free(shifted_rec);
                        volk_free(mean);
                        volk_free(stddev);
                        volk_free(magnitude);
                    }
                    std::map<float, double>::iterator it;

                    it = amps.begin();

                    d_phase = it->second;
                    d_phase_set = true;
                }

                lv_32fc_t* shifted_rec = (lv_32fc_t*)volk_malloc(sizeof(gr_complex)*d_acq_samples_count, alignment);
                lv_32fc_t* out = (lv_32fc_t*)volk_malloc(sizeof(gr_complex)*d_acq_samples_count, alignment);

                gr_complex multiplier = gr_complex(cos(d_phase), sin(d_phase));

                volk_32fc_s32fc_multiply_32fc(shifted_rec, rec, multiplier, d_acq_samples_count);

                volk_32fc_x2_add_32fc(out, in, shifted_rec, d_acq_samples_count);

                memcpy(&d_data_buffer[0], out, sizeof(gr_complex) * d_acq_samples_count);

                // Recovery stuff ------------------------------------------------------------------------
               
                // Set acq state for re-acquisition

                // Reset d_magnitude grid
                d_magnitude_grid = volk_gnsssdr::vector<volk_gnsssdr::vector<float>>(d_num_doppler_bins, volk_gnsssdr::vector<float>(d_fft_size));
  
                for (uint32_t doppler_index = 0; doppler_index < d_num_doppler_bins; doppler_index++)
                {
                    std::fill(d_magnitude_grid[doppler_index].begin(), d_magnitude_grid[doppler_index].end(), 0.0);
                }

                update_grid_doppler_wipeoffs();

                d_acq_state = 2;
                d_repeat_acq = true;
                
            }
        }
       
    }
    while(d_repeat_acq);

    return 0;
}


//////////////////////////////////////////////////// Recovery signal generator function

void semperfi_pcps_acquisition::signal_gen_init()
{
    //std::cout << "1473 sig gen int";

    num_sats_ = 1;

    prn = d_gnss_synchro->PRN;
    fs_in_ = d_acq_parameters.fs_in;

    complex_phase_.reserve(d_acq_samples_count);

    // True if Galileo satellites are present
    bool galileo_signal = false;

    start_phase_rad = 0;
    current_data_bit_int_ = 0;

    ms_counter_ = 0;

    samples_per_code_ = round(static_cast<float>(fs_in_) / (GPS_L1_CA_CODE_RATE_CPS / GPS_L1_CA_CODE_LENGTH_CHIPS));

    num_of_codes_per_vector_ = (galileo_signal ? 4 * static_cast<int>(GALILEO_E1_C_SECONDARY_CODE_LENGTH) : 1);
    data_bit_duration_ms_ = (1e3 / GPS_CA_TELEMETRY_RATE_BITS_SECOND);
}


void semperfi_pcps_acquisition::generate_codes()
{
    //std::cout << "1499 code gen";
    auto delay_samples = d_code_phase;
    int delay_chips = static_cast<int>(delay_samples * static_cast<int>(GPS_L1_CA_CODE_LENGTH_CHIPS)) / samples_per_code_;

    sampled_code_data_ = std::vector<gr_complex>(std::vector<gr_complex>(d_acq_samples_count));

    std::array<gr_complex, 64000> code{};

    // Generate one code-period of 1C signal
    gps_l1_ca_code_gen_complex_sampled(code, prn, fs_in_, static_cast<int>(GPS_L1_CA_CODE_LENGTH_CHIPS) - delay_chips);

    // Concatenate "num_of_codes_per_vector_" codes
    for (unsigned int i = 0; i < num_of_codes_per_vector_; i++)
        {
            memcpy(&(sampled_code_data_[i * samples_per_code_]),
                code.data(), sizeof(gr_complex) * samples_per_code_);
        }
}


void semperfi_pcps_acquisition::generate_signal(gr_complex data_bit)
{
    // Set delay samples and doppler
    auto doppler_Hz = d_gnss_synchro->Acq_doppler_hz;

    // ONLY FOR DEBUGGING - REMOVE AFTER IMPLEMENTING AMP ESTIMATION
    if (d_acq_parameters.amp != 0)
    {
        d_amp_est = d_acq_parameters.amp;
    }

    unsigned int alignment = volk_get_alignment();

    lv_32fc_t* out = (lv_32fc_t*)volk_malloc(sizeof(gr_complex)*d_acq_samples_count, alignment);

    float phase_step_rad = -static_cast<float>(TWO_PI) * doppler_Hz / static_cast<float>(fs_in_);

    std::array<float, 1> _phase{};
    _phase[0] = -start_phase_rad;
    volk_gnsssdr_s32f_sincos_32fc(complex_phase_.data(), -phase_step_rad, _phase.data(), d_acq_samples_count);
    start_phase_rad += static_cast<float>(d_acq_samples_count) * phase_step_rad;

    unsigned int out_idx = 0;
    unsigned int i = 0;
    unsigned int k = 0;

    for (out_idx = 0; out_idx < samples_per_code_; out_idx++)
    {
        out[out_idx] = gr_complex(0.0, 0.0);
    }

    out_idx = 0;
    for (i = 0; i < num_of_codes_per_vector_; i++)
    {
        for (k = 0; k < samples_per_code_; k++)
            {
                out[out_idx] = sampled_code_data_[out_idx] * data_bit * complex_phase_[out_idx];
                out_idx++;
            }

        ms_counter_ = (ms_counter_ + static_cast<int>(round(1e3 * GPS_L1_CA_CODE_PERIOD_S))) % data_bit_duration_ms_;
    }

    lv_32fc_t* temp_out = (lv_32fc_t*)volk_malloc(sizeof(gr_complex)*d_acq_samples_count, alignment);

    gr_complex amp_est = gr_complex(-d_amp_est, 0);
    volk_32fc_s32fc_multiply_32fc(temp_out, out, amp_est, d_acq_samples_count);
    memcpy(&d_recovery_signal_buff[0], temp_out, sizeof(gr_complex) * d_acq_samples_count);
}