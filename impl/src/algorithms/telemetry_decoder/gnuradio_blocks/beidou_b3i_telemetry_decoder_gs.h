/*!
 * \file beidou_b3i_telemetry_decoder_gs.h
 * \brief Implementation of a BEIDOU B3I DNAV data decoder block
 * \author Damian Miralles, 2019. dmiralles2009(at)gmail.com
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

#ifndef GNSS_SDR_BEIDOU_B3I_TELEMETRY_DECODER_GS_H
#define GNSS_SDR_BEIDOU_B3I_TELEMETRY_DECODER_GS_H

#include "beidou_dnav_navigation_message.h"
#include "gnss_satellite.h"
#include <boost/circular_buffer.hpp>
#include <gnuradio/block.h>  // for block
#include <gnuradio/types.h>  // for gr_vector_const_void_star
#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#if GNURADIO_USES_STD_POINTERS
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

class beidou_b3i_telemetry_decoder_gs;

#if GNURADIO_USES_STD_POINTERS
using beidou_b3i_telemetry_decoder_gs_sptr =
    std::shared_ptr<beidou_b3i_telemetry_decoder_gs>;
#else
using beidou_b3i_telemetry_decoder_gs_sptr =
    boost::shared_ptr<beidou_b3i_telemetry_decoder_gs>;
#endif

beidou_b3i_telemetry_decoder_gs_sptr beidou_b3i_make_telemetry_decoder_gs(
    const Gnss_Satellite &satellite,
    bool dump);

/*!
 * \brief This class implements a block that decodes the BeiDou DNAV data.
 */
class beidou_b3i_telemetry_decoder_gs : public gr::block
{
public:
    ~beidou_b3i_telemetry_decoder_gs();                   //!< Class destructor
    void set_satellite(const Gnss_Satellite &satellite);  //!< Set satellite PRN
    void set_channel(int channel);                        //!< Set receiver's channel
    void reset();

    /*!
     * \brief This is where all signal processing takes place
     */
    int general_work(int noutput_items, gr_vector_int &ninput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items);

private:
    friend beidou_b3i_telemetry_decoder_gs_sptr beidou_b3i_make_telemetry_decoder_gs(
        const Gnss_Satellite &satellite,
        bool dump);

    beidou_b3i_telemetry_decoder_gs(const Gnss_Satellite &satellite, bool dump);

    void decode_subframe(float *symbols);
    void decode_word(int32_t word_counter, const float *enc_word_symbols,
        int32_t *dec_word_symbols);
    void decode_bch15_11_01(const int32_t *bits, std::array<int32_t, 15> &decbits);

    // Preamble decoding
    std::array<int32_t, BEIDOU_DNAV_PREAMBLE_LENGTH_SYMBOLS> d_preamble_samples{};
    std::array<float, BEIDOU_DNAV_PREAMBLE_PERIOD_SYMBOLS> d_subframe_symbols{};

    // Storage for incoming data
    boost::circular_buffer<float> d_symbol_history;

    // Navigation Message variable
    Beidou_Dnav_Navigation_Message d_nav;

    Gnss_Satellite d_satellite;

    std::string d_dump_filename;
    std::ofstream d_dump_file;

    uint64_t d_sample_counter;  // Sample counter as an index (1,2,3,..etc) indicating number of samples processed
    uint64_t d_preamble_index;  // Index of sample number where preamble was found
    uint32_t d_required_symbols;
    uint32_t d_stat;  // Status of decoder

    int32_t d_channel;
    int32_t d_CRC_error_counter;  // Number of failed CRC operations
    int32_t d_symbols_per_preamble;
    int32_t d_samples_per_preamble;
    int32_t d_preamble_period_samples;

    // Values to populate gnss synchronization structure
    uint64_t d_last_valid_preamble;
    uint32_t d_symbol_duration_ms;
    uint32_t d_TOW_at_Preamble_ms;
    uint32_t d_TOW_at_current_symbol_ms;

    bool flag_SOW_set;       // Indicates when time of week is set
    bool d_flag_frame_sync;  // Indicate when a frame sync is achieved
    bool d_flag_preamble;    // Flag indicating when preamble was found
    bool d_flag_valid_word;
    bool d_sent_tlm_failed_msg;
    bool Flag_valid_word;
    bool d_dump;
};

#endif  // GNSS_SDR_BEIDOU_B3I_TELEMETRY_DECODER_GS_H