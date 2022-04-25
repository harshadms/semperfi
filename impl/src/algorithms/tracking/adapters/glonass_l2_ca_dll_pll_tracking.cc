/*!
 * \file glonass_l2_ca_dll_pll_tracking.cc
 * \brief  Interface of an adapter of a DLL+PLL tracking loop block
 * for Glonass L2 C/A to a TrackingInterface
 * \author Damian Miralles, 2018, dmiralles2009(at)gmail.com *
 *
 * Code DLL + carrier PLL according to the algorithms described in:
 * K.Borre, D.M.Akos, N.Bertelsen, P.Rinder, and S.H.Jensen,
 * A Software-Defined GPS and Galileo Receiver. A Single-Frequency
 * Approach, Birkha user, 2007
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

#include "glonass_l2_ca_dll_pll_tracking.h"
#include "GLONASS_L1_L2_CA.h"
#include "configuration_interface.h"
#include "gnss_sdr_flags.h"
#include <glog/logging.h>


GlonassL2CaDllPllTracking::GlonassL2CaDllPllTracking(
    const ConfigurationInterface* configuration, const std::string& role,
    unsigned int in_streams, unsigned int out_streams) : role_(role), in_streams_(in_streams), out_streams_(out_streams)
{
    DLOG(INFO) << "role " << role;
    // ################# CONFIGURATION PARAMETERS ########################
    const std::string default_item_type("gr_complex");
    std::string item_type = configuration->property(role + ".item_type", default_item_type);
    int fs_in_deprecated = configuration->property("GNSS-SDR.internal_fs_hz", 2048000);
    int fs_in = configuration->property("GNSS-SDR.internal_fs_sps", fs_in_deprecated);
    bool dump = configuration->property(role + ".dump", false);
    float pll_bw_hz = configuration->property(role + ".pll_bw_hz", static_cast<float>(50.0));
    if (FLAGS_pll_bw_hz != 0.0)
        {
            pll_bw_hz = static_cast<float>(FLAGS_pll_bw_hz);
        }
    float dll_bw_hz = configuration->property(role + ".dll_bw_hz", static_cast<float>(2.0));
    if (FLAGS_dll_bw_hz != 0.0)
        {
            dll_bw_hz = static_cast<float>(FLAGS_dll_bw_hz);
        }
    float early_late_space_chips = configuration->property(role + ".early_late_space_chips", static_cast<float>(0.5));
    const std::string default_dump_filename("./track_ch");
    std::string dump_filename = configuration->property(role + ".dump_filename", default_dump_filename);
    const auto vector_length = static_cast<int>(std::round(fs_in / (GLONASS_L2_CA_CODE_RATE_CPS / GLONASS_L2_CA_CODE_LENGTH_CHIPS)));

    // ################# MAKE TRACKING GNURadio object ###################
    if (item_type == "gr_complex")
        {
            item_size_ = sizeof(gr_complex);
            tracking_ = glonass_l2_ca_dll_pll_make_tracking_cc(
                fs_in,
                vector_length,
                dump,
                dump_filename,
                pll_bw_hz,
                dll_bw_hz,
                early_late_space_chips);
        }
    else
        {
            item_size_ = sizeof(gr_complex);
            LOG(WARNING) << item_type << " unknown tracking item type.";
        }
    channel_ = 0;
    DLOG(INFO) << "tracking(" << tracking_->unique_id() << ")";
    if (in_streams_ > 1)
        {
            LOG(ERROR) << "This implementation only supports one input stream";
        }
    if (out_streams_ > 1)
        {
            LOG(ERROR) << "This implementation only supports one output stream";
        }
}


void GlonassL2CaDllPllTracking::stop_tracking()
{
}


void GlonassL2CaDllPllTracking::start_tracking()
{
    tracking_->start_tracking();
}


/*
 * Set tracking channel unique ID
 */
void GlonassL2CaDllPllTracking::set_channel(unsigned int channel)
{
    channel_ = channel;
    tracking_->set_channel(channel);
}


void GlonassL2CaDllPllTracking::set_gnss_synchro(Gnss_Synchro* p_gnss_synchro)
{
    tracking_->set_gnss_synchro(p_gnss_synchro);
}


void GlonassL2CaDllPllTracking::connect(gr::top_block_sptr top_block)
{
    if (top_block)
        { /* top_block is not null */
        };
    // nothing to connect, now the tracking uses gr_sync_decimator
}


void GlonassL2CaDllPllTracking::disconnect(gr::top_block_sptr top_block)
{
    if (top_block)
        { /* top_block is not null */
        };
    // nothing to disconnect, now the tracking uses gr_sync_decimator
}


gr::basic_block_sptr GlonassL2CaDllPllTracking::get_left_block()
{
    return tracking_;
}


gr::basic_block_sptr GlonassL2CaDllPllTracking::get_right_block()
{
    return tracking_;
}