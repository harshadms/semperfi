/*!
 * \file gnss_block_factory_test.cc
 * \brief This class implements a Unit Test for the GNSSBlockFactory class.
 * \authors <ul>
 *          <li> Carlos Aviles, 2010. carlos.avilesr(at)googlemail.com
 *          <li> Luis Esteve, 2012. luis(at)epsilon-formacion.com
 *          </ul>
 *
 * This class test the instantiation of all blocks in gnss_block_factory
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

#include "acquisition_interface.h"
#include "channel.h"
#include "concurrent_queue.h"
#include "gnss_block_factory.h"
#include "gnss_block_interface.h"
#include "gnss_sdr_make_unique.h"
#include "in_memory_configuration.h"
#include "observables_interface.h"
#include "pvt_interface.h"
#include "telemetry_decoder_interface.h"
#include "tracking_interface.h"
#include <gtest/gtest.h>
#include <pmt/pmt.h>
#include <vector>

TEST(GNSSBlockFactoryTest, InstantiateFileSignalSource)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("SignalSource.implementation", "File_Signal_Source");
    std::string path = std::string(TEST_PATH);
    std::string filename = path + "signal_samples/GPS_L1_CA_ID_1_Fs_4Msps_2ms.dat";
    configuration->set_property("SignalSource.filename", filename);
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue = std::make_shared<Concurrent_Queue<pmt::pmt_t>>();
    // Example of a factory as a shared_ptr
    std::shared_ptr<GNSSBlockFactory> factory = std::make_shared<GNSSBlockFactory>();
    // Example of a block as a shared_ptr
    std::shared_ptr<GNSSBlockInterface> signal_source = factory->GetSignalSource(configuration.get(), queue.get());
    EXPECT_STREQ("SignalSource", signal_source->role().c_str());
    EXPECT_STREQ("File_Signal_Source", signal_source->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateWrongSignalSource)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("SignalSource.implementation", "Pepito");
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue = std::make_shared<Concurrent_Queue<pmt::pmt_t>>();
    // Example of a factory as a unique_ptr
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    // Example of a block as a unique_ptr
    std::unique_ptr<GNSSBlockInterface> signal_source = factory->GetSignalSource(configuration.get(), queue.get());
    EXPECT_EQ(nullptr, signal_source);
}


TEST(GNSSBlockFactoryTest, InstantiateSignalConditioner)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("SignalConditioner.implementation", "Signal_Conditioner");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<GNSSBlockInterface> signal_conditioner = factory->GetSignalConditioner(configuration.get());
    EXPECT_STREQ("SignalConditioner", signal_conditioner->role().c_str());
    EXPECT_STREQ("Signal_Conditioner", signal_conditioner->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateFIRFilter)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue = std::make_shared<Concurrent_Queue<pmt::pmt_t>>();

    configuration->set_property("InputFilter.implementation", "Fir_Filter");

    configuration->set_property("InputFilter.number_of_taps", "4");
    configuration->set_property("InputFilter.number_of_bands", "2");

    configuration->set_property("InputFilter.band1_begin", "0.0");
    configuration->set_property("InputFilter.band1_end", "0.45");
    configuration->set_property("InputFilter.band2_begin", "0.55");
    configuration->set_property("InputFilter.band2_end", "1.0");

    configuration->set_property("InputFilter.ampl1_begin", "1.0");
    configuration->set_property("InputFilter.ampl1_end", "1.0");
    configuration->set_property("InputFilter.ampl2_begin", "0.0");
    configuration->set_property("InputFilter.ampl2_end", "0.0");

    configuration->set_property("InputFilter.band1_error", "1.0");
    configuration->set_property("InputFilter.band2_error", "1.0");

    configuration->set_property("InputFilter.filter_type", "bandpass");
    configuration->set_property("InputFilter.grid_density", "16");

    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<GNSSBlockInterface> input_filter = factory->GetBlock(configuration.get(), "InputFilter", 1, 1);

    EXPECT_STREQ("InputFilter", input_filter->role().c_str());
    EXPECT_STREQ("Fir_Filter", input_filter->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateFreqXlatingFIRFilter)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue = std::make_shared<Concurrent_Queue<pmt::pmt_t>>();

    configuration->set_property("InputFilter.implementation", "Freq_Xlating_Fir_Filter");

    configuration->set_property("InputFilter.number_of_taps", "4");
    configuration->set_property("InputFilter.number_of_bands", "2");

    configuration->set_property("InputFilter.band1_begin", "0.0");
    configuration->set_property("InputFilter.band1_end", "0.45");
    configuration->set_property("InputFilter.band2_begin", "0.55");
    configuration->set_property("InputFilter.band2_end", "1.0");

    configuration->set_property("InputFilter.ampl1_begin", "1.0");
    configuration->set_property("InputFilter.ampl1_end", "1.0");
    configuration->set_property("InputFilter.ampl2_begin", "0.0");
    configuration->set_property("InputFilter.ampl2_end", "0.0");

    configuration->set_property("InputFilter.band1_error", "1.0");
    configuration->set_property("InputFilter.band2_error", "1.0");

    configuration->set_property("InputFilter.filter_type", "bandpass");
    configuration->set_property("InputFilter.grid_density", "16");

    configuration->set_property("InputFilter.sampling_frequency", "4000000");
    configuration->set_property("InputFilter.IF", "34000");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<GNSSBlockInterface> input_filter = factory->GetBlock(configuration.get(), "InputFilter", 1, 1);

    EXPECT_STREQ("InputFilter", input_filter->role().c_str());
    EXPECT_STREQ("Freq_Xlating_Fir_Filter", input_filter->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiatePulseBlankingFilter)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue = std::make_shared<Concurrent_Queue<pmt::pmt_t>>();
    configuration->set_property("InputFilter.implementation", "Pulse_Blanking_Filter");

    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<GNSSBlockInterface> input_filter = factory->GetBlock(configuration.get(), "InputFilter", 1, 1);

    EXPECT_STREQ("InputFilter", input_filter->role().c_str());
    EXPECT_STREQ("Pulse_Blanking_Filter", input_filter->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateNotchFilter)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue = std::make_shared<Concurrent_Queue<pmt::pmt_t>>();
    configuration->set_property("InputFilter.implementation", "Notch_Filter");

    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<GNSSBlockInterface> input_filter = factory->GetBlock(configuration.get(), "InputFilter", 1, 1);

    EXPECT_STREQ("InputFilter", input_filter->role().c_str());
    EXPECT_STREQ("Notch_Filter", input_filter->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateNotchFilterLite)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue = std::make_shared<Concurrent_Queue<pmt::pmt_t>>();
    configuration->set_property("InputFilter.implementation", "Notch_Filter_Lite");

    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<GNSSBlockInterface> input_filter = factory->GetBlock(configuration.get(), "InputFilter", 1, 1);

    EXPECT_STREQ("InputFilter", input_filter->role().c_str());
    EXPECT_STREQ("Notch_Filter_Lite", input_filter->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateDirectResampler)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Resampler.implementation", "Direct_Resampler");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<GNSSBlockInterface> resampler = factory->GetBlock(configuration.get(), "Resampler", 1, 1);
    EXPECT_STREQ("Resampler", resampler->role().c_str());
    EXPECT_STREQ("Direct_Resampler", resampler->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateGpsL1CaPcpsAcquisition)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Acquisition.implementation", "GPS_L1_CA_PCPS_Acquisition");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> acq_ = factory->GetBlock(configuration.get(), "Acquisition", 1, 0);
    std::shared_ptr<AcquisitionInterface> acquisition = std::dynamic_pointer_cast<AcquisitionInterface>(acq_);
    EXPECT_STREQ("Acquisition", acquisition->role().c_str());
    EXPECT_STREQ("GPS_L1_CA_PCPS_Acquisition", acquisition->implementation().c_str());
}

TEST(GNSSBlockFactoryTest, InstantiateGpsL1CaPcpsQuickSyncAcquisition)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Acquisition.implementation", "GPS_L1_CA_PCPS_QuickSync_Acquisition");
    std::shared_ptr<GNSSBlockFactory> factory = std::make_shared<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> acq_ = factory->GetBlock(configuration.get(), "Acquisition", 1, 0);
    std::shared_ptr<AcquisitionInterface> acquisition = std::dynamic_pointer_cast<AcquisitionInterface>(acq_);
    EXPECT_STREQ("Acquisition", acquisition->role().c_str());
    EXPECT_STREQ("GPS_L1_CA_PCPS_QuickSync_Acquisition", acquisition->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateGalileoE1PcpsQuickSyncAmbiguousAcquisition)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Acquisition.implementation", "Galileo_E1_PCPS_QuickSync_Ambiguous_Acquisition");
    std::shared_ptr<GNSSBlockFactory> factory = std::make_shared<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> acq_ = factory->GetBlock(configuration.get(), "Acquisition", 1, 0);
    std::shared_ptr<AcquisitionInterface> acquisition = std::dynamic_pointer_cast<AcquisitionInterface>(acq_);
    EXPECT_STREQ("Acquisition", acquisition->role().c_str());
    EXPECT_STREQ("Galileo_E1_PCPS_QuickSync_Ambiguous_Acquisition", acquisition->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateGalileoE1PcpsAmbiguousAcquisition)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Acquisition.implementation", "Galileo_E1_PCPS_Ambiguous_Acquisition");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> acq_ = factory->GetBlock(configuration.get(), "Acquisition", 1, 0);
    std::shared_ptr<AcquisitionInterface> acquisition = std::dynamic_pointer_cast<AcquisitionInterface>(acq_);
    EXPECT_STREQ("Acquisition", acquisition->role().c_str());
    EXPECT_STREQ("Galileo_E1_PCPS_Ambiguous_Acquisition", acquisition->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateGpsL1CaDllPllTracking)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Tracking.implementation", "GPS_L1_CA_DLL_PLL_Tracking");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> trk_ = factory->GetBlock(configuration.get(), "Tracking", 1, 1);
    std::shared_ptr<TrackingInterface> tracking = std::dynamic_pointer_cast<TrackingInterface>(trk_);
    EXPECT_STREQ("Tracking", tracking->role().c_str());
    EXPECT_STREQ("GPS_L1_CA_DLL_PLL_Tracking", tracking->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateGpsL1CaTcpConnectorTracking)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Tracking.implementation", "GPS_L1_CA_TCP_CONNECTOR_Tracking");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> trk_ = factory->GetBlock(configuration.get(), "Tracking", 1, 1);
    std::shared_ptr<TrackingInterface> tracking = std::dynamic_pointer_cast<TrackingInterface>(trk_);
    EXPECT_STREQ("Tracking", tracking->role().c_str());
    EXPECT_STREQ("GPS_L1_CA_TCP_CONNECTOR_Tracking", tracking->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateGalileoE1DllPllVemlTracking)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Tracking.implementation", "Galileo_E1_DLL_PLL_VEML_Tracking");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> trk_ = factory->GetBlock(configuration.get(), "Tracking", 1, 1);
    std::shared_ptr<TrackingInterface> tracking = std::dynamic_pointer_cast<TrackingInterface>(trk_);
    EXPECT_STREQ("Tracking", tracking->role().c_str());
    EXPECT_STREQ("Galileo_E1_DLL_PLL_VEML_Tracking", tracking->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateGpsL1CaTelemetryDecoder)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("TelemetryDecoder.implementation", "GPS_L1_CA_Telemetry_Decoder");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> telemetry_decoder = factory->GetBlock(configuration.get(), "TelemetryDecoder", 1, 1);
    EXPECT_STREQ("TelemetryDecoder", telemetry_decoder->role().c_str());
    EXPECT_STREQ("GPS_L1_CA_Telemetry_Decoder", telemetry_decoder->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateChannels)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Channels_1C.count", "2");
    configuration->set_property("Channels_1E.count", "0");
    configuration->set_property("Channels.in_acquisition", "2");
    configuration->set_property("Tracking_1C.implementation", "GPS_L1_CA_DLL_PLL_Tracking");
    configuration->set_property("TelemetryDecoder_1C.implementation", "GPS_L1_CA_Telemetry_Decoder");
    configuration->set_property("Channel0.item_type", "gr_complex");
    configuration->set_property("Acquisition_1C.implementation", "GPS_L1_CA_PCPS_Acquisition");
    configuration->set_property("Channel1.item_type", "gr_complex");
    std::shared_ptr<Concurrent_Queue<pmt::pmt_t>> queue = std::make_shared<Concurrent_Queue<pmt::pmt_t>>();
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<std::vector<std::unique_ptr<GNSSBlockInterface>>> channels = factory->GetChannels(configuration.get(), queue.get());
    EXPECT_EQ(static_cast<unsigned int>(2), channels->size());
    channels->erase(channels->begin(), channels->end());
}


TEST(GNSSBlockFactoryTest, InstantiateObservables)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Observables.implementation", "Pass_Through");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    auto observables = factory->GetObservables(configuration.get());
    EXPECT_STREQ("Observables", observables->role().c_str());
    EXPECT_STREQ("Pass_Through", observables->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateGpsL1CaObservables)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("Observables.implementation", "Hybrid_Observables");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::unique_ptr<GNSSBlockInterface> observables = factory->GetObservables(configuration.get());
    EXPECT_STREQ("Observables", observables->role().c_str());
    EXPECT_STREQ("Hybrid_Observables", observables->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateRTKLIBPvt)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("PVT.implementation", "RTKLIB_PVT");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> pvt_ = factory->GetPVT(configuration.get());
    std::shared_ptr<PvtInterface> pvt = std::dynamic_pointer_cast<PvtInterface>(pvt_);
    EXPECT_STREQ("PVT", pvt->role().c_str());
    EXPECT_STREQ("RTKLIB_PVT", pvt->implementation().c_str());
}


TEST(GNSSBlockFactoryTest, InstantiateWrongPvt)
{
    std::shared_ptr<InMemoryConfiguration> configuration = std::make_shared<InMemoryConfiguration>();
    configuration->set_property("PVT.implementation", "Pepito");
    std::unique_ptr<GNSSBlockFactory> factory = std::make_unique<GNSSBlockFactory>();
    std::shared_ptr<GNSSBlockInterface> pvt_ = factory->GetPVT(configuration.get());
    std::shared_ptr<PvtInterface> pvt = std::dynamic_pointer_cast<PvtInterface>(pvt_);
    EXPECT_EQ(nullptr, pvt);
}
