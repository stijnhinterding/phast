/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

#include "fake_timetag_device_plugin.h"

#include "ui_fake_device_settings.h"

#include <iostream>

fake_timetag_device_plugin::fake_timetag_device_plugin() :
    comm(new FakeCommunicator(1e-12, 1/50e-9, 1e-6, 0.3, 128))
{
}

fake_timetag_device_plugin::~fake_timetag_device_plugin()
{
    delete this->comm;
}

QString fake_timetag_device_plugin::TimeTagPluginName() const
{
    return "fake-timetag-device plugin";
}

QString fake_timetag_device_plugin::TimeTagPluginDescriptor() const
{
    return "Fake device to ease testing of the time-tag software";
}

QString fake_timetag_device_plugin::TimeTagMenuDescriptor() const
{
    return "Fake device settings";
}

ITimeTaggerCommunicator* fake_timetag_device_plugin::GetComm() const
{
    return this->comm;
}

void fake_timetag_device_plugin::DisableChan(chan_id ID)
{
    this->comm->DisableChan(ID);
}

bool fake_timetag_device_plugin::TimeTagHasSettingsMenu() const
{
    return true;
}

void fake_timetag_device_plugin::ShowTriggerSettingsUI(QWidget *parent, const std::vector<chan_id> &chan_IDs)
{
    ui_fake_device_settings::fake_device_settings settings;
    settings.detection_probability = this->comm->GetDetectionProbability();
    settings.lifetime = 1.0/this->comm->GetDecayRate();
    settings.pulse_period = this->comm->GetPulsePeriod();
    settings.time_unit = this->comm->TimeUnit();
    settings.sync_divider = this->comm->GetSyncDivider(0);
    settings.n_emitters = this->comm->GetNumEmitters();
    settings.sim_speedup_factor = this->comm->GetSimSpeedupFactor();

    ui_fake_device_settings ui(parent, settings, false);

    ui.exec();
    this->comm->SetDecayRate(1.0/ui.Settings().lifetime);
    this->comm->SetDetectionProbability(ui.Settings().detection_probability);
    this->comm->SetPulsePeriod(ui.Settings().pulse_period);
    this->comm->SetSyncDivider(ui.Settings().sync_divider);
    this->comm->SetNumEmitters(ui.Settings().n_emitters);
    this->comm->SetSimSpeedupFactor(ui.Settings().sim_speedup_factor);
}

void fake_timetag_device_plugin::ShowInitDialog(QWidget* parent)
{
    std::cout << "You are using a FAKE time tage device, intended for debugging and testing purposes. "
                 "This device has three channels: chan 0: pulse channel. chan 1: detector 1. chan 2: detector 2." << std::endl;

#if defined(_WIN32) && !defined(NDEBUG)
    std::cout << "Currently the fake time tag device is OFF. Press the F9 key to start generating fake "
                 "timestamps. Press F10 to stop generating fake timestamps.";
#endif

    ui_fake_device_settings::fake_device_settings settings;
    settings.detection_probability = 0.3;
    settings.lifetime = 50e-9;
    settings.pulse_period = 1e-6;
    settings.time_unit = 1e-12;
    settings.sync_divider = 128;
    settings.n_emitters = 1;
    settings.sim_speedup_factor = 1.0;
    ui_fake_device_settings ui(parent, settings);

    ui.exec();

    delete this->comm;

    this->comm = new FakeCommunicator(ui.Settings().time_unit,
                                      1.0/ui.Settings().lifetime,
                                      ui.Settings().pulse_period,
                                      ui.Settings().detection_probability,
                                      ui.Settings().sync_divider,
                                      ui.Settings().chan1_delay,
                                      ui.Settings().chan2_delay,
                                      ui.Settings().n_emitters,
                                      ui.Settings().sim_speedup_factor);
}
