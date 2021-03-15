/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

#include "ui_fake_device_settings.h"
#include "ui_ui_fake_device_settings.h"

#include <qwt/qwt_counter.h>

ui_fake_device_settings::ui_fake_device_settings(QWidget *parent, fake_device_settings initial_settings, bool can_change_timeunit) :
    QDialog(parent),
    ui(new Ui::ui_fake_device_settings),
    settings(initial_settings)
{
    ui->setupUi(this);


    this->setWindowFlags(Qt::Window | Qt::WindowTitleHint | Qt::CustomizeWindowHint);

    this->ui->cnt_time_unit->setValue(this->settings.time_unit);
    this->ui->cnt_decay_rate->setValue(this->settings.lifetime);
    this->ui->cnt_pulse_period->setValue(this->settings.pulse_period);
    this->ui->cnt_det_prob->setValue(this->settings.detection_probability);
    this->ui->cnt_sync_divider->setValue(this->settings.sync_divider);
    this->ui->cnt_n_emitters->setValue(this->settings.n_emitters);
    this->ui->cnt_sim_speedup->setValue(this->settings.sim_speedup_factor);
    if (!can_change_timeunit) {
        this->ui->cnt_time_unit->setEnabled(false);
    }
}

ui_fake_device_settings::~ui_fake_device_settings()
{
    delete ui;
}

void ui_fake_device_settings::on_cnt_time_unit_valueChanged(double value)
{
    this->settings.time_unit = value;
}

ui_fake_device_settings::fake_device_settings ui_fake_device_settings::Settings() const
{
    return this->settings;
}

void ui_fake_device_settings::on_cnt_decay_rate_valueChanged(double value)
{
    this->settings.lifetime = value;
}

void ui_fake_device_settings::on_cnt_pulse_period_valueChanged(double value)
{
    this->settings.pulse_period = value;
}

void ui_fake_device_settings::on_cnt_det_prob_valueChanged(double value)
{
    this->settings.detection_probability = value;
}

void ui_fake_device_settings::on_cnt_sync_divider_valueChanged(double value)
{
    this->settings.sync_divider = value;
}

void ui_fake_device_settings::on_cnt_chan1_delay_valueChanged(double value)
{
    this->settings.chan1_delay = (int64_t)value;
}

void ui_fake_device_settings::on_cnt_chan2_delay_valueChanged(double value)
{
    this->settings.chan2_delay = (int64_t)value;
}

void ui_fake_device_settings::on_cnt_n_emitters_valueChanged(double value)
{
    this->settings.n_emitters = (int64_t)value;
}

void ui_fake_device_settings::on_cnt_sim_speedup_valueChanged(double value)
{
    this->settings.sim_speedup_factor = value;
}
