/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

#ifndef UI_FAKE_DEVICE_SETTINGS_H
#define UI_FAKE_DEVICE_SETTINGS_H

#include <QDialog>

namespace Ui {
class ui_fake_device_settings;
}

class ui_fake_device_settings : public QDialog
{
    Q_OBJECT

public:
    struct fake_device_settings
    {
      public:
        double time_unit;
        double lifetime;
        double pulse_period;
        double detection_probability;
        uint64_t sync_divider;
        int64_t chan1_delay;
        int64_t chan2_delay;
        int64_t n_emitters;
        double sim_speedup_factor;
    };

    fake_device_settings settings;

public:
    explicit ui_fake_device_settings(QWidget *parent, fake_device_settings initial_settings, bool can_change_timeunit=true);
    ~ui_fake_device_settings();

    fake_device_settings Settings() const;

private slots:
    void on_cnt_time_unit_valueChanged(double value);

    void on_cnt_decay_rate_valueChanged(double value);

    void on_cnt_pulse_period_valueChanged(double value);

    void on_cnt_det_prob_valueChanged(double value);

    void on_cnt_sync_divider_valueChanged(double value);

    void on_cnt_chan1_delay_valueChanged(double value);

    void on_cnt_chan2_delay_valueChanged(double value);

    void on_cnt_n_emitters_valueChanged(double value);

    void on_cnt_sim_speedup_valueChanged(double value);

private:
    Ui::ui_fake_device_settings *ui;
};

#endif // UI_FAKE_DEVICE_SETTINGS_H
