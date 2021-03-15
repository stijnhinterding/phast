/* Copyright (c) 2020 Stijn Hinterding, Utrecht University
 * This sofware is licensed under the MIT license (see the LICENSE file)	
*/

#ifndef CORRELATION_CALC_H
#define CORRELATION_CALC_H

#include <phast_gui/base/base_ui_calc.h>

class correlation_calc : public Base_ui_calc
{
private:
    std::vector<int64_t> xs;
    std::vector<int64_t> ys;

    chan_id leftchannel;
    chan_id rightchannel;

    timestamp t_min;
    timestamp t_max;

    double bin_low;
    double bin_high;

    uint64_t n_photons_left;
    uint64_t n_photons_right;

    bool correlate_FCS;

    int64_t max_num_bins;

    bool lin_FCS;
    int64_t num_lin_FCS_bins;

    void update_tmin(DetectionUpdate event);
    void update_tmax(DetectionUpdate event);

    void update_n_photons_leftchan(DetectionUpdate event);
    void update_n_photons_rightchan(DetectionUpdate event);

protected:
    virtual chan_data_map get_times_to_save(DetectionUpdate event) const;

    virtual DataView::data_update arrange_update(DetectionUpdate event, bool rebin_only=false);

public:
    correlation_calc(const std::vector<chan_id>& chans_to_correlate,
                   uint64_t bin_width,
                   opaque_ptr our_ref,
                   double bin_low,
                   double bin_high, bool display_only_initial_data, double time_unit, bool correlate_FCS, int64_t max_num_bins=50000000, int64_t num_lin_FCS_bins=1000);

    virtual void ClearData();
};
#endif // CORRELATION_CALC_H
